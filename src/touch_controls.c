#include "touch_controls.h"
#include "camera.h"
#include "sokol_app.h"
#include "util/sokol_gl.h"
#include "util/sokol_debugtext.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- layout constants ---------- */

/* Joysticks */
#define STICK_OUTER_RADIUS      70.0f   /* outer ring radius */
#define STICK_KNOB_RADIUS       24.0f   /* inner thumb knob */
#define STICK_DEADZONE          0.15f
#define STICK_MARGIN_X          120.0f  /* center of joystick from screen edge */
#define STICK_MARGIN_Y          160.0f  /* center of joystick from bottom */
#define TOUCH_LOOK_SENSITIVITY  0.005f

/* Face buttons (diamond / ABXY layout) */
#define FACE_BTN_RADIUS         28.0f
#define FACE_CENTER_MARGIN_X    100.0f  /* diamond center from right edge */
#define FACE_CENTER_MARGIN_Y    200.0f  /* diamond center from bottom edge */
#define FACE_SPREAD             52.0f   /* distance from diamond center to each button center */

/* Menu / start button */
#define MENU_BTN_W              70.0f
#define MENU_BTN_H              30.0f
#define MENU_BTN_RADIUS         35.0f   /* hit area radius */
#define MENU_TOP_MARGIN         30.0f   /* from top edge */

/* Hotbar */
#define HOTBAR_VISUAL_SLOTS     8
#define HOTBAR_SLOT_SIZE        40.0f
#define HOTBAR_PADDING          4.0f
#define HOTBAR_BOTTOM_MARGIN    12.0f

/* ---------- helpers ---------- */

static bool point_in_circle(float px, float py,
                            float cx, float cy, float r) {
    float dx = px - cx, dy = py - cy;
    return dx * dx + dy * dy <= r * r;
}

static bool point_in_rect(float px, float py,
                           float rx, float ry, float rw, float rh) {
    return px >= rx && px <= rx + rw && py >= ry && py <= ry + rh;
}

/* Compute button positions from current screen size */
static void update_layout(TouchControls* tc) {
    float sw = sapp_widthf();
    float sh = sapp_heightf();

    /* Face buttons — diamond centered at (sw - margin, sh - margin) */
    float dcx = sw - FACE_CENTER_MARGIN_X;
    float dcy = sh - FACE_CENTER_MARGIN_Y;

    /* Top = Jump (A), Bottom = Crouch (B), Left = Break (X), Right = Place (Y) */
    tc->btn_jump.x   = dcx;                tc->btn_jump.y   = dcy - FACE_SPREAD;
    tc->btn_crouch.x = dcx;                tc->btn_crouch.y = dcy + FACE_SPREAD;
    tc->btn_break.x  = dcx - FACE_SPREAD;  tc->btn_break.y  = dcy;
    tc->btn_place.x  = dcx + FACE_SPREAD;  tc->btn_place.y  = dcy;

    tc->btn_jump.radius   = FACE_BTN_RADIUS;
    tc->btn_crouch.radius = FACE_BTN_RADIUS;
    tc->btn_break.radius  = FACE_BTN_RADIUS;
    tc->btn_place.radius  = FACE_BTN_RADIUS;

    /* Menu button — top center */
    tc->btn_menu.x = sw * 0.5f;
    tc->btn_menu.y = MENU_TOP_MARGIN + MENU_BTN_H * 0.5f;
    tc->btn_menu.radius = MENU_BTN_RADIUS;
}

static int hotbar_hit_test(float px, float py, int hotbar_count) {
    float sw = sapp_widthf();
    float sh = sapp_heightf();
    float total_w = HOTBAR_VISUAL_SLOTS * HOTBAR_SLOT_SIZE
                  + (HOTBAR_VISUAL_SLOTS - 1) * HOTBAR_PADDING;
    float x0 = (sw - total_w) * 0.5f;
    float y0 = sh - HOTBAR_SLOT_SIZE - HOTBAR_BOTTOM_MARGIN;

    if (py < y0 || py > y0 + HOTBAR_SLOT_SIZE) return -1;

    for (int i = 0; i < hotbar_count; i++) {
        float sx = x0 + i * (HOTBAR_SLOT_SIZE + HOTBAR_PADDING);
        if (px >= sx && px < sx + HOTBAR_SLOT_SIZE) return i;
    }
    return -1;
}

static TouchButton* find_button_hit(TouchControls* tc, float px, float py) {
    float forgive = 1.2f;
    if (point_in_circle(px, py, tc->btn_jump.x,   tc->btn_jump.y,   tc->btn_jump.radius   * forgive)) return &tc->btn_jump;
    if (point_in_circle(px, py, tc->btn_crouch.x,  tc->btn_crouch.y, tc->btn_crouch.radius * forgive)) return &tc->btn_crouch;
    if (point_in_circle(px, py, tc->btn_break.x,   tc->btn_break.y,  tc->btn_break.radius  * forgive)) return &tc->btn_break;
    if (point_in_circle(px, py, tc->btn_place.x,   tc->btn_place.y,  tc->btn_place.radius  * forgive)) return &tc->btn_place;
    if (point_in_circle(px, py, tc->btn_menu.x,    tc->btn_menu.y,   tc->btn_menu.radius   * forgive)) return &tc->btn_menu;
    return NULL;
}

static void release_button_by_id(TouchControls* tc, uintptr_t tid) {
    TouchButton* btns[] = {
        &tc->btn_jump, &tc->btn_crouch, &tc->btn_break, &tc->btn_place, &tc->btn_menu
    };
    for (int i = 0; i < 5; i++) {
        if (btns[i]->pressed && btns[i]->touch_id == tid) {
            btns[i]->pressed = false;
            btns[i]->just_released = true;
            return;
        }
    }
}

/* ---------- public API ---------- */

void touch_init(TouchControls* tc) {
    memset(tc, 0, sizeof(*tc));
    tc->hotbar_touch_slot = -1;
    tc->move_stick.radius = STICK_OUTER_RADIUS;
    tc->look_stick.radius = STICK_OUTER_RADIUS;
}

bool touch_handle_event(TouchControls* tc, const sapp_event* ev, int hotbar_count) {
    if (ev->type != SAPP_EVENTTYPE_TOUCHES_BEGAN &&
        ev->type != SAPP_EVENTTYPE_TOUCHES_MOVED &&
        ev->type != SAPP_EVENTTYPE_TOUCHES_ENDED &&
        ev->type != SAPP_EVENTTYPE_TOUCHES_CANCELLED) {
        return false;
    }

    tc->is_touch_device = true;
    tc->hotbar_touch_slot = -1;
    update_layout(tc);

    float sw = sapp_widthf();

    for (int t = 0; t < ev->num_touches; t++) {
        if (!ev->touches[t].changed) continue;

        uintptr_t tid = ev->touches[t].identifier;
        float px = ev->touches[t].pos_x;
        float py = ev->touches[t].pos_y;

        if (ev->type == SAPP_EVENTTYPE_TOUCHES_BEGAN) {
            /* 1. Hotbar */
            int slot = hotbar_hit_test(px, py, hotbar_count);
            if (slot >= 0) {
                tc->hotbar_touch_slot = slot;
                continue;
            }

            /* 2. Action buttons (face buttons + menu) */
            TouchButton* btn = find_button_hit(tc, px, py);
            if (btn) {
                btn->pressed  = true;
                btn->touch_id = tid;
                continue;
            }

            /* 3. Left half → move joystick (floating) */
            if (px < sw * 0.4f) {
                tc->move_stick.active    = true;
                tc->move_stick.touch_id  = tid;
                tc->move_stick.origin_x  = px;
                tc->move_stick.origin_y  = py;
                tc->move_stick.current_x = px;
                tc->move_stick.current_y = py;
                tc->move_stick.dx = 0;
                tc->move_stick.dy = 0;
                continue;
            }

            /* 4. Right half (middle area) → look joystick (floating) */
            if (!tc->look_stick.active) {
                tc->look_stick.active    = true;
                tc->look_stick.touch_id  = tid;
                tc->look_stick.origin_x  = px;
                tc->look_stick.origin_y  = py;
                tc->look_stick.current_x = px;
                tc->look_stick.current_y = py;
                tc->look_stick.dx = 0;
                tc->look_stick.dy = 0;
                continue;
            }
        }
        else if (ev->type == SAPP_EVENTTYPE_TOUCHES_MOVED) {
            /* Move joystick */
            if (tc->move_stick.active && tc->move_stick.touch_id == tid) {
                tc->move_stick.current_x = px;
                tc->move_stick.current_y = py;
                float raw_dx =  (px - tc->move_stick.origin_x) / tc->move_stick.radius;
                float raw_dy = -(py - tc->move_stick.origin_y) / tc->move_stick.radius;
                float len = sqrtf(raw_dx * raw_dx + raw_dy * raw_dy);
                if (len > 1.0f) { raw_dx /= len; raw_dy /= len; }
                tc->move_stick.dx = raw_dx;
                tc->move_stick.dy = raw_dy;
                continue;
            }

            /* Look joystick — accumulate delta */
            if (tc->look_stick.active && tc->look_stick.touch_id == tid) {
                tc->look_stick.dx += px - tc->look_stick.current_x;
                tc->look_stick.dy += py - tc->look_stick.current_y;
                tc->look_stick.current_x = px;
                tc->look_stick.current_y = py;
                continue;
            }
        }
        else {
            /* TOUCHES_ENDED / CANCELLED */
            if (tc->move_stick.active && tc->move_stick.touch_id == tid) {
                tc->move_stick.active = false;
                tc->move_stick.dx = 0;
                tc->move_stick.dy = 0;
                continue;
            }
            if (tc->look_stick.active && tc->look_stick.touch_id == tid) {
                tc->look_stick.active = false;
                tc->look_stick.dx = 0;
                tc->look_stick.dy = 0;
                continue;
            }
            release_button_by_id(tc, tid);
        }
    }

    return true;
}

void touch_apply_to_camera(TouchControls* tc, Camera* cam, float dt) {
    if (!tc->is_touch_device) return;
    (void)dt;

    /* Clear movement keys — touch is the sole input source each frame */
    cam->key_w = cam->key_s = cam->key_a = cam->key_d = false;
    cam->key_space = cam->key_ctrl = false;

    /* Move stick → movement keys */
    if (tc->move_stick.dy >  STICK_DEADZONE) cam->key_w = true;
    if (tc->move_stick.dy < -STICK_DEADZONE) cam->key_s = true;
    if (tc->move_stick.dx >  STICK_DEADZONE) cam->key_d = true;
    if (tc->move_stick.dx < -STICK_DEADZONE) cam->key_a = true;

    /* Look joystick delta → yaw / pitch */
    cam->yaw   -= tc->look_stick.dx * TOUCH_LOOK_SENSITIVITY;
    cam->pitch -= tc->look_stick.dy * TOUCH_LOOK_SENSITIVITY;
    float limit = (float)M_PI / 2.0f - 0.01f;
    if (cam->pitch >  limit) cam->pitch =  limit;
    if (cam->pitch < -limit) cam->pitch = -limit;
    tc->look_stick.dx = 0;
    tc->look_stick.dy = 0;

    /* Buttons → keys */
    if (tc->btn_jump.pressed)   cam->key_space = true;
    if (tc->btn_crouch.pressed) cam->key_ctrl  = true;

    cam->mouse_locked = true;
}

/* ---------- rendering helpers ---------- */

static void draw_filled_circle(float cx, float cy, float r,
                               float cr, float cg, float cb, float ca,
                               int segs) {
    sgl_begin_triangles();
    sgl_c4f(cr, cg, cb, ca);
    for (int i = 0; i < segs; i++) {
        float a0 = 2.0f * (float)M_PI * (float)i / (float)segs;
        float a1 = 2.0f * (float)M_PI * (float)(i + 1) / (float)segs;
        sgl_v2f(cx, cy);
        sgl_v2f(cx + cosf(a0) * r, cy + sinf(a0) * r);
        sgl_v2f(cx + cosf(a1) * r, cy + sinf(a1) * r);
    }
    sgl_end();
}

static void draw_ring(float cx, float cy, float r,
                      float cr, float cg, float cb, float ca,
                      int segs, float thickness) {
    float inner = r - thickness * 0.5f;
    float outer = r + thickness * 0.5f;
    sgl_begin_triangle_strip();
    sgl_c4f(cr, cg, cb, ca);
    for (int i = 0; i <= segs; i++) {
        float a = 2.0f * (float)M_PI * (float)i / (float)segs;
        float cs = cosf(a), sn = sinf(a);
        sgl_v2f(cx + cs * inner, cy + sn * inner);
        sgl_v2f(cx + cs * outer, cy + sn * outer);
    }
    sgl_end();
}

static void draw_rounded_rect(float cx, float cy, float w, float h,
                               float cr, float cg, float cb, float ca) {
    float hw = w * 0.5f, hh = h * 0.5f;
    sgl_begin_triangles();
    sgl_c4f(cr, cg, cb, ca);
    /* Two triangles for the rect */
    sgl_v2f(cx - hw, cy - hh);
    sgl_v2f(cx + hw, cy - hh);
    sgl_v2f(cx + hw, cy + hh);
    sgl_v2f(cx - hw, cy - hh);
    sgl_v2f(cx + hw, cy + hh);
    sgl_v2f(cx - hw, cy + hh);
    sgl_end();
}

static void setup_sgl_ortho(void) {
    sgl_defaults();
    sgl_matrix_mode_projection();
    sgl_load_identity();
    sgl_ortho(0.0f, sapp_widthf(), sapp_heightf(), 0.0f, -1.0f, 1.0f);
    sgl_matrix_mode_modelview();
    sgl_load_identity();
}

/* ---------- public render functions ---------- */

void touch_render_crosshair(void) {
    float cx = sapp_widthf() * 0.5f;
    float cy = sapp_heightf() * 0.5f;

    setup_sgl_ortho();
    draw_filled_circle(cx, cy, 3.0f, 1.0f, 1.0f, 1.0f, 1.0f, 12);
    sgl_draw();
}

void touch_render(const TouchControls* tc, int hotbar_selected, int hotbar_count) {
    if (!tc->is_touch_device) return;
    (void)hotbar_selected;
    (void)hotbar_count;

    float sw = sapp_widthf();
    float sh = sapp_heightf();

    setup_sgl_ortho();

    /* ---- Left joystick base (always visible as ghost) ---- */
    {
        float base_x = STICK_MARGIN_X;
        float base_y = sh - STICK_MARGIN_Y;
        /* Ghost ring — always shown so player knows where to touch */
        draw_ring(base_x, base_y, STICK_OUTER_RADIUS,
                  1.0f, 1.0f, 1.0f, 0.15f, 24, 2.0f);
    }

    /* ---- Left joystick active ---- */
    if (tc->move_stick.active) {
        float ox = tc->move_stick.origin_x;
        float oy = tc->move_stick.origin_y;
        float r  = tc->move_stick.radius;

        draw_ring(ox, oy, r, 1.0f, 1.0f, 1.0f, 0.3f, 24, 2.0f);
        float knob_x = ox + tc->move_stick.dx * r;
        float knob_y = oy - tc->move_stick.dy * r;
        draw_filled_circle(knob_x, knob_y, STICK_KNOB_RADIUS,
                           1.0f, 1.0f, 1.0f, 0.6f, 16);
    }

    /* ---- Right joystick base (always visible as ghost) ---- */
    {
        float base_x = sw - STICK_MARGIN_X;
        float base_y = sh - STICK_MARGIN_Y;
        draw_ring(base_x, base_y, STICK_OUTER_RADIUS,
                  1.0f, 1.0f, 1.0f, 0.15f, 24, 2.0f);
    }

    /* ---- Right joystick active ---- */
    if (tc->look_stick.active) {
        float ox = tc->look_stick.origin_x;
        float oy = tc->look_stick.origin_y;
        float r  = tc->look_stick.radius;

        draw_ring(ox, oy, r, 1.0f, 1.0f, 1.0f, 0.3f, 24, 2.0f);
        float knob_x = tc->look_stick.current_x;
        float knob_y = tc->look_stick.current_y;
        /* Clamp knob to circle */
        float dx = knob_x - ox, dy = knob_y - oy;
        float dist = sqrtf(dx * dx + dy * dy);
        if (dist > r) {
            knob_x = ox + dx / dist * r;
            knob_y = oy + dy / dist * r;
        }
        draw_filled_circle(knob_x, knob_y, STICK_KNOB_RADIUS,
                           1.0f, 1.0f, 1.0f, 0.6f, 16);
    }

    /* ---- Face buttons (diamond layout) ---- */
    struct {
        const TouchButton* btn;
        float cr, cg, cb;
        char label;
    } faces[] = {
        { &tc->btn_jump,   0.2f, 0.8f, 0.2f, 'A' },  /* Top — green */
        { &tc->btn_crouch, 0.3f, 0.5f, 0.9f, 'B' },  /* Bottom — blue */
        { &tc->btn_break,  0.9f, 0.2f, 0.2f, 'X' },  /* Left — red */
        { &tc->btn_place,  0.9f, 0.8f, 0.2f, 'Y' },  /* Right — yellow */
    };
    for (int i = 0; i < 4; i++) {
        float bcx = faces[i].btn->x;
        float bcy = faces[i].btn->y;
        float br  = faces[i].btn->radius;
        float a   = faces[i].btn->pressed ? 0.85f : 0.45f;
        /* Dark background circle */
        draw_filled_circle(bcx, bcy, br, 0.0f, 0.0f, 0.0f, 0.3f, 16);
        /* Colored ring */
        draw_ring(bcx, bcy, br, faces[i].cr, faces[i].cg, faces[i].cb, a, 16, 3.0f);
        /* Colored fill when pressed */
        if (faces[i].btn->pressed) {
            draw_filled_circle(bcx, bcy, br - 3.0f,
                               faces[i].cr, faces[i].cg, faces[i].cb, 0.3f, 16);
        }
    }

    /* ---- Menu / Start button (top center, pill shape) ---- */
    {
        float mx = tc->btn_menu.x;
        float my = tc->btn_menu.y;
        float a  = tc->btn_menu.pressed ? 0.7f : 0.35f;
        draw_rounded_rect(mx, my, MENU_BTN_W, MENU_BTN_H,
                          0.4f, 0.4f, 0.4f, a);
        /* Three horizontal lines (hamburger icon) */
        float lw = 16.0f, lh = 2.0f, gap = 5.0f;
        for (int i = -1; i <= 1; i++) {
            draw_rounded_rect(mx, my + (float)i * gap, lw, lh,
                              1.0f, 1.0f, 1.0f, a + 0.2f);
        }
    }

    sgl_draw();

    /* ---- Button labels via sdtx ---- */
    {
        float canvas_w = sapp_widthf()  * 0.5f;
        float canvas_h = sapp_heightf() * 0.5f;
        sdtx_canvas(canvas_w, canvas_h);
        sdtx_font(0);

        float sx = canvas_w / sapp_widthf()  / 8.0f;
        float sy = canvas_h / sapp_heightf() / 8.0f;

        for (int i = 0; i < 4; i++) {
            float bcx = faces[i].btn->x;
            float bcy = faces[i].btn->y;
            sdtx_origin(bcx * sx - 0.5f, bcy * sy - 0.5f);
            sdtx_home();
            sdtx_color3f(1.0f, 1.0f, 1.0f);
            sdtx_putc(faces[i].label);
        }
        sdtx_draw();
    }
}

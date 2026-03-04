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

/* ---------- layout constants (hotbar must match render.c) ---------- */
#define HOTBAR_VISUAL_SLOTS     8
#define HOTBAR_SLOT_SIZE        40.0f
#define HOTBAR_PADDING          4.0f
#define HOTBAR_BOTTOM_MARGIN    12.0f

/* buttons */
#define BTN_RADIUS              30.0f
#define BTN_SPACING             70.0f   /* diameter + 10 px gap */
#define BTN_RIGHT_MARGIN        50.0f
#define BTN_BASE_Y_OFFSET       200.0f  /* from bottom edge to center of Place btn */

/* joystick */
#define STICK_RADIUS            60.0f
#define STICK_DEADZONE          0.3f
#define TOUCH_LOOK_SENSITIVITY  0.005f

/* ---------- helpers ---------- */

static bool point_in_circle(float px, float py,
                            float cx, float cy, float r) {
    float dx = px - cx, dy = py - cy;
    return dx * dx + dy * dy <= r * r;
}

static void update_button_layout(TouchControls* tc) {
    float sw = sapp_widthf();
    float sh = sapp_heightf();
    float r  = BTN_RADIUS;
    float d  = r * 2.0f;
    float cx = sw - BTN_RIGHT_MARGIN;
    float base_y = sh - BTN_BASE_Y_OFFSET;

    /* From bottom up: Place, Break, Crouch, Jump */
    tc->btn_place.x  = cx - r;  tc->btn_place.y  = base_y - r;
    tc->btn_place.w  = d;       tc->btn_place.h  = d;

    tc->btn_break.x  = cx - r;  tc->btn_break.y  = base_y - BTN_SPACING - r;
    tc->btn_break.w  = d;       tc->btn_break.h  = d;

    tc->btn_crouch.x = cx - r;  tc->btn_crouch.y = base_y - 2 * BTN_SPACING - r;
    tc->btn_crouch.w = d;       tc->btn_crouch.h = d;

    tc->btn_jump.x   = cx - r;  tc->btn_jump.y   = base_y - 3 * BTN_SPACING - r;
    tc->btn_jump.w   = d;       tc->btn_jump.h   = d;
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
    float r = BTN_RADIUS;
    float hit_r = r * 1.2f;  /* slightly forgiving hit area */
    if (point_in_circle(px, py, tc->btn_place.x  + r, tc->btn_place.y  + r, hit_r)) return &tc->btn_place;
    if (point_in_circle(px, py, tc->btn_break.x  + r, tc->btn_break.y  + r, hit_r)) return &tc->btn_break;
    if (point_in_circle(px, py, tc->btn_crouch.x + r, tc->btn_crouch.y + r, hit_r)) return &tc->btn_crouch;
    if (point_in_circle(px, py, tc->btn_jump.x   + r, tc->btn_jump.y   + r, hit_r)) return &tc->btn_jump;
    return NULL;
}

static void release_button_by_id(TouchControls* tc, uintptr_t tid) {
    TouchButton* btns[] = {
        &tc->btn_place, &tc->btn_break, &tc->btn_crouch, &tc->btn_jump
    };
    for (int i = 0; i < 4; i++) {
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
    tc->move_stick.radius = STICK_RADIUS;
    tc->look_stick.radius = STICK_RADIUS;
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
    update_button_layout(tc);

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

            /* 2. Action buttons */
            TouchButton* btn = find_button_hit(tc, px, py);
            if (btn) {
                btn->pressed  = true;
                btn->touch_id = tid;
                continue;
            }

            /* 3. Left half → move joystick */
            if (px < sw * 0.5f) {
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

            /* 4. Right half → look zone */
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
                /* clamp to unit circle */
                float len = sqrtf(raw_dx * raw_dx + raw_dy * raw_dy);
                if (len > 1.0f) { raw_dx /= len; raw_dy /= len; }
                tc->move_stick.dx = raw_dx;
                tc->move_stick.dy = raw_dy;  /* positive = up/forward */
                continue;
            }

            /* Look zone — accumulate delta */
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

    /* Move stick → movement keys (dy > 0 = forward = up on screen) */
    if (tc->move_stick.dy >  STICK_DEADZONE) cam->key_w = true;
    if (tc->move_stick.dy < -STICK_DEADZONE) cam->key_s = true;
    if (tc->move_stick.dx >  STICK_DEADZONE) cam->key_d = true;
    if (tc->move_stick.dx < -STICK_DEADZONE) cam->key_a = true;

    /* Look zone delta → yaw / pitch */
    cam->yaw   -= tc->look_stick.dx * TOUCH_LOOK_SENSITIVITY;
    cam->pitch  -= tc->look_stick.dy * TOUCH_LOOK_SENSITIVITY;
    float limit = (float)M_PI / 2.0f - 0.01f;
    if (cam->pitch >  limit) cam->pitch =  limit;
    if (cam->pitch < -limit) cam->pitch = -limit;
    /* clear accumulated delta */
    tc->look_stick.dx = 0;
    tc->look_stick.dy = 0;

    /* Buttons → keys */
    if (tc->btn_jump.pressed)   cam->key_space = true;
    if (tc->btn_crouch.pressed) cam->key_ctrl  = true;

    /* Auto-lock mouse so the game logic treats us as "playing" */
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

    setup_sgl_ortho();

    /* ---- Move joystick (only when active) ---- */
    if (tc->move_stick.active) {
        float ox = tc->move_stick.origin_x;
        float oy = tc->move_stick.origin_y;
        float r  = tc->move_stick.radius;

        /* outer ring */
        draw_ring(ox, oy, r, 1.0f, 1.0f, 1.0f, 0.3f, 24, 2.0f);

        /* inner knob — note dy is flipped (positive = screen-up) */
        float knob_x = ox + tc->move_stick.dx * r;
        float knob_y = oy - tc->move_stick.dy * r;
        draw_filled_circle(knob_x, knob_y, 20.0f, 1.0f, 1.0f, 1.0f, 0.6f, 16);
    }

    /* ---- Action buttons (always visible) ---- */
    float r = BTN_RADIUS;

    /* Jump — green */
    {
        float cx = tc->btn_jump.x + r, cy = tc->btn_jump.y + r;
        float a  = tc->btn_jump.pressed ? 0.8f : 0.5f;
        draw_filled_circle(cx, cy, r, 0.2f, 0.8f, 0.2f, a, 16);
    }
    /* Crouch — blue */
    {
        float cx = tc->btn_crouch.x + r, cy = tc->btn_crouch.y + r;
        float a  = tc->btn_crouch.pressed ? 0.8f : 0.5f;
        draw_filled_circle(cx, cy, r, 0.2f, 0.4f, 0.9f, a, 16);
    }
    /* Break — red */
    {
        float cx = tc->btn_break.x + r, cy = tc->btn_break.y + r;
        float a  = tc->btn_break.pressed ? 0.8f : 0.5f;
        draw_filled_circle(cx, cy, r, 0.9f, 0.2f, 0.2f, a, 16);
    }
    /* Place — yellow */
    {
        float cx = tc->btn_place.x + r, cy = tc->btn_place.y + r;
        float a  = tc->btn_place.pressed ? 0.8f : 0.5f;
        draw_filled_circle(cx, cy, r, 0.9f, 0.8f, 0.2f, a, 16);
    }

    sgl_draw();

    /* ---- Button labels (single char each, via sdtx) ---- */
    {
        float canvas_w = sapp_widthf()  * 0.5f;
        float canvas_h = sapp_heightf() * 0.5f;
        sdtx_canvas(canvas_w, canvas_h);
        sdtx_font(0);

        /* scale: screen-px → sdtx char-cell position */
        float sx = canvas_w / sapp_widthf()  / 8.0f;   /* 0.5 / 8 */
        float sy = canvas_h / sapp_heightf() / 8.0f;

        struct { const TouchButton* btn; char ch; } labels[] = {
            { &tc->btn_jump,   'J' },
            { &tc->btn_crouch, 'C' },
            { &tc->btn_break,  'X' },
            { &tc->btn_place,  'P' },
        };

        for (int i = 0; i < 4; i++) {
            float bcx = labels[i].btn->x + r;  /* screen center */
            float bcy = labels[i].btn->y + r;
            sdtx_origin(bcx * sx - 0.5f, bcy * sy - 0.5f);
            sdtx_home();
            sdtx_color3f(1.0f, 1.0f, 1.0f);
            sdtx_putc(labels[i].ch);
        }
        sdtx_draw();
    }
}

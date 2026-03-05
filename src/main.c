#include "sokol_app.h"
#include "sokol_gfx.h"
#include "sokol_glue.h"
#include "sokol_log.h"
#include "sokol_time.h"
#include "util/sokol_debugtext.h"
#include "util/sokol_gl.h"

#define HANDMADE_MATH_IMPLEMENTATION
#include "HandmadeMath.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#pragma warning(disable: 4996)  // freopen deprecation
#else
#include <pthread.h>
#include <sys/stat.h>
#endif

#include "camera.h"
#include "render.h"
#include "planet.h"
#include "screenshot.h"
#include "log_config.h"
#include "crash_handler.h"
#include "world.h"
#include "lobby.h"
#include "touch_controls.h"

typedef enum {
    STATE_MENU,
    STATE_WORLD_SELECT,
    STATE_SPAWN_SELECT,
    STATE_MULTI_MENU,
    STATE_MULTI_HOST,
    STATE_MULTI_JOIN,
    STATE_LOADING,
    STATE_PLAYING,
} GameState;

// ---- Background thread for planet_init ----
typedef struct {
    Planet* planet;
    int subdivision;
    volatile int completed;  // 0 = running, 1 = done
} PlanetInitTask;

#ifdef _WIN32
static DWORD WINAPI planet_init_thread(LPVOID param) {
    PlanetInitTask* task = (PlanetInitTask*)param;
    planet_init(task->planet, task->subdivision);
    task->completed = 1;
    return 0;
}
#else
static void* planet_init_thread(void* param) {
    PlanetInitTask* task = (PlanetInitTask*)param;
    planet_init(task->planet, task->subdivision);
    task->completed = 1;
    return NULL;
}
#endif

static struct {
    Camera camera;
    Planet planet;
    Renderer renderer;
    uint64_t last_time;
    GameState state;
    int loading_phase;  // 0 = spawn thread, 1 = wait for planet_init, 2 = render_init, 3 = LOD loading
    bool screenshot_requested;
    uint64_t load_start_time;   // Benchmark: time when loading started
    int load_frame_count;       // Benchmark: frames spent in loading phase 3

    // Block selector hotbar
    uint8_t selected_block_type;
    int hotbar_slot;

    // World system
    WorldList world_list;
    int active_world_idx;       // index into world_list.worlds[]
    char active_edits_dir[256]; // edits dir for active world
    bool is_multiplayer;
    bool is_host;
    int spawn_type;             // 0 = Summit (origin), 1 = Twilight

    // Touch controls
    TouchControls touch;

    // Lobby
    LobbyRequest lobby_req;

    // Background loading
    PlanetInitTask init_task;
#ifdef _WIN32
    HANDLE init_thread;
#else
    pthread_t init_thread;
#endif
} app;

bool log_verbose = false;

// Menu state
static int menu_selected = 0;
static int menu_hover = -1;
static int menu_pressed = -1;
static float menu_mouse_x = 0;
static float menu_mouse_y = 0;

// World select state
static int ws_hover = -1;
static int ws_pressed = -1;

// Spawn select state
static int sp_hover = -1;
static int sp_pressed = -1;
static int sp_selected = 0;

// Multiplayer menu state
static int multi_hover = -1;
static int multi_pressed = -1;

// Room code input state (STATE_MULTI_JOIN)
static char room_code[8] = "";
static int room_code_len = 0;
static int room_code_cursor_blink = 0;

// Block interaction: use hex terrain selection from renderer
static void interact_break(void) {
    if (!app.renderer.hex_placement.valid) return;
    HexHitResult* hit = &app.renderer.hex_placement;

    // Check if breaking a torch — need to remove visual instance
    uint8_t broken_type = hex_terrain_get_voxel(&app.renderer.hex_terrain,
        hit->gcol, hit->grow, hit->layer);

    if (hex_terrain_break(&app.renderer.hex_terrain, hit)) {
        if (broken_type == VOXEL_TORCH) {
            torch_remove(&app.renderer.torch_system, hit->gcol, hit->grow, hit->layer);
        }
    }
}

// Hotbar: placeable block types in order
static const uint8_t hotbar_types[] = {
    VOXEL_STONE, VOXEL_DIRT, VOXEL_GRASS, VOXEL_SAND,
    VOXEL_WATER, VOXEL_ICE, VOXEL_TORCH
};
#define HOTBAR_COUNT 7

static void interact_place(void) {
    if (!app.renderer.hex_placement.valid) return;
    HexHitResult* hit = &app.renderer.hex_placement;

    if (hex_terrain_place(&app.renderer.hex_terrain, hit, app.selected_block_type)) {
        if (app.selected_block_type == VOXEL_TORCH) {
            // Compute world position for the torch model
            float wx, wy, wz, ux, uy, uz;
            hex_terrain_hex_to_world(&app.renderer.hex_terrain,
                hit->place_gcol, hit->place_grow, hit->place_layer,
                &wx, &wy, &wz, &ux, &uy, &uz);

            int cx = (int)floorf((float)hit->place_gcol / HEX_CHUNK_SIZE);
            int cz = (int)floorf((float)hit->place_grow / HEX_CHUNK_SIZE);
            torch_add(&app.renderer.torch_system, wx, wy, wz, ux, uy, uz,
                      cx, cz, hit->place_gcol, hit->place_grow, hit->place_layer);
        }
    }
}

// ---- Helper to get the active world's player.dat path ----
static void get_player_path(char* out, int size) {
    if (app.active_world_idx >= 0 && app.active_world_idx < app.world_list.count) {
        world_get_player_path(&app.world_list.worlds[app.active_world_idx], out, size);
    } else {
        snprintf(out, size, "cache/player.dat");
    }
}

static void init(void) {
    // Allocate a console window for logging (Windows GUI apps have no console by default)
    #ifdef _WIN32
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    #endif

    crash_handler_install();

    // macOS: disable mouse event coalescing to prevent input lag
    #ifdef __APPLE__
    extern void platform_macos_init(void);
    platform_macos_init();
    #endif

    printf("[GAME] Hex Planets starting...\n");
    fflush(stdout);

    printf("[GAME] init: sg_setup...\n"); fflush(stdout);
    sg_setup(&(sg_desc){
        .environment = sglue_environment(),
        .logger.func = crash_log_func,
        .buffer_pool_size = 16384,  // LOD tree needs thousands of vertex buffers
    });

    sgl_setup(&(sgl_desc_t){ .logger.func = crash_log_func });

    printf("[GAME] init: stm_setup...\n"); fflush(stdout);
    stm_setup();

    printf("[GAME] init: sdtx_setup...\n"); fflush(stdout);
    // Setup debug text for menu
    sdtx_setup(&(sdtx_desc_t){
        .fonts[0] = sdtx_font_cpc(),
        .logger.func = crash_log_func,
    });

    printf("[GAME] init: camera_init...\n"); fflush(stdout);
    camera_init(&app.camera);

    // World system: load world list, migrate legacy saves
    app.active_world_idx = -1;
    world_migrate_legacy(&app.world_list);
    if (!world_list_load(&app.world_list)) {
        app.world_list.count = 0;
        app.world_list.active_index = -1;
    }

    // Lobby system
    lobby_init("https://hex-planets.vercel.app");

    // Touch controls
    touch_init(&app.touch);
    app.renderer.touch = &app.touch;

    app.state = STATE_MENU;
    app.last_time = stm_now();
    printf("[GAME] init: done, entering menu.\n"); fflush(stdout);
}

static void player_save(void);
static bool player_load(void);
static int autosave_counter;

static void start_loading(void) {
    app.state = STATE_LOADING;
    app.loading_phase = 0;
    app.selected_block_type = VOXEL_STONE;
    app.hotbar_slot = 0;

    // Set up edits dir from active world
    if (app.active_world_idx >= 0 && app.active_world_idx < app.world_list.count) {
        WorldMeta* w = &app.world_list.worlds[app.active_world_idx];
        world_get_edits_dir(w, app.active_edits_dir, sizeof(app.active_edits_dir));
        w->last_played = (uint64_t)time(NULL);
        app.world_list.active_index = app.active_world_idx;
        world_list_save(&app.world_list);
    } else {
        snprintf(app.active_edits_dir, sizeof(app.active_edits_dir), "cache/edits");
    }
    printf("[GAME] Loading with edits dir: %s\n", app.active_edits_dir);
    fflush(stdout);
}

static void draw_loading_screen(const char* status, int patches_ready, int patches_total, int verts) {
    sg_begin_pass(&(sg_pass){
        .action = {
            .colors[0] = {
                .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f}
            }
        },
        .swapchain = sglue_swapchain()
    });

    sdtx_canvas(sapp_widthf() * 0.5f, sapp_heightf() * 0.5f);
    sdtx_origin(3.0f, 5.0f);
    sdtx_font(0);

    sdtx_color3f(0.8f, 0.9f, 1.0f);
    sdtx_puts("=== HEX PLANETS ===\n\n");

    sdtx_color3f(1.0f, 1.0f, 0.5f);
    sdtx_printf("%s\n\n", status);

    if (patches_total > 0) {
        sdtx_color3f(0.6f, 0.8f, 0.6f);
        sdtx_printf("LOD patches: %d / %d\n", patches_ready, patches_total);
        sdtx_printf("Vertices:    %d\n", verts);

        // Progress bar
        int bar_width = 30;
        int filled = patches_total > 0 ? (patches_ready * bar_width / patches_total) : 0;
        sdtx_color3f(0.5f, 0.5f, 0.5f);
        sdtx_putc('[');
        for (int i = 0; i < bar_width; i++) {
            if (i < filled) {
                sdtx_color3f(0.3f, 1.0f, 0.3f);
                sdtx_putc('#');
            } else {
                sdtx_color3f(0.3f, 0.3f, 0.3f);
                sdtx_putc('.');
            }
        }
        sdtx_color3f(0.5f, 0.5f, 0.5f);
        sdtx_putc(']');
        sdtx_putc('\n');
    }

    sdtx_draw();
    sg_end_pass();
    sg_commit();
}

// ---- Reusable button drawing ----
// Draws a list of buttons with hover/press state, returns button dimensions for text overlay
static void draw_buttons(int count, float btn_w, float btn_h, float btn_gap,
                          float margin_x, float btn_y0, float win_w, float win_h,
                          int hover, int pressed_btn) {
    sgl_defaults();
    sgl_matrix_mode_projection();
    sgl_load_identity();
    sgl_ortho(0.0f, win_w, win_h, 0.0f, -1.0f, 1.0f);
    sgl_matrix_mode_modelview();
    sgl_load_identity();

    for (int i = 0; i < count; i++) {
        bool is_pressed = (pressed_btn == i && hover == i);
        bool is_hovered = (hover == i);
        float x0 = margin_x, y0 = btn_y0 + i * (btn_h + btn_gap);
        float x1 = x0 + btn_w, y1 = y0 + btn_h;

        float fr, fg, fb, fa;
        if (is_pressed)      { fr = 0.25f; fg = 0.25f; fb = 0.3f; fa = 0.9f; }
        else if (is_hovered) { fr = 0.12f; fg = 0.14f; fb = 0.2f; fa = 0.8f; }
        else                 { fr = 0.08f; fg = 0.08f; fb = 0.12f; fa = 0.6f; }
        sgl_begin_quads();
        sgl_c4f(fr, fg, fb, fa);
        sgl_v2f(x0, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y1); sgl_v2f(x0, y1);
        sgl_end();

        float br, bg_, bb;
        if (is_pressed)      { br = 1.0f; bg_ = 1.0f; bb = 1.0f; }
        else if (is_hovered) { br = 1.0f; bg_ = 1.0f; bb = 0.3f; }
        else                 { br = 0.3f; bg_ = 0.3f; bb = 0.4f; }
        float bw = 2.0f;
        sgl_begin_quads();
        sgl_c3f(br, bg_, bb);
        sgl_v2f(x0, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y0 + bw); sgl_v2f(x0, y0 + bw);
        sgl_v2f(x0, y1 - bw); sgl_v2f(x1, y1 - bw); sgl_v2f(x1, y1); sgl_v2f(x0, y1);
        sgl_v2f(x0, y0); sgl_v2f(x0 + bw, y0); sgl_v2f(x0 + bw, y1); sgl_v2f(x0, y1);
        sgl_v2f(x1 - bw, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y1); sgl_v2f(x1 - bw, y1);
        sgl_end();
    }
    sgl_draw();
}

// Hit-test helper: returns button index or -1
static int hit_test_buttons(int count, float btn_w, float btn_h, float btn_gap,
                             float margin_x, float btn_y0, float mx, float my) {
    for (int i = 0; i < count; i++) {
        float by = btn_y0 + i * (btn_h + btn_gap);
        if (mx >= margin_x && mx < margin_x + btn_w &&
            my >= by && my < by + btn_h) {
            return i;
        }
    }
    return -1;
}

static int frame_count = 0;

static void frame(void) {
    frame_count++;

    uint64_t now = stm_now();
    float dt = (float)stm_sec(stm_diff(now, app.last_time));
    app.last_time = now;
    if (dt > 0.1f) dt = 0.1f;

    // Shared layout values for all menu screens
    float win_w = sapp_widthf();
    float win_h = sapp_heightf();
    float canvas_w = win_w * 0.5f;
    float canvas_h = win_h * 0.5f;
    float char_px = win_w / (canvas_w / 8.0f);
    float margin_x = 3.0f * char_px;
    float title_y = 3.0f * char_px;

    // ---- STATE_MENU: Main menu ----
    if (app.state == STATE_MENU) {
        if (frame_count % 120 == 0) {
            printf("[MENU] frame %d alive\n", frame_count);
            fflush(stdout);
        }
        sg_begin_pass(&(sg_pass){
            .action = { .colors[0] = { .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f} } },
            .swapchain = sglue_swapchain()
        });

        #define MAIN_MENU_BTN_COUNT 2
        const char* btn_labels[MAIN_MENU_BTN_COUNT] = {
            "[1] Single Player",
            "[2] Multiplayer",
        };
        int max_label_len = 0;
        for (int i = 0; i < MAIN_MENU_BTN_COUNT; i++) {
            int len = (int)strlen(btn_labels[i]);
            if (len > max_label_len) max_label_len = len;
        }
        float btn_w = (max_label_len + 4) * char_px;
        float btn_h = 3.0f * char_px;
        float btn_gap = 10.0f;
        float btn_y0 = title_y + 2.0f * char_px;

        menu_hover = hit_test_buttons(MAIN_MENU_BTN_COUNT, btn_w, btn_h, btn_gap,
                                       margin_x, btn_y0, menu_mouse_x, menu_mouse_y);
        draw_buttons(MAIN_MENU_BTN_COUNT, btn_w, btn_h, btn_gap, margin_x, btn_y0,
                     win_w, win_h, menu_hover, menu_pressed);

        sdtx_canvas(canvas_w, canvas_h);
        sdtx_origin(0.0f, 0.0f);
        sdtx_font(0);
        sdtx_pos(margin_x / char_px, title_y / char_px);
        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("=== HEX PLANETS ===");

        for (int i = 0; i < MAIN_MENU_BTN_COUNT; i++) {
            bool pressed = (menu_pressed == i && menu_hover == i);
            bool hovered = (menu_hover == i);
            if (pressed)      sdtx_color3f(1.0f, 1.0f, 1.0f);
            else if (hovered) sdtx_color3f(1.0f, 1.0f, 0.3f);
            else              sdtx_color3f(0.6f, 0.6f, 0.6f);
            int label_len = (int)strlen(btn_labels[i]);
            float by = btn_y0 + i * (btn_h + btn_gap);
            float tx = margin_x + (btn_w - label_len * char_px) * 0.5f;
            float ty = by + (btn_h - char_px) * 0.5f;
            sdtx_pos(tx / char_px, ty / char_px);
            sdtx_puts(btn_labels[i]);
        }

        float help_y = btn_y0 + MAIN_MENU_BTN_COUNT * (btn_h + btn_gap) + char_px;
        sdtx_origin(margin_x / char_px, help_y / char_px);
        sdtx_home();
        sdtx_color3f(0.4f, 0.4f, 0.5f);
        sdtx_puts("WASD to move, Mouse to look\n");
        sdtx_puts("Click to lock mouse\n");
        sdtx_puts("LMB = break, RMB = place\n");
        sdtx_puts("Space = jump, ESC = unlock\n");
        sdtx_puts("Double-tap Space = jetpack\n");
        sdtx_puts("Scroll wheel = jetpack speed\n");
#ifdef __EMSCRIPTEN__
        sdtx_puts("C = jetpack descend\n");
#else
        sdtx_puts("Ctrl = jetpack descend\n");
#endif
        sdtx_puts("Alt+P = screenshot\n");
        sdtx_puts("L = toggle LOD debug colors\n");
        sdtx_puts("V = toggle verbose logs\n");
        sdtx_puts("R = reload config.yaml\n");

        sdtx_draw();
        sg_end_pass();
        sg_commit();
        #undef MAIN_MENU_BTN_COUNT
        return;
    }

    // ---- STATE_WORLD_SELECT: World list + [New World] ----
    if (app.state == STATE_WORLD_SELECT) {
        sg_begin_pass(&(sg_pass){
            .action = { .colors[0] = { .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f} } },
            .swapchain = sglue_swapchain()
        });

        int total_btns = app.world_list.count + 1;  // worlds + [New World]
        float btn_h = 3.0f * char_px;
        float btn_gap = 8.0f;
        float btn_w = 28.0f * char_px;
        float btn_y0 = title_y + 2.0f * char_px;

        ws_hover = hit_test_buttons(total_btns, btn_w, btn_h, btn_gap,
                                     margin_x, btn_y0, menu_mouse_x, menu_mouse_y);

        // Draw button rects (reuse helper but with green tint for [New World])
        sgl_defaults();
        sgl_matrix_mode_projection();
        sgl_load_identity();
        sgl_ortho(0.0f, win_w, win_h, 0.0f, -1.0f, 1.0f);
        sgl_matrix_mode_modelview();
        sgl_load_identity();

        for (int i = 0; i < total_btns; i++) {
            bool is_pressed = (ws_pressed == i && ws_hover == i);
            bool is_hovered = (ws_hover == i);
            bool is_new = (i == app.world_list.count);
            float x0 = margin_x, y0 = btn_y0 + i * (btn_h + btn_gap);
            float x1 = x0 + btn_w, y1 = y0 + btn_h;

            float fr, fg, fb, fa;
            if (is_pressed)      { fr = 0.25f; fg = 0.25f; fb = 0.3f; fa = 0.9f; }
            else if (is_hovered) { fr = is_new ? 0.08f : 0.12f; fg = is_new ? 0.18f : 0.14f; fb = 0.2f; fa = 0.8f; }
            else                 { fr = 0.08f; fg = 0.08f; fb = 0.12f; fa = 0.6f; }
            sgl_begin_quads();
            sgl_c4f(fr, fg, fb, fa);
            sgl_v2f(x0, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y1); sgl_v2f(x0, y1);
            sgl_end();

            float br, bg_, bb;
            if (is_new) {
                if (is_pressed)      { br = 0.5f; bg_ = 1.0f; bb = 0.5f; }
                else if (is_hovered) { br = 0.3f; bg_ = 1.0f; bb = 0.3f; }
                else                 { br = 0.2f; bg_ = 0.5f; bb = 0.2f; }
            } else {
                if (is_pressed)      { br = 1.0f; bg_ = 1.0f; bb = 1.0f; }
                else if (is_hovered) { br = 1.0f; bg_ = 1.0f; bb = 0.3f; }
                else                 { br = 0.3f; bg_ = 0.3f; bb = 0.4f; }
            }
            float bw = 2.0f;
            sgl_begin_quads();
            sgl_c3f(br, bg_, bb);
            sgl_v2f(x0, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y0 + bw); sgl_v2f(x0, y0 + bw);
            sgl_v2f(x0, y1 - bw); sgl_v2f(x1, y1 - bw); sgl_v2f(x1, y1); sgl_v2f(x0, y1);
            sgl_v2f(x0, y0); sgl_v2f(x0 + bw, y0); sgl_v2f(x0 + bw, y1); sgl_v2f(x0, y1);
            sgl_v2f(x1 - bw, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y1); sgl_v2f(x1 - bw, y1);
            sgl_end();
        }
        sgl_draw();

        // Text
        sdtx_canvas(canvas_w, canvas_h);
        sdtx_origin(0.0f, 0.0f);
        sdtx_font(0);
        sdtx_pos(margin_x / char_px, title_y / char_px);
        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("Select World    [ESC] Back");

        for (int i = 0; i < total_btns; i++) {
            bool pressed = (ws_pressed == i && ws_hover == i);
            bool hovered = (ws_hover == i);
            float by = btn_y0 + i * (btn_h + btn_gap);
            float tx = margin_x + 2.0f * char_px;
            float ty = by + (btn_h - char_px) * 0.5f;
            sdtx_pos(tx / char_px, ty / char_px);

            if (i < app.world_list.count) {
                if (pressed)      sdtx_color3f(1.0f, 1.0f, 1.0f);
                else if (hovered) sdtx_color3f(1.0f, 1.0f, 0.3f);
                else              sdtx_color3f(0.6f, 0.6f, 0.6f);
                sdtx_puts(app.world_list.worlds[i].name);
            } else {
                if (pressed)      sdtx_color3f(0.5f, 1.0f, 0.5f);
                else if (hovered) sdtx_color3f(0.3f, 1.0f, 0.3f);
                else              sdtx_color3f(0.2f, 0.7f, 0.2f);
                sdtx_puts("[+] New World");
            }
        }

        sdtx_draw();
        sg_end_pass();
        sg_commit();
        return;
    }

    // ---- STATE_SPAWN_SELECT: Summit / Twilight ----
    if (app.state == STATE_SPAWN_SELECT) {
        sg_begin_pass(&(sg_pass){
            .action = { .colors[0] = { .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f} } },
            .swapchain = sglue_swapchain()
        });

        const int sp_count = 2;
        float btn_h = 3.0f * char_px;
        float btn_gap = 8.0f;
        float btn_w = 28.0f * char_px;
        float btn_y0 = title_y + 2.0f * char_px;

        sp_hover = hit_test_buttons(sp_count, btn_w, btn_h, btn_gap,
                                     margin_x, btn_y0, menu_mouse_x, menu_mouse_y);

        sgl_defaults();
        sgl_matrix_mode_projection();
        sgl_load_identity();
        sgl_ortho(0.0f, win_w, win_h, 0.0f, -1.0f, 1.0f);
        sgl_matrix_mode_modelview();
        sgl_load_identity();

        const char* sp_labels[] = { "[1] Summit", "[2] Twilight" };
        float sp_colors[][3] = {
            { 1.0f, 0.85f, 0.3f },  // Summit: warm gold
            { 0.4f, 0.5f, 1.0f },   // Twilight: cool blue
        };

        for (int i = 0; i < sp_count; i++) {
            bool is_pressed = (sp_pressed == i && sp_hover == i);
            bool is_hovered = (sp_hover == i);
            float x0 = margin_x, y0 = btn_y0 + i * (btn_h + btn_gap);
            float x1 = x0 + btn_w, y1 = y0 + btn_h;

            float fr, fg, fb, fa;
            if (is_pressed)      { fr = 0.25f; fg = 0.25f; fb = 0.3f; fa = 0.9f; }
            else if (is_hovered) { fr = 0.12f; fg = 0.12f; fb = 0.18f; fa = 0.8f; }
            else                 { fr = 0.08f; fg = 0.08f; fb = 0.12f; fa = 0.6f; }
            sgl_begin_quads();
            sgl_c4f(fr, fg, fb, fa);
            sgl_v2f(x0, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y1); sgl_v2f(x0, y1);
            sgl_end();

            float br = sp_colors[i][0], bg_ = sp_colors[i][1], bb = sp_colors[i][2];
            if (is_pressed)      { br *= 0.5f; bg_ *= 0.5f; bb *= 0.5f; }
            else if (!is_hovered) { br *= 0.4f; bg_ *= 0.4f; bb *= 0.4f; }
            float bw = 2.0f;
            sgl_begin_quads();
            sgl_c3f(br, bg_, bb);
            sgl_v2f(x0, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y0 + bw); sgl_v2f(x0, y0 + bw);
            sgl_v2f(x0, y1 - bw); sgl_v2f(x1, y1 - bw); sgl_v2f(x1, y1); sgl_v2f(x0, y1);
            sgl_v2f(x0, y0); sgl_v2f(x0 + bw, y0); sgl_v2f(x0 + bw, y1); sgl_v2f(x0, y1);
            sgl_v2f(x1 - bw, y0); sgl_v2f(x1, y0); sgl_v2f(x1, y1); sgl_v2f(x1 - bw, y1);
            sgl_end();
        }
        sgl_draw();

        sdtx_canvas(canvas_w, canvas_h);
        sdtx_origin(0.0f, 0.0f);
        sdtx_font(0);
        sdtx_pos(margin_x / char_px, title_y / char_px);
        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("Choose Starting Location    [ESC] Back");

        for (int i = 0; i < sp_count; i++) {
            bool pressed = (sp_pressed == i && sp_hover == i);
            bool hovered = (sp_hover == i);
            float by = btn_y0 + i * (btn_h + btn_gap);
            float tx = margin_x + 2.0f * char_px;
            float ty = by + (btn_h - char_px) * 0.5f;
            sdtx_pos(tx / char_px, ty / char_px);

            float cr = sp_colors[i][0], cg = sp_colors[i][1], cb = sp_colors[i][2];
            if (pressed)      sdtx_color3f(1.0f, 1.0f, 1.0f);
            else if (hovered) sdtx_color3f(cr, cg, cb);
            else              sdtx_color3f(cr * 0.6f, cg * 0.6f, cb * 0.6f);
            sdtx_puts(sp_labels[i]);
        }

        sdtx_draw();
        sg_end_pass();
        sg_commit();
        return;
    }

    // ---- STATE_MULTI_MENU: Host / Join ----
    if (app.state == STATE_MULTI_MENU) {
        sg_begin_pass(&(sg_pass){
            .action = { .colors[0] = { .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f} } },
            .swapchain = sglue_swapchain()
        });

        #define MULTI_BTN_COUNT 2
        const char* btn_labels[MULTI_BTN_COUNT] = { "[1] Host Game", "[2] Join Game" };
        float btn_h = 3.0f * char_px;
        float btn_gap = 10.0f;
        float btn_w = 22.0f * char_px;
        float btn_y0 = title_y + 2.0f * char_px;

        multi_hover = hit_test_buttons(MULTI_BTN_COUNT, btn_w, btn_h, btn_gap,
                                        margin_x, btn_y0, menu_mouse_x, menu_mouse_y);
        draw_buttons(MULTI_BTN_COUNT, btn_w, btn_h, btn_gap, margin_x, btn_y0,
                     win_w, win_h, multi_hover, multi_pressed);

        sdtx_canvas(canvas_w, canvas_h);
        sdtx_origin(0.0f, 0.0f);
        sdtx_font(0);
        sdtx_pos(margin_x / char_px, title_y / char_px);
        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("Multiplayer    [ESC] Back");

        for (int i = 0; i < MULTI_BTN_COUNT; i++) {
            bool pressed = (multi_pressed == i && multi_hover == i);
            bool hovered = (multi_hover == i);
            if (pressed)      sdtx_color3f(1.0f, 1.0f, 1.0f);
            else if (hovered) sdtx_color3f(1.0f, 1.0f, 0.3f);
            else              sdtx_color3f(0.6f, 0.6f, 0.6f);
            int label_len = (int)strlen(btn_labels[i]);
            float by = btn_y0 + i * (btn_h + btn_gap);
            float tx = margin_x + (btn_w - label_len * char_px) * 0.5f;
            float ty = by + (btn_h - char_px) * 0.5f;
            sdtx_pos(tx / char_px, ty / char_px);
            sdtx_puts(btn_labels[i]);
        }

        sdtx_draw();
        sg_end_pass();
        sg_commit();
        #undef MULTI_BTN_COUNT
        return;
    }

    // ---- STATE_MULTI_HOST: show room code, wait for players ----
    if (app.state == STATE_MULTI_HOST) {
        lobby_poll(&app.lobby_req);

        sg_begin_pass(&(sg_pass){
            .action = { .colors[0] = { .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f} } },
            .swapchain = sglue_swapchain()
        });
        sdtx_canvas(canvas_w, canvas_h);
        sdtx_origin(0.0f, 0.0f);
        sdtx_font(0);
        sdtx_pos(margin_x / char_px, title_y / char_px);
        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("Host Game    [ESC] Cancel");

        float content_y = title_y + 3.0f * char_px;
        sdtx_pos(margin_x / char_px, content_y / char_px);

        if (app.lobby_req.state == LOBBY_PENDING) {
            sdtx_color3f(1.0f, 1.0f, 0.5f);
            sdtx_puts("Creating lobby...");
        } else if (app.lobby_req.state == LOBBY_SUCCESS) {
            sdtx_color3f(0.5f, 1.0f, 0.5f);
            sdtx_puts("Room Code:\n\n");
            sdtx_pos((margin_x + 2.0f * char_px) / char_px, (content_y + 3.0f * char_px) / char_px);
            sdtx_color3f(1.0f, 1.0f, 1.0f);
            for (int i = 0; i < (int)strlen(app.lobby_req.lobby_id); i++) {
                sdtx_putc(app.lobby_req.lobby_id[i]);
                sdtx_putc(' ');
            }
            sdtx_pos(margin_x / char_px, (content_y + 6.0f * char_px) / char_px);
            sdtx_color3f(0.5f, 0.5f, 0.6f);
            sdtx_puts("Waiting for players...");
            sdtx_pos(margin_x / char_px, (content_y + 8.0f * char_px) / char_px);
            sdtx_color3f(0.4f, 0.4f, 0.5f);
            sdtx_puts("Share this code with friends");
        } else if (app.lobby_req.state == LOBBY_ERROR) {
            sdtx_color3f(1.0f, 0.3f, 0.3f);
            sdtx_printf("Error: %s", app.lobby_req.error);
        }

        sdtx_draw();
        sg_end_pass();
        sg_commit();
        return;
    }

    // ---- STATE_MULTI_JOIN: room code text input ----
    if (app.state == STATE_MULTI_JOIN) {
        lobby_poll(&app.lobby_req);
        room_code_cursor_blink++;

        sg_begin_pass(&(sg_pass){
            .action = { .colors[0] = { .load_action = SG_LOADACTION_CLEAR,
                .clear_value = {0.05f, 0.05f, 0.1f, 1.0f} } },
            .swapchain = sglue_swapchain()
        });
        sdtx_canvas(canvas_w, canvas_h);
        sdtx_origin(0.0f, 0.0f);
        sdtx_font(0);
        sdtx_pos(margin_x / char_px, title_y / char_px);
        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("Join Game    [ESC] Back");

        float content_y = title_y + 3.0f * char_px;

        if (app.lobby_req.state == LOBBY_PENDING) {
            sdtx_pos(margin_x / char_px, content_y / char_px);
            sdtx_color3f(1.0f, 1.0f, 0.5f);
            sdtx_puts("Joining lobby...");
        } else if (app.lobby_req.state == LOBBY_SUCCESS) {
            sdtx_pos(margin_x / char_px, content_y / char_px);
            sdtx_color3f(0.5f, 1.0f, 0.5f);
            sdtx_puts("Joined! Loading world...");
        } else if (app.lobby_req.state == LOBBY_ERROR) {
            sdtx_pos(margin_x / char_px, content_y / char_px);
            sdtx_color3f(1.0f, 0.3f, 0.3f);
            sdtx_printf("Error: %s", app.lobby_req.error);
            sdtx_pos(margin_x / char_px, (content_y + 2.0f * char_px) / char_px);
            sdtx_color3f(0.5f, 0.5f, 0.6f);
            sdtx_puts("Press any key to retry");
        } else {
            // LOBBY_IDLE: show input field
            sdtx_pos(margin_x / char_px, content_y / char_px);
            sdtx_color3f(0.7f, 0.7f, 0.8f);
            sdtx_puts("Enter Room Code:");

            // Input field background
            sgl_defaults();
            sgl_matrix_mode_projection();
            sgl_load_identity();
            sgl_ortho(0.0f, win_w, win_h, 0.0f, -1.0f, 1.0f);
            sgl_matrix_mode_modelview();
            sgl_load_identity();

            float input_y = content_y + 2.5f * char_px;
            float input_w = 14.0f * char_px;
            float input_h = 2.5f * char_px;
            sgl_begin_quads();
            sgl_c4f(0.1f, 0.1f, 0.15f, 0.8f);
            sgl_v2f(margin_x, input_y);
            sgl_v2f(margin_x + input_w, input_y);
            sgl_v2f(margin_x + input_w, input_y + input_h);
            sgl_v2f(margin_x, input_y + input_h);
            sgl_end();
            float bw = 2.0f;
            sgl_begin_quads();
            sgl_c3f(0.4f, 0.4f, 0.6f);
            sgl_v2f(margin_x, input_y); sgl_v2f(margin_x + input_w, input_y);
            sgl_v2f(margin_x + input_w, input_y + bw); sgl_v2f(margin_x, input_y + bw);
            sgl_v2f(margin_x, input_y + input_h - bw); sgl_v2f(margin_x + input_w, input_y + input_h - bw);
            sgl_v2f(margin_x + input_w, input_y + input_h); sgl_v2f(margin_x, input_y + input_h);
            sgl_end();
            sgl_draw();

            // Room code characters
            float text_x = margin_x + 1.5f * char_px;
            float text_y = input_y + (input_h - char_px) * 0.5f;
            sdtx_pos(text_x / char_px, text_y / char_px);
            sdtx_color3f(1.0f, 1.0f, 1.0f);
            for (int i = 0; i < room_code_len; i++) {
                sdtx_putc(room_code[i]);
                sdtx_putc(' ');
            }
            // Blinking cursor
            if ((room_code_cursor_blink / 30) % 2 == 0 && room_code_len < 6) {
                sdtx_color3f(0.5f, 0.5f, 0.8f);
                sdtx_putc('_');
            }

            sdtx_pos(margin_x / char_px, (input_y + input_h + 2.0f * char_px) / char_px);
            sdtx_color3f(0.4f, 0.4f, 0.5f);
            sdtx_puts("Type 6-character code, Enter to join");
        }

        sdtx_draw();
        sg_end_pass();
        sg_commit();
        return;
    }

    // ---- STATE_LOADING ----
    if (app.state == STATE_LOADING) {
        if (app.loading_phase == 0) {
            // Phase 0: spawn background thread for planet_init
            printf("[GAME] PHASE 0 ENTER: spawning planet_init thread\n");
            fflush(stdout);

            app.init_task.planet = &app.planet;
            app.init_task.subdivision = 64;
            app.init_task.completed = 0;

#ifdef __EMSCRIPTEN__
            // Emscripten: no pthreads — run planet_init synchronously
            planet_init(app.init_task.planet, app.init_task.subdivision);
            app.init_task.completed = 1;
            app.loading_phase = 2;  // Skip phase 1 (thread polling)
            printf("[GAME] PHASE 0 DONE: planet_init completed synchronously, moving to phase 2\n"); fflush(stdout);
            draw_loading_screen("Generating planet...", 0, 0, 0);
            return;
#elif defined(_WIN32)
            app.init_thread = CreateThread(NULL, 0, planet_init_thread, &app.init_task, 0, NULL);
#else
            pthread_create(&app.init_thread, NULL, planet_init_thread, &app.init_task);
#endif
            app.loading_phase = 1;
            printf("[GAME] PHASE 0 DONE: thread spawned, moving to phase 1\n"); fflush(stdout);
            draw_loading_screen("Generating planet...", 0, 0, 0);
            return;
        }
        if (app.loading_phase == 1) {
            // Phase 1: poll for planet_init completion (non-blocking)
            if (app.init_task.completed) {
                // Join thread
#ifdef _WIN32
                WaitForSingleObject(app.init_thread, INFINITE);
                CloseHandle(app.init_thread);
#else
                pthread_join(app.init_thread, NULL);
#endif
                printf("[GAME] PHASE 1 DONE: planet ready, %d cells, moving to phase 2\n", app.planet.cell_count);
                fflush(stdout);
                app.loading_phase = 2;
            }
            draw_loading_screen("Generating planet...", 0, 0, 0);
            return;
        }
        if (app.loading_phase == 2) {
            // Phase 2: render_init (GPU resources — must be on main thread)
            printf("[GAME] PHASE 2 ENTER: calling render_init...\n"); fflush(stdout);
            render_init(&app.renderer, &app.planet, &app.camera, app.active_edits_dir);
            printf("[GAME] PHASE 2: render_init done\n"); fflush(stdout);

            // Initialize solar system (moons + Tenebris fallback mesh)
            solar_system_init(&app.renderer.solar_system);
            {
                float surface_r = app.planet.radius + app.planet.sea_level * app.planet.layer_thickness;
                solar_system_generate_planet_mesh(&app.renderer.solar_system, surface_r, 42);
            }
            printf("[GAME] Solar system initialized: %d moons\n",
                   app.renderer.solar_system.moon_count);
            fflush(stdout);

            // Try loading saved player position; fall back to spawn type selection
            if (!player_load()) {
                HMM_Vec3 spawn_dir;
                if (app.spawn_type == 0) {
                    // Summit: planet origin (top)
                    spawn_dir = (HMM_Vec3){{0.0f, 1.0f, 0.0f}};
                } else {
                    // Twilight zone, grass biome (lat 68.55, lon -106.52)
                    spawn_dir = HMM_NormV3((HMM_Vec3){{-0.103983f, 0.930767f, -0.350514f}});
                }
                double surface_r = lod_tree_terrain_height(&app.renderer.lod_tree, spawn_dir) + 10.0;
                app.camera.pos_d[0] = (double)spawn_dir.X * surface_r;
                app.camera.pos_d[1] = (double)spawn_dir.Y * surface_r;
                app.camera.pos_d[2] = (double)spawn_dir.Z * surface_r;
                app.camera.position = (HMM_Vec3){{
                    (float)app.camera.pos_d[0],
                    (float)app.camera.pos_d[1],
                    (float)app.camera.pos_d[2]
                }};
            }

            // Set floating origin to spawn/loaded position BEFORE any LOD meshes are generated.
            app.renderer.lod_tree.world_origin[0] = app.camera.pos_d[0];
            app.renderer.lod_tree.world_origin[1] = app.camera.pos_d[1];
            app.renderer.lod_tree.world_origin[2] = app.camera.pos_d[2];

            double cam_r = sqrt(app.camera.pos_d[0]*app.camera.pos_d[0] +
                               app.camera.pos_d[1]*app.camera.pos_d[1] +
                               app.camera.pos_d[2]*app.camera.pos_d[2]);
            printf("[GAME] PHASE 2 DONE: camera r=%.0f, origin set, moving to phase 3\n", cam_r);
            fflush(stdout);

            app.loading_phase = 3;
            draw_loading_screen("Building terrain...", 0, 20, 0);
            return;
        }
        if (app.loading_phase == 3) {
            // Phase 3: one update+upload per frame, tree builds progressively
            app.load_frame_count++;

            HMM_Mat4 vp = HMM_MulM4(app.camera.proj, app.camera.view);
            lod_tree_update(&app.renderer.lod_tree, app.camera.position, vp);
            lod_tree_upload_meshes(&app.renderer.lod_tree);

            int active = lod_tree_active_leaves(&app.renderer.lod_tree);
            int total_verts = app.renderer.lod_tree.total_vertex_count;
            int pending = job_system_pending(app.renderer.lod_tree.jobs);

            hex_terrain_update(&app.renderer.hex_terrain, app.camera.position,
                               app.renderer.lod_tree.world_origin);
            hex_terrain_upload_meshes(&app.renderer.hex_terrain);

            bool hex_ready = true;
            for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
                HexChunk* chunk = &app.renderer.hex_terrain.chunks[i];
                if (chunk->active && (chunk->generating || chunk->dirty)) {
                    hex_ready = false;
                    break;
                }
            }

            printf("[LOAD] frame %d: %d patches, %d pending, hex=%s, %dk verts\n",
                app.load_frame_count, active, pending,
                hex_ready ? "ready" : "loading", total_verts / 1000);
            fflush(stdout);

            int display_total = active + pending;
            if (display_total < 100) display_total = 100;
            double elapsed_ms = stm_ms(stm_diff(stm_now(), app.load_start_time));
            char status_buf[128];
            snprintf(status_buf, sizeof(status_buf),
                "Building terrain... (%d patches, %d pending, %.1fs)",
                active, pending, elapsed_ms / 1000.0);
            draw_loading_screen(status_buf, active, display_total, total_verts);

            bool tree_stable = pending == 0 && app.renderer.lod_tree.splits_this_frame == 0;
            bool patches_ready = active >= 20 && tree_stable;
            if (patches_ready && hex_ready) {
                printf("[GAME] Terrain ready: %d patches, %d vertices. Load time: %.0fms (%d frames). Starting game.\n",
                    active, total_verts, elapsed_ms, app.load_frame_count);
                fflush(stdout);
                app.state = STATE_PLAYING;
            }
            return;
        }
    }

    // ---- STATE_PLAYING ----
    {
        uint64_t t_frame_start = stm_now();

        solar_system_update(&app.renderer.solar_system, (double)dt);
        touch_apply_to_camera(&app.touch, &app.camera, dt);
        camera_update(&app.camera, &app.planet, &app.renderer.lod_tree, &app.renderer.hex_terrain, dt,
                      &app.renderer.solar_system);
        uint64_t t_cam = stm_now();

        app.renderer.profile.accum_camera_ms += (float)stm_ms(stm_diff(t_cam, t_frame_start));

        if (++autosave_counter >= 300) {
            autosave_counter = 0;
            player_save();
        }

        // ---- LOD retarget: switch body when gravity changes ----
        if (app.camera.gravity_body != app.renderer.lod_current_body) {
            SolarSystem* ss = &app.renderer.solar_system;
            if (app.camera.gravity_body == -1) {
                // Retarget to Tenebris — moon is no longer pinned
                double center[3] = {0, 0, 0};
                lod_tree_retarget(&app.renderer.lod_tree, LOD_BODY_PLANET,
                    center, app.planet.radius, app.planet.layer_thickness,
                    app.planet.sea_level, 42, NULL, NULL);
                ss->pinned_body = -1;
            } else {
                // Retarget to moon — pin it as stationary reference frame
                const CelestialBody* moon = &ss->moons[app.camera.gravity_body];
                lod_tree_retarget(&app.renderer.lod_tree, LOD_BODY_MOON,
                    moon->pos_d, moon->radius, 1.0f, 0,
                    moon->shape.noise_seed,
                    &moon->shape, &moon->palette);
                ss->pinned_body = app.camera.gravity_body;
                ss->pinned_center_d[0] = moon->pos_d[0];
                ss->pinned_center_d[1] = moon->pos_d[1];
                ss->pinned_center_d[2] = moon->pos_d[2];
            }
            app.renderer.lod_current_body = app.camera.gravity_body;
        }

        // Moon reference frame: sync body_center and pinned_center to the moon's
        // current Kepler position. The camera tracks orbital motion via camera_update
        // (pos_d += delta), so camera and body center stay close together.
        // world_origin is NOT shifted here — the floating-origin recenter system
        // handles that when camera drifts far enough, ensuring mesh jobs are never
        // falsely discarded as stale by the origin check.
        if (app.renderer.lod_current_body >= 0) {
            SolarSystem* ss = &app.renderer.solar_system;
            const CelestialBody* moon = &ss->moons[app.renderer.lod_current_body];

            // Sync body center to actual Kepler position (mesh generation uses this)
            app.renderer.lod_tree.body_center_d[0] = moon->pos_d[0];
            app.renderer.lod_tree.body_center_d[1] = moon->pos_d[1];
            app.renderer.lod_tree.body_center_d[2] = moon->pos_d[2];

            // Sync pinned center for SOI checks
            ss->pinned_center_d[0] = moon->pos_d[0];
            ss->pinned_center_d[1] = moon->pos_d[1];
            ss->pinned_center_d[2] = moon->pos_d[2];
        }

        lod_tree_update_origin(&app.renderer.lod_tree, app.camera.pos_d);
        render_update_mesh(&app.renderer, &app.planet, &app.camera);

        // Scan newly-loaded chunks for torch voxels and create visual instances
        {
            HexTerrain* ht = &app.renderer.hex_terrain;
            TorchSystem* ts = &app.renderer.torch_system;
            for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
                HexChunk* chunk = &ht->chunks[i];
                if (!chunk->active || !chunk->voxels || chunk->torch_scanned) continue;
                chunk->torch_scanned = true;

                torch_remove_chunk(ts, chunk->cx, chunk->cz);

                for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
                    for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
                        int min_s = chunk->col_min_solid[col][row];
                        int max_s = chunk->col_max_solid[col][row];
                        for (int l = min_s; l <= max_s; l++) {
                            if (HEX_VOXEL(chunk->voxels, col, row, l) == VOXEL_TORCH) {
                                int gcol = chunk->cx * HEX_CHUNK_SIZE + col;
                                int grow = chunk->cz * HEX_CHUNK_SIZE + row;
                                int world_layer = l + chunk->base_layer;
                                float wx, wy, wz, ux, uy, uz;
                                hex_terrain_hex_to_world(ht, gcol, grow, world_layer,
                                    &wx, &wy, &wz, &ux, &uy, &uz);
                                torch_add(ts, wx, wy, wz, ux, uy, uz,
                                          chunk->cx, chunk->cz, gcol, grow, world_layer);
                            }
                        }
                    }
                }
            }
        }

        app.renderer.hotbar_selected_slot = app.hotbar_slot;
        render_frame(&app.renderer, &app.camera, dt);
    }

    // Screenshot: capture after render_frame (which calls sg_commit)
    if (app.screenshot_requested) {
        app.screenshot_requested = false;
        screenshot_capture();
    }
}

// ---- Helper: activate world select choice ----
static void spawn_select_activate(int slot);

static void world_select_activate(int idx) {
    if (idx == app.world_list.count) {
        // [New World] button — go to spawn location selection
        app.state = STATE_SPAWN_SELECT;
        sp_hover = -1;
        sp_pressed = -1;
        sp_selected = 0;
        return;
    }
    app.active_world_idx = idx;
    if (app.is_multiplayer && app.is_host) {
        lobby_create("Host", &app.lobby_req);
        app.state = STATE_MULTI_HOST;
    } else {
        printf("[GAME] Starting single player, world=%s\n",
            app.world_list.worlds[idx].name);
        fflush(stdout);
        app.load_start_time = stm_now();
        app.load_frame_count = 0;
        start_loading();
    }
}

static void spawn_select_activate(int slot) {
    app.spawn_type = slot;
    int new_idx = world_create_new(&app.world_list);
    if (new_idx < 0) {
        app.state = STATE_WORLD_SELECT;
        return;
    }
    app.active_world_idx = new_idx;
    if (app.is_multiplayer && app.is_host) {
        lobby_create("Host", &app.lobby_req);
        app.state = STATE_MULTI_HOST;
    } else {
        printf("[GAME] Starting single player (spawn=%s), world=%s\n",
            slot == 0 ? "Summit" : "Twilight",
            app.world_list.worlds[new_idx].name);
        fflush(stdout);
        app.load_start_time = stm_now();
        app.load_frame_count = 0;
        start_loading();
    }
}

static void event(const sapp_event* ev) {
    // Track mouse/touch position for all menu states
    if (ev->type == SAPP_EVENTTYPE_MOUSE_MOVE) {
        menu_mouse_x = ev->mouse_x;
        menu_mouse_y = ev->mouse_y;
    }
    // Track touch position so hover detection works on mobile
    if (ev->type == SAPP_EVENTTYPE_TOUCHES_BEGAN ||
        ev->type == SAPP_EVENTTYPE_TOUCHES_MOVED) {
        for (int t = 0; t < ev->num_touches; t++) {
            if (ev->touches[t].changed) {
                menu_mouse_x = ev->touches[t].pos_x;
                menu_mouse_y = ev->touches[t].pos_y;
                break;
            }
        }
    }

    // Helper: detect touch began/ended as mouse-like press/release
    bool touch_began = false;
    bool touch_ended = false;
    if (ev->type == SAPP_EVENTTYPE_TOUCHES_BEGAN) touch_began = true;
    if (ev->type == SAPP_EVENTTYPE_TOUCHES_ENDED ||
        ev->type == SAPP_EVENTTYPE_TOUCHES_CANCELLED) touch_ended = true;

    // ---- STATE_MENU ----
    if (app.state == STATE_MENU) {
        #define MAIN_MENU_ACTIVATE(slot) do { \
            if ((slot) == 0) { \
                app.is_multiplayer = false; \
                app.is_host = false; \
                app.state = STATE_WORLD_SELECT; \
            } else if ((slot) == 1) { \
                app.state = STATE_MULTI_MENU; \
            } \
        } while(0)

        if ((ev->type == SAPP_EVENTTYPE_MOUSE_DOWN &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_began) {
            menu_pressed = menu_hover;
            return;
        }
        if ((ev->type == SAPP_EVENTTYPE_MOUSE_UP &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_ended) {
            if (menu_pressed >= 0 && menu_pressed == menu_hover) {
                MAIN_MENU_ACTIVATE(menu_pressed);
            }
            menu_pressed = -1;
            return;
        }
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
            if (ev->key_code == SAPP_KEYCODE_1) { MAIN_MENU_ACTIVATE(0); return; }
            if (ev->key_code == SAPP_KEYCODE_2) { MAIN_MENU_ACTIVATE(1); return; }
            if (ev->key_code == SAPP_KEYCODE_UP || ev->key_code == SAPP_KEYCODE_W) {
                menu_selected = (menu_selected - 1 + 2) % 2; return;
            }
            if (ev->key_code == SAPP_KEYCODE_DOWN || ev->key_code == SAPP_KEYCODE_S) {
                menu_selected = (menu_selected + 1) % 2; return;
            }
            if (ev->key_code == SAPP_KEYCODE_ENTER || ev->key_code == SAPP_KEYCODE_SPACE) {
                MAIN_MENU_ACTIVATE(menu_selected); return;
            }
        }
        #undef MAIN_MENU_ACTIVATE
        return;
    }

    // ---- STATE_WORLD_SELECT ----
    if (app.state == STATE_WORLD_SELECT) {
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ESCAPE) {
            app.state = app.is_multiplayer ? STATE_MULTI_MENU : STATE_MENU;
            return;
        }
        if ((ev->type == SAPP_EVENTTYPE_MOUSE_DOWN &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_began) {
            ws_pressed = ws_hover;
            return;
        }
        if ((ev->type == SAPP_EVENTTYPE_MOUSE_UP &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_ended) {
            if (ws_pressed >= 0 && ws_pressed == ws_hover) {
                world_select_activate(ws_pressed);
            }
            ws_pressed = -1;
            return;
        }
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ENTER) {
            // Quick-select: if worlds exist, pick last played; else create new
            if (app.world_list.count > 0) {
                int idx = app.world_list.active_index >= 0 ? app.world_list.active_index : 0;
                world_select_activate(idx);
            } else {
                world_select_activate(0);  // triggers [New World]
            }
            return;
        }
        return;
    }

    // ---- STATE_SPAWN_SELECT ----
    if (app.state == STATE_SPAWN_SELECT) {
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ESCAPE) {
            app.state = STATE_WORLD_SELECT;
            return;
        }

        #define SPAWN_ACTIVATE(slot) do { \
            if ((slot) >= 0 && (slot) < 2) spawn_select_activate(slot); \
        } while(0)

        if ((ev->type == SAPP_EVENTTYPE_MOUSE_DOWN &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_began) {
            sp_pressed = sp_hover;
            return;
        }
        if ((ev->type == SAPP_EVENTTYPE_MOUSE_UP &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_ended) {
            if (sp_pressed >= 0 && sp_pressed == sp_hover) {
                SPAWN_ACTIVATE(sp_pressed);
            }
            sp_pressed = -1;
            return;
        }
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
            if (ev->key_code == SAPP_KEYCODE_1) { SPAWN_ACTIVATE(0); return; }
            if (ev->key_code == SAPP_KEYCODE_2) { SPAWN_ACTIVATE(1); return; }
            if (ev->key_code == SAPP_KEYCODE_UP || ev->key_code == SAPP_KEYCODE_W) {
                sp_selected = (sp_selected - 1 + 2) % 2; return;
            }
            if (ev->key_code == SAPP_KEYCODE_DOWN || ev->key_code == SAPP_KEYCODE_S) {
                sp_selected = (sp_selected + 1) % 2; return;
            }
            if (ev->key_code == SAPP_KEYCODE_ENTER || ev->key_code == SAPP_KEYCODE_SPACE) {
                SPAWN_ACTIVATE(sp_selected); return;
            }
        }
        #undef SPAWN_ACTIVATE
        return;
    }

    // ---- STATE_MULTI_MENU ----
    if (app.state == STATE_MULTI_MENU) {
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ESCAPE) {
            app.state = STATE_MENU;
            return;
        }

        #define MULTI_ACTIVATE(slot) do { \
            if ((slot) == 0) { \
                app.is_multiplayer = true; \
                app.is_host = true; \
                app.state = STATE_WORLD_SELECT; \
            } else if ((slot) == 1) { \
                app.is_multiplayer = true; \
                app.is_host = false; \
                room_code_len = 0; \
                room_code[0] = '\0'; \
                memset(&app.lobby_req, 0, sizeof(app.lobby_req)); \
                app.state = STATE_MULTI_JOIN; \
            } \
        } while(0)

        if ((ev->type == SAPP_EVENTTYPE_MOUSE_DOWN &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_began) {
            multi_pressed = multi_hover;
            return;
        }
        if ((ev->type == SAPP_EVENTTYPE_MOUSE_UP &&
             ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) || touch_ended) {
            if (multi_pressed >= 0 && multi_pressed == multi_hover) {
                MULTI_ACTIVATE(multi_pressed);
            }
            multi_pressed = -1;
            return;
        }
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
            if (ev->key_code == SAPP_KEYCODE_1) { MULTI_ACTIVATE(0); return; }
            if (ev->key_code == SAPP_KEYCODE_2) { MULTI_ACTIVATE(1); return; }
            if (ev->key_code == SAPP_KEYCODE_ENTER || ev->key_code == SAPP_KEYCODE_SPACE) {
                MULTI_ACTIVATE(multi_hover >= 0 ? multi_hover : 0); return;
            }
        }
        #undef MULTI_ACTIVATE
        return;
    }

    // ---- STATE_MULTI_HOST ----
    if (app.state == STATE_MULTI_HOST) {
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ESCAPE) {
            app.state = STATE_MULTI_MENU;
            return;
        }
        return;
    }

    // ---- STATE_MULTI_JOIN: room code text input ----
    if (app.state == STATE_MULTI_JOIN) {
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ESCAPE) {
            app.state = STATE_MULTI_MENU;
            return;
        }

        // If error state, any key resets to idle for retry
        if (app.lobby_req.state == LOBBY_ERROR) {
            if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
                memset(&app.lobby_req, 0, sizeof(app.lobby_req));
                room_code_len = 0;
                room_code[0] = '\0';
            }
            return;
        }

        // Don't accept input while request is in flight
        if (app.lobby_req.state != LOBBY_IDLE) return;

        // Backspace
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_BACKSPACE) {
            if (room_code_len > 0) {
                room_code_len--;
                room_code[room_code_len] = '\0';
            }
            return;
        }

        // Enter: submit when 6 chars entered
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_ENTER) {
            if (room_code_len == 6) {
                lobby_join(room_code, "Player", &app.lobby_req);
            }
            return;
        }

        // Character input: A-Z, 2-9 (matching server's lobby ID charset)
        if (ev->type == SAPP_EVENTTYPE_CHAR && room_code_len < 6) {
            uint32_t c = ev->char_code;
            // Auto-uppercase
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
            // Accept A-Z (skip I, O) and 2-9
            bool valid = false;
            if (c >= 'A' && c <= 'Z' && c != 'I' && c != 'O') valid = true;
            if (c >= '2' && c <= '9') valid = true;
            if (valid) {
                room_code[room_code_len++] = (char)c;
                room_code[room_code_len] = '\0';
            }
            return;
        }
        return;
    }

    // ---- STATE_PLAYING: Touch controls ----
    if (touch_handle_event(&app.touch, ev, HOTBAR_COUNT)) {
        if (app.touch.hotbar_touch_slot >= 0 && app.touch.hotbar_touch_slot < HOTBAR_COUNT) {
            app.hotbar_slot = app.touch.hotbar_touch_slot;
            app.selected_block_type = hotbar_types[app.hotbar_slot];
        }
        if (app.touch.btn_break.just_released) {
            interact_break();
            app.touch.btn_break.just_released = false;
        }
        if (app.touch.btn_place.just_released) {
            interact_place();
            app.touch.btn_place.just_released = false;
        }
        if (app.touch.btn_menu.just_released) {
            app.touch.btn_menu.just_released = false;
            app.state = app.is_multiplayer ? STATE_MULTI_MENU : STATE_MENU;
            sapp_show_mouse(true);
            return;
        }
        return;
    }

    // ---- STATE_PLAYING ----
    // Screenshot with Alt+P
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN &&
        ev->key_code == SAPP_KEYCODE_P && (ev->modifiers & SAPP_MODIFIER_ALT)) {
        app.screenshot_requested = true;
        return;
    }

    // Toggle mouse diagnostics with Alt+M
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN &&
        ev->key_code == SAPP_KEYCODE_M && (ev->modifiers & SAPP_MODIFIER_ALT)) {
        app.camera.mouse_diag_enabled = !app.camera.mouse_diag_enabled;
        printf("[MOUSE DIAG] %s\n", app.camera.mouse_diag_enabled ? "ENABLED" : "DISABLED");
        fflush(stdout);
        // Reset counters
        app.camera.diag_frame_count = 0;
        app.camera.diag_total_events = 0;
        app.camera.diag_zero_event_frames = 0;
        app.camera.diag_max_gap_ms = 0.0f;
        app.camera.diag_max_delta = 0.0f;
        app.camera.diag_max_accum = 0.0f;
        app.camera.diag_max_frame_ms = 0.0f;
        return;
    }

    // Toggle LOD debug view with 'L' key
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_L) {
        app.renderer.show_lod_debug = !app.renderer.show_lod_debug;
        LOG(GAME, "LOD debug: %s\n", app.renderer.show_lod_debug ? "ON" : "OFF");
        return;
    }

    // Toggle verbose logging with 'V' key
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_V) {
        log_verbose = !log_verbose;
        LOG(GAME, "Verbose logging: %s\n", log_verbose ? "ON" : "OFF");
        return;
    }

    // Toggle wireframe overlay with P (physics debug)
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN &&
        ev->key_code == SAPP_KEYCODE_P && !(ev->modifiers & SAPP_MODIFIER_ALT)) {
        app.renderer.show_wireframe = !app.renderer.show_wireframe;
        LOG(GAME, "Wireframe: %s\n", app.renderer.show_wireframe ? "ON" : "OFF");
        return;
    }

    // Toggle profiler overlay with F3
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_F3) {
        app.renderer.show_profiler = !app.renderer.show_profiler;
        return;
    }

    // Hot-reload visual config with R key
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN && ev->key_code == SAPP_KEYCODE_R) {
        VisualConfig new_cfg = config_defaults();
        if (config_load(&new_cfg, "config.yaml")) {
            app.renderer.visual_config = new_cfg;
            float surface_r = app.planet.radius +
                app.planet.sea_level * app.planet.layer_thickness;
            app.renderer.atmosphere.config.rayleigh_scale = new_cfg.rayleigh_scale;
            app.renderer.atmosphere.config.mie_scale = new_cfg.mie_scale;
            app.renderer.atmosphere.config.mie_g = new_cfg.mie_g;
            app.renderer.atmosphere.config.sun_intensity = new_cfg.sun_intensity;
            app.renderer.atmosphere.config.atmosphere_radius =
                surface_r + new_cfg.atmosphere_height;
            app.renderer.lod_tree.split_factor = new_cfg.lod_split_factor;
        }
        return;
    }

    // Track ctrl state for inverted placement mode
    app.renderer.ctrl_mode = (ev->modifiers & SAPP_MODIFIER_CTRL) != 0 || app.camera.key_ctrl;

    // Block selector: keys 1-8 select hotbar slot
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
        int slot = -1;
        if (ev->key_code >= SAPP_KEYCODE_1 && ev->key_code <= SAPP_KEYCODE_8)
            slot = ev->key_code - SAPP_KEYCODE_1;
        if (slot >= 0 && slot < HOTBAR_COUNT) {
            app.hotbar_slot = slot;
            app.selected_block_type = hotbar_types[slot];
        }
    }

    // Scroll wheel cycles hotbar when mouse is locked (but not while holding space for jetpack speed)
    if (ev->type == SAPP_EVENTTYPE_MOUSE_SCROLL && app.camera.mouse_locked
        && !app.camera.key_space && !app.camera.jetpack_active && !app.camera.space_mode) {
        int dir = (ev->scroll_y > 0.0f) ? -1 : 1;
        app.hotbar_slot = (app.hotbar_slot + dir + HOTBAR_COUNT) % HOTBAR_COUNT;
        app.selected_block_type = hotbar_types[app.hotbar_slot];
    }

    // Block interactions: while mouse is locked, left=break, right=place
    if (app.camera.mouse_locked) {
        if (ev->type == SAPP_EVENTTYPE_MOUSE_DOWN) {
            if (ev->mouse_button == SAPP_MOUSEBUTTON_LEFT) {
                interact_break();
                return;
            }
            if (ev->mouse_button == SAPP_MOUSEBUTTON_RIGHT) {
                interact_place();
                return;
            }
        }
    }
    camera_handle_event(&app.camera, ev);
}

// ---- Player state persistence ----
#define PLAYER_SAVE_MAGIC 0x504C5952  // "PLYR"
typedef struct {
    uint32_t magic;
    double pos_d[3];
    float yaw;
    float pitch;
    uint8_t hotbar_slot;
    uint8_t pad[3];
} PlayerSave;

static void player_save(void) {
    char path[256];
    get_player_path(path, sizeof(path));

    // Ensure parent directory exists
    #ifdef _WIN32
    CreateDirectoryA("cache", NULL);
    if (app.active_world_idx >= 0) {
        char dir[256];
        world_get_dir(&app.world_list.worlds[app.active_world_idx], dir, sizeof(dir));
        CreateDirectoryA(dir, NULL);
    }
    #else
    mkdir("cache", 0755);
    if (app.active_world_idx >= 0) {
        char dir[256];
        world_get_dir(&app.world_list.worlds[app.active_world_idx], dir, sizeof(dir));
        mkdir(dir, 0755);
    }
    #endif

    FILE* f = fopen(path, "wb");
    if (!f) {
        printf("[SAVE] ERROR: could not open %s for writing\n", path); fflush(stdout);
        return;
    }
    PlayerSave save = {
        .magic = PLAYER_SAVE_MAGIC,
        .pos_d = { app.camera.pos_d[0], app.camera.pos_d[1], app.camera.pos_d[2] },
        .yaw = app.camera.yaw,
        .pitch = app.camera.pitch,
        .hotbar_slot = (uint8_t)app.hotbar_slot,
    };
    fwrite(&save, sizeof(save), 1, f);
    fclose(f);
    double r = sqrt(save.pos_d[0]*save.pos_d[0] + save.pos_d[1]*save.pos_d[1] + save.pos_d[2]*save.pos_d[2]);
    printf("[SAVE] Player saved to %s: r=%.0f yaw=%.2f pitch=%.2f slot=%d\n",
           path, r, save.yaw, save.pitch, save.hotbar_slot);
    fflush(stdout);
}

static bool player_load(void) {
    char path[256];
    get_player_path(path, sizeof(path));

    FILE* f = fopen(path, "rb");
    if (!f) return false;
    PlayerSave save;
    if (fread(&save, sizeof(save), 1, f) != 1 || save.magic != PLAYER_SAVE_MAGIC) {
        fclose(f);
        return false;
    }
    fclose(f);

    // Validate: position must be a reasonable distance from planet center
    double r = sqrt(save.pos_d[0]*save.pos_d[0] + save.pos_d[1]*save.pos_d[1] + save.pos_d[2]*save.pos_d[2]);
    if (r < 700000.0 || r > 900000.0) return false;  // Sanity check (~800km planet)

    app.camera.pos_d[0] = save.pos_d[0];
    app.camera.pos_d[1] = save.pos_d[1];
    app.camera.pos_d[2] = save.pos_d[2];
    app.camera.position = (HMM_Vec3){{
        (float)save.pos_d[0], (float)save.pos_d[1], (float)save.pos_d[2]
    }};
    app.camera.yaw = save.yaw;
    app.camera.pitch = save.pitch;
    app.hotbar_slot = save.hotbar_slot % HOTBAR_COUNT;
    app.selected_block_type = hotbar_types[app.hotbar_slot];

    printf("[SAVE] Player position loaded from %s (r=%.0f)\n", path, r); fflush(stdout);
    return true;
}

static void cleanup(void) {
    if (app.state == STATE_PLAYING || app.state == STATE_LOADING) {
        player_save();
        // Update last_played timestamp
        if (app.active_world_idx >= 0 && app.active_world_idx < app.world_list.count) {
            app.world_list.worlds[app.active_world_idx].last_played = (uint64_t)time(NULL);
            world_list_save(&app.world_list);
        }
        render_shutdown(&app.renderer);
        planet_destroy(&app.planet);
    }
    lobby_shutdown();
    sgl_shutdown();
    sdtx_shutdown();
}

sapp_desc sokol_main(int argc, char* argv[]) {
    (void)argc; (void)argv;
    return (sapp_desc){
        .init_cb = init,
        .frame_cb = frame,
        .cleanup_cb = cleanup,
        .event_cb = event,
        .width = 1280,
        .height = 720,
        .window_title = "Hex Planets",
        .logger.func = crash_log_func,
    };
}

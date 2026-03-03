#include "sokol_app.h"
#include "sokol_gfx.h"
#include "sokol_glue.h"
#include "sokol_log.h"
#include "sokol_time.h"
#include "util/sokol_debugtext.h"

#define HANDMADE_MATH_IMPLEMENTATION
#include "HandmadeMath.h"

#include <stdio.h>

#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#pragma warning(disable: 4996)  // freopen deprecation
#else
#include <pthread.h>
#endif

#include "camera.h"
#include "render.h"
#include "planet.h"
#include "screenshot.h"
#include "log_config.h"
#include "crash_handler.h"

typedef enum {
    STATE_MENU,
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

    // Background loading
    PlanetInitTask init_task;
#ifdef _WIN32
    HANDLE init_thread;
#else
    pthread_t init_thread;
#endif
} app;

bool log_verbose = false;

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
    VOXEL_WATER, VOXEL_ICE, VOXEL_BEDROCK, VOXEL_TORCH
};
#define HOTBAR_COUNT 8

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

static void init(void) {
    // Allocate a console window for logging (Windows GUI apps have no console by default)
    #ifdef _WIN32
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    #endif

    crash_handler_install();

    printf("[GAME] Hex Planets starting...\n");
    fflush(stdout);

    printf("[GAME] init: sg_setup...\n"); fflush(stdout);
    sg_setup(&(sg_desc){
        .environment = sglue_environment(),
        .logger.func = crash_log_func,
        .buffer_pool_size = 16384,  // LOD tree needs thousands of vertex buffers
    });

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
    app.state = STATE_MENU;
    app.last_time = stm_now();
    printf("[GAME] init: done, entering menu.\n"); fflush(stdout);
}

static void start_loading(void) {
    app.state = STATE_LOADING;
    app.loading_phase = 0;
    app.selected_block_type = VOXEL_STONE;
    app.hotbar_slot = 0;
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

static int frame_count = 0;

static void frame(void) {
    frame_count++;

    uint64_t now = stm_now();
    float dt = (float)stm_sec(stm_diff(now, app.last_time));
    app.last_time = now;
    if (dt > 0.1f) dt = 0.1f;

    if (app.state == STATE_MENU) {
        if (frame_count % 120 == 0) {
            printf("[MENU] frame %d alive\n", frame_count);
            fflush(stdout);
        }
        // Clear screen
        sg_begin_pass(&(sg_pass){
            .action = {
                .colors[0] = {
                    .load_action = SG_LOADACTION_CLEAR,
                    .clear_value = {0.05f, 0.05f, 0.1f, 1.0f}
                }
            },
            .swapchain = sglue_swapchain()
        });

        // Draw menu text
        sdtx_canvas(sapp_widthf() * 0.5f, sapp_heightf() * 0.5f);
        sdtx_origin(3.0f, 3.0f);
        sdtx_font(0);

        sdtx_color3f(0.8f, 0.9f, 1.0f);
        sdtx_puts("=== HEX PLANETS ===\n\n");

        sdtx_color3f(0.7f, 0.7f, 0.7f);
        sdtx_puts("[1] Single Player\n\n");
        sdtx_puts("[2] Multiplayer\n\n\n");

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
        return;
    }

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
            render_init(&app.renderer, &app.planet, &app.camera);
            printf("[GAME] PHASE 2: render_init done\n"); fflush(stdout);

            // Position camera above the actual terrain at spawn point
            // Twilight zone, grass biome (lat 68.55, lon -106.52)
            HMM_Vec3 spawn_dir = HMM_NormV3((HMM_Vec3){{-0.103983f, 0.930767f, -0.350514f}});
            double surface_r = lod_tree_terrain_height(&app.renderer.lod_tree, spawn_dir) + 10.0;
            app.camera.pos_d[0] = (double)spawn_dir.X * surface_r;
            app.camera.pos_d[1] = (double)spawn_dir.Y * surface_r;
            app.camera.pos_d[2] = (double)spawn_dir.Z * surface_r;
            app.camera.position = (HMM_Vec3){{
                (float)app.camera.pos_d[0],
                (float)app.camera.pos_d[1],
                (float)app.camera.pos_d[2]
            }};

            // Set floating origin to spawn position BEFORE any LOD meshes are generated.
            // Without this, the first frame of gameplay would recenter (802km > 50km threshold)
            // and destroy all meshes that were just loaded during the loading screen.
            app.renderer.lod_tree.world_origin[0] = app.camera.pos_d[0];
            app.renderer.lod_tree.world_origin[1] = app.camera.pos_d[1];
            app.renderer.lod_tree.world_origin[2] = app.camera.pos_d[2];

            printf("[GAME] PHASE 2 DONE: camera at (0, %.0f, 0), origin set, moving to phase 3\n", surface_r);
            fflush(stdout);

            app.loading_phase = 3;
            draw_loading_screen("Building terrain...", 0, 20, 0);
            return;
        }
        if (app.loading_phase == 3) {
            // Phase 3: one update+upload per frame, tree builds progressively
            // No camera movement during loading — only update projection
            app.load_frame_count++;

            // LOD tree
            HMM_Mat4 vp = HMM_MulM4(app.camera.proj, app.camera.view);
            lod_tree_update(&app.renderer.lod_tree, app.camera.position, vp);
            lod_tree_upload_meshes(&app.renderer.lod_tree);

            int active = lod_tree_active_leaves(&app.renderer.lod_tree);
            int total_verts = app.renderer.lod_tree.total_vertex_count;
            int pending = job_system_pending(app.renderer.lod_tree.jobs);

            // Hex terrain (close-range voxel grid)
            hex_terrain_update(&app.renderer.hex_terrain, app.camera.position,
                               app.renderer.lod_tree.world_origin);
            hex_terrain_upload_meshes(&app.renderer.hex_terrain);

            // Check hex terrain readiness: no chunks still generating or dirty
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

            // Show loading progress with detail
            int display_total = active + pending;
            if (display_total < 100) display_total = 100;
            double elapsed_ms = stm_ms(stm_diff(stm_now(), app.load_start_time));
            char status_buf[128];
            snprintf(status_buf, sizeof(status_buf),
                "Building terrain... (%d patches, %d pending, %.1fs)",
                active, pending, elapsed_ms / 1000.0);
            draw_loading_screen(status_buf, active, display_total, total_verts);

            // Transition once both LOD tree and hex terrain have stabilized
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

    // STATE_PLAYING
    {
        uint64_t t_frame_start = stm_now();

        camera_update(&app.camera, &app.planet, &app.renderer.lod_tree, &app.renderer.hex_terrain, dt);
        uint64_t t_cam = stm_now();

        // Feed camera timing into profiler
        app.renderer.profile.accum_camera_ms += (float)stm_ms(stm_diff(t_cam, t_frame_start));

        // Floating origin: recenter if camera has drifted far from current origin
        lod_tree_update_origin(&app.renderer.lod_tree, app.camera.pos_d);

        // Re-upload mesh if planet was modified or player moved far enough
        render_update_mesh(&app.renderer, &app.planet, &app.camera);

        // Scan newly-loaded chunks for torch voxels and create visual instances
        {
            HexTerrain* ht = &app.renderer.hex_terrain;
            TorchSystem* ts = &app.renderer.torch_system;
            for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
                HexChunk* chunk = &ht->chunks[i];
                if (!chunk->active || !chunk->voxels || chunk->torch_scanned) continue;
                chunk->torch_scanned = true;

                // Remove any old instances for this chunk position (handles re-activation)
                torch_remove_chunk(ts, chunk->cx, chunk->cz);

                // Scan for VOXEL_TORCH in this chunk's voxel data
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

static void event(const sapp_event* ev) {
    // if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
    //     printf("[EVENT] key_down: code=%d state=%d\n", ev->key_code, app.state);
    //     fflush(stdout);
    // }
    if (app.state == STATE_MENU) {
        if (ev->type == SAPP_EVENTTYPE_KEY_DOWN) {
            if (ev->key_code == SAPP_KEYCODE_1) {
                printf("[GAME] KEY '1' PRESSED - starting single player\n"); fflush(stdout);
                app.load_start_time = stm_now();
                app.load_frame_count = 0;
                start_loading();
                return;
            }
            if (ev->key_code == SAPP_KEYCODE_2) {
                LOG(GAME, "Starting multiplayer (placeholder)...\n");
                start_loading();
                return;
            }
        }
        return;
    }

    // STATE_PLAYING
    // Screenshot with Alt+P
    if (ev->type == SAPP_EVENTTYPE_KEY_DOWN &&
        ev->key_code == SAPP_KEYCODE_P && (ev->modifiers & SAPP_MODIFIER_ALT)) {
        app.screenshot_requested = true;
        return;
    }

    // Toggle LOD debug view with 'L' key (color patches by depth + show stats)
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
            // Update atmosphere config from new values
            float surface_r = app.planet.radius +
                app.planet.sea_level * app.planet.layer_thickness;
            app.renderer.atmosphere.config.rayleigh_scale = new_cfg.rayleigh_scale;
            app.renderer.atmosphere.config.mie_scale = new_cfg.mie_scale;
            app.renderer.atmosphere.config.mie_g = new_cfg.mie_g;
            app.renderer.atmosphere.config.sun_intensity = new_cfg.sun_intensity;
            app.renderer.atmosphere.config.atmosphere_radius =
                surface_r + new_cfg.atmosphere_height;
            // Update LOD split factor
            app.renderer.lod_tree.split_factor = new_cfg.lod_split_factor;
        }
        return;
    }

    // Track ctrl state for inverted placement mode
    // C key serves as Ctrl alternative (browsers intercept Ctrl for shortcuts)
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

    // Scroll wheel cycles hotbar when mouse is locked
    if (ev->type == SAPP_EVENTTYPE_MOUSE_SCROLL && app.camera.mouse_locked) {
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

static void cleanup(void) {
    if (app.state == STATE_PLAYING || app.state == STATE_LOADING) {
        render_shutdown(&app.renderer);
        planet_destroy(&app.planet);
    }
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

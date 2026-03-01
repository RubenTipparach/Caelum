#include "render.h"
#include "sokol_gfx.h"
#include "sokol_app.h"
#include "sokol_glue.h"
#include "sokol_time.h"
#include "util/sokol_debugtext.h"
#include "planet.glsl.h"
#include "sky.glsl.h"
#include "highlight.glsl.h"
#include "hex_terrain.glsl.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// LOD debug pre-draw callback: sets per-patch lod_debug uniform
typedef struct {
    planet_fs_params_t fs_params;
    float max_depth;
} LodDebugState;

static void lod_debug_pre_draw(int depth, void* user_data) {
    LodDebugState* state = (LodDebugState*)user_data;
    state->fs_params.lod_debug = (HMM_Vec4){{(float)depth, state->max_depth, 0.0f, 0.0f}};
    sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(state->fs_params));
}

void render_init(Renderer* r, Planet* planet, const Camera* cam) {
    // ---- Planet pipeline ----
    printf("[RENDER] render_init: creating planet shader...\n"); fflush(stdout);
    r->show_lod_debug = false;

    sg_shader planet_shd = sg_make_shader(planet_shader_desc(sg_query_backend()));
    printf("[RENDER] render_init: creating planet pipeline...\n"); fflush(stdout);
    r->pip = sg_make_pipeline(&(sg_pipeline_desc){
        .shader = planet_shd,
        .layout = {
            .attrs = {
                [ATTR_planet_a_position] = { .format = SG_VERTEXFORMAT_FLOAT3 },
                [ATTR_planet_a_normal]   = { .format = SG_VERTEXFORMAT_FLOAT3 },
                [ATTR_planet_a_color]    = { .format = SG_VERTEXFORMAT_FLOAT3 },
            }
        },
        .depth = {
            .compare = SG_COMPAREFUNC_LESS_EQUAL,
            .write_enabled = true,
        },
        .cull_mode = SG_CULLMODE_BACK,
        .face_winding = SG_FACEWINDING_CCW,
        .label = "planet-pipeline",
    });

    // ---- Sky pipeline (fullscreen triangle) ----
    printf("[RENDER] render_init: creating sky vertex buffer...\n"); fflush(stdout);
    float sky_verts[] = {
        -1.0f, -1.0f,
         3.0f, -1.0f,
        -1.0f,  3.0f,
    };
    r->sky_bind.vertex_buffers[0] = sg_make_buffer(&(sg_buffer_desc){
        .data = SG_RANGE(sky_verts),
        .label = "sky-vertices",
    });

    printf("[RENDER] render_init: creating sky shader...\n"); fflush(stdout);
    sg_shader sky_shd = sg_make_shader(sky_shader_desc(sg_query_backend()));
    printf("[RENDER] render_init: creating sky pipeline...\n"); fflush(stdout);
    r->sky_pip = sg_make_pipeline(&(sg_pipeline_desc){
        .shader = sky_shd,
        .layout = {
            .attrs = {
                [ATTR_sky_a_pos] = { .format = SG_VERTEXFORMAT_FLOAT2 },
            }
        },
        .depth = {
            .compare = SG_COMPAREFUNC_ALWAYS,
            .write_enabled = false,
        },
        .cull_mode = SG_CULLMODE_NONE,
        .label = "sky-pipeline",
    });

    // ---- Visual config (config.yaml) ----
    r->visual_config = config_defaults();
    config_load(&r->visual_config, "config.yaml");

    // ---- Atmosphere ----
    printf("[RENDER] render_init: creating atmosphere...\n"); fflush(stdout);
    float surface_r = planet->radius + planet->sea_level * planet->layer_thickness;
    AtmosphereConfig atmos_config = atmosphere_default_config(surface_r);
    atmos_config.atmosphere_radius = surface_r + r->visual_config.atmosphere_height;
    atmos_config.rayleigh_scale = r->visual_config.rayleigh_scale;
    atmos_config.mie_scale = r->visual_config.mie_scale;
    atmos_config.mie_g = r->visual_config.mie_g;
    atmos_config.sun_intensity = r->visual_config.sun_intensity;
    atmosphere_init(&r->atmosphere, atmos_config);

    // ---- LOD tree ----
    printf("[RENDER] render_init: creating LOD tree...\n"); fflush(stdout);
    lod_tree_init(&r->lod_tree, planet->radius, planet->layer_thickness,
                  planet->sea_level, 42);
    printf("[RENDER] render_init: LOD tree done\n"); fflush(stdout);

    // ---- Hex terrain (close-range voxel grid) ----
    printf("[RENDER] render_init: creating hex terrain...\n"); fflush(stdout);
    hex_terrain_init(&r->hex_terrain, planet->radius, planet->layer_thickness,
                     planet->sea_level, 42, r->lod_tree.jobs);
    printf("[RENDER] render_init: hex terrain done\n"); fflush(stdout);

    r->lod_tree.suppress_range = 0.0f;  // Updated dynamically each frame

    r->pass_action = (sg_pass_action){
        .colors[0] = {
            .load_action = SG_LOADACTION_CLEAR,
            .clear_value = {0.0f, 0.0f, 0.0f, 1.0f}
        }
    };

    // Sun state
    r->sun_angle = 0.0f;
    r->sun_direction = HMM_NormV3((HMM_Vec3){{1.0f, 0.3f, 0.5f}});

    // FPS counter
    r->fps_accumulator = 0.0f;
    r->fps_frame_count = 0;
    r->display_fps = 0.0f;

    // ---- Highlight wireframe pipeline ----
    printf("[RENDER] render_init: creating highlight pipeline...\n"); fflush(stdout);
    sg_shader highlight_shd = sg_make_shader(highlight_shader_desc(sg_query_backend()));
    r->highlight_pip = sg_make_pipeline(&(sg_pipeline_desc){
        .shader = highlight_shd,
        .layout = {
            .attrs = {
                [ATTR_highlight_a_position] = { .format = SG_VERTEXFORMAT_FLOAT3 },
            }
        },
        .primitive_type = SG_PRIMITIVETYPE_LINES,
        .depth = {
            .compare = SG_COMPAREFUNC_LESS_EQUAL,
            .write_enabled = false,
        },
        .cull_mode = SG_CULLMODE_NONE,
        .colors[0] = {
            .blend = {
                .enabled = true,
                .src_factor_rgb = SG_BLENDFACTOR_SRC_ALPHA,
                .dst_factor_rgb = SG_BLENDFACTOR_ONE_MINUS_SRC_ALPHA,
                .src_factor_alpha = SG_BLENDFACTOR_ONE,
                .dst_factor_alpha = SG_BLENDFACTOR_ZERO,
            }
        },
        .label = "highlight-pipeline",
    });
    r->highlight_buf = sg_make_buffer(&(sg_buffer_desc){
        .size = 12 * 3 * sizeof(float),  // 12 vertices * 3 floats
        .usage.vertex_buffer = true,
        .usage.stream_update = true,
        .label = "highlight-vertices",
    });
    r->hex_selection.valid = false;

    // ---- Hex terrain textured pipeline ----
    printf("[RENDER] render_init: loading hex terrain textures...\n"); fflush(stdout);
    {
        // Load 9 terrain textures and build atlas
        // Layout: water, sand, dirt, grass, stone, ice, snow, dirt_grass, dirt_snow
        const char* tex_files[9] = {
            "assets/textures/water.png",
            "assets/textures/sand.png",
            "assets/textures/dirt.png",
            "assets/textures/grass.png",
            "assets/textures/rocks.png",
            "assets/textures/ice.png",
            "assets/textures/snow.png",
            "assets/textures/dirt_grass.png",
            NULL,  // dirt_snow: generated from dirt+snow blend
        };
        const int TILE_SIZE = 16;
        const int NUM_TILES = 9;
        const int atlas_w = TILE_SIZE * NUM_TILES;
        const int atlas_h = TILE_SIZE;
        unsigned char* atlas_data = (unsigned char*)calloc(atlas_w * atlas_h * 4, 1);

        for (int i = 0; i < NUM_TILES; i++) {
            if (tex_files[i] == NULL) continue;  // skip generated tiles
            int w, h, channels;
            unsigned char* img = stbi_load(tex_files[i], &w, &h, &channels, 4);
            if (img) {
                // Blit into atlas at column i*TILE_SIZE
                int src_w = (w < TILE_SIZE) ? w : TILE_SIZE;
                int src_h = (h < TILE_SIZE) ? h : TILE_SIZE;
                for (int y = 0; y < src_h; y++) {
                    for (int x = 0; x < src_w; x++) {
                        int dst_idx = ((y * atlas_w) + (i * TILE_SIZE + x)) * 4;
                        int src_idx = (y * w + x) * 4;
                        atlas_data[dst_idx + 0] = img[src_idx + 0];
                        atlas_data[dst_idx + 1] = img[src_idx + 1];
                        atlas_data[dst_idx + 2] = img[src_idx + 2];
                        atlas_data[dst_idx + 3] = img[src_idx + 3];
                    }
                }
                stbi_image_free(img);
                printf("[RENDER]   loaded %s (%dx%d)\n", tex_files[i], w, h);
            } else {
                printf("[RENDER]   FAILED to load %s: %s\n", tex_files[i], stbi_failure_reason());
                // Fill with magenta for missing textures
                for (int y = 0; y < TILE_SIZE; y++) {
                    for (int x = 0; x < TILE_SIZE; x++) {
                        int idx = ((y * atlas_w) + (i * TILE_SIZE + x)) * 4;
                        atlas_data[idx + 0] = 255;
                        atlas_data[idx + 1] = 0;
                        atlas_data[idx + 2] = 255;
                        atlas_data[idx + 3] = 255;
                    }
                }
            }
        }

        // Generate dirt_snow (tile 8) by blending dirt (tile 2) and snow (tile 6)
        // Top half is snowy (70% snow), bottom half is dirty (70% dirt)
        {
            int dirt_tile = 2;
            int snow_tile = 6;
            int out_tile = 8;
            for (int y = 0; y < TILE_SIZE; y++) {
                float snow_blend = (float)y / (float)(TILE_SIZE - 1);  // 0 at top, 1 at bottom
                snow_blend = 1.0f - snow_blend;  // 1 at top (snowy), 0 at bottom (dirty)
                snow_blend = snow_blend * 0.6f + 0.2f;  // range [0.2, 0.8]
                for (int x = 0; x < TILE_SIZE; x++) {
                    int dirt_idx = ((y * atlas_w) + (dirt_tile * TILE_SIZE + x)) * 4;
                    int snow_idx = ((y * atlas_w) + (snow_tile * TILE_SIZE + x)) * 4;
                    int out_idx  = ((y * atlas_w) + (out_tile * TILE_SIZE + x)) * 4;
                    for (int c = 0; c < 4; c++) {
                        float d = (float)atlas_data[dirt_idx + c];
                        float s = (float)atlas_data[snow_idx + c];
                        atlas_data[out_idx + c] = (unsigned char)(d * (1.0f - snow_blend) + s * snow_blend);
                    }
                }
            }
            printf("[RENDER]   generated dirt_snow (blended dirt+snow)\n");
        }

        r->hex_atlas_img = sg_make_image(&(sg_image_desc){
            .width = atlas_w,
            .height = atlas_h,
            .pixel_format = SG_PIXELFORMAT_RGBA8,
            .data.mip_levels[0] = (sg_range){ atlas_data, (size_t)(atlas_w * atlas_h * 4) },
            .label = "hex-atlas-image",
        });
        free(atlas_data);

        r->hex_atlas_view = sg_make_view(&(sg_view_desc){
            .texture.image = r->hex_atlas_img,
            .label = "hex-atlas-view",
        });

        r->hex_atlas_smp = sg_make_sampler(&(sg_sampler_desc){
            .min_filter = SG_FILTER_NEAREST,
            .mag_filter = SG_FILTER_NEAREST,
            .wrap_u = SG_WRAP_CLAMP_TO_EDGE,
            .wrap_v = SG_WRAP_CLAMP_TO_EDGE,
            .label = "hex-atlas-sampler",
        });

        sg_shader hex_shd = sg_make_shader(hex_terrain_shader_desc(sg_query_backend()));
        r->hex_pip = sg_make_pipeline(&(sg_pipeline_desc){
            .shader = hex_shd,
            .layout = {
                .attrs = {
                    [ATTR_hex_terrain_a_position] = { .format = SG_VERTEXFORMAT_FLOAT3 },
                    [ATTR_hex_terrain_a_normal]   = { .format = SG_VERTEXFORMAT_FLOAT3 },
                    [ATTR_hex_terrain_a_uv]       = { .format = SG_VERTEXFORMAT_FLOAT2 },
                }
            },
            .depth = {
                .compare = SG_COMPAREFUNC_LESS_EQUAL,
                .write_enabled = true,
            },
            .cull_mode = SG_CULLMODE_BACK,
            .face_winding = SG_FACEWINDING_CCW,
            .label = "hex-terrain-pipeline",
        });
    }

    printf("[RENDER] All pipelines created\n");
    fflush(stdout);
}

void render_update_mesh(Renderer* r, Planet* planet, const Camera* cam) {
    (void)planet;

    uint64_t t0 = stm_now();

    HMM_Mat4 vp = HMM_MulM4(cam->proj, cam->view);
    lod_tree_update(&r->lod_tree, cam->position, vp);
    uint64_t t1 = stm_now();

    lod_tree_upload_meshes(&r->lod_tree);
    uint64_t t2 = stm_now();

    // Update hex terrain (close-range voxel grid).
    // Only run after LOD tree has built enough patches (avoid competing during loading).
    double hex_update_ms = 0, hex_upload_ms = 0;
    if (lod_tree_active_leaves(&r->lod_tree) >= 20) {
        hex_terrain_update(&r->hex_terrain, cam->position, r->lod_tree.world_origin);
        uint64_t t3 = stm_now();
        hex_terrain_upload_meshes(&r->hex_terrain);
        uint64_t t4 = stm_now();
        hex_update_ms = stm_ms(stm_diff(t3, t2));
        hex_upload_ms = stm_ms(stm_diff(t4, t3));

        // Suppress LOD patches where hex terrain has coverage (prevents overlap)
        r->lod_tree.suppress_range = hex_terrain_effective_range(&r->hex_terrain);
    } else {
        r->lod_tree.suppress_range = 0.0f;
    }

    // Hex selection raycast (every frame)
    r->hex_selection = hex_terrain_raycast(&r->hex_terrain,
        cam->position, cam->forward, 10.0f);

    // Build wireframe vertices if selection is valid
    if (r->hex_selection.valid) {
        float verts[36];  // 12 vertices * 3 floats
        if (hex_terrain_build_highlight(&r->hex_terrain, &r->hex_selection,
                                         r->lod_tree.world_origin, verts)) {
            sg_update_buffer(r->highlight_buf, &(sg_range){ verts, sizeof(verts) });
        } else {
            r->hex_selection.valid = false;
        }
    }

    // Accumulate profiler timings
    r->profile.accum_lod_update_ms += (float)stm_ms(stm_diff(t1, t0));
    r->profile.accum_lod_upload_ms += (float)stm_ms(stm_diff(t2, t1));
    r->profile.accum_hex_update_ms += (float)hex_update_ms;
    r->profile.accum_hex_upload_ms += (float)hex_upload_ms;
}

void render_frame(Renderer* r, const Camera* cam, float dt) {
    uint64_t render_start = stm_now();

    // FPS counter update (every 0.5 seconds)
    r->fps_accumulator += dt;
    r->fps_frame_count++;
    if (r->fps_accumulator >= 0.5f) {
        r->display_fps = (float)r->fps_frame_count / r->fps_accumulator;
        r->fps_accumulator = 0.0f;
        r->fps_frame_count = 0;
    }

    // Compute matrices
    HMM_Mat4 vp = HMM_MulM4(cam->proj, cam->view);

    // Rotation-only VP for camera-relative terrain rendering
    HMM_Mat4 view_rot = cam->view;
    view_rot.Elements[3][0] = 0.0f;
    view_rot.Elements[3][1] = 0.0f;
    view_rot.Elements[3][2] = 0.0f;
    HMM_Mat4 vp_rot = HMM_MulM4(cam->proj, view_rot);

    // Extract camera basis vectors from view matrix (column-major: Elements[col][row])
    // These are the S, U, F vectors used to build the view matrix in camera.c
    HMM_Vec3 cam_S = {{ cam->view.Elements[0][0], cam->view.Elements[1][0], cam->view.Elements[2][0] }};
    HMM_Vec3 cam_U = {{ cam->view.Elements[0][1], cam->view.Elements[1][1], cam->view.Elements[2][1] }};
    HMM_Vec3 cam_F = {{ -cam->view.Elements[0][2], -cam->view.Elements[1][2], -cam->view.Elements[2][2] }};
    float tan_half_fov = tanf(HMM_ToRad(70.0f) * 0.5f);
    float aspect = sapp_widthf() / sapp_heightf();

    // Begin render pass
    sg_begin_pass(&(sg_pass){
        .action = r->pass_action,
        .swapchain = sglue_swapchain()
    });

    // ---- 1. Draw sky (fullscreen, no depth) ----
    sky_vs_params_t sky_vs = {
        .cam_right   = (HMM_Vec4){{ cam_S.X, cam_S.Y, cam_S.Z, tan_half_fov }},
        .cam_up      = (HMM_Vec4){{ cam_U.X, cam_U.Y, cam_U.Z, aspect }},
        .cam_forward = (HMM_Vec4){{ cam_F.X, cam_F.Y, cam_F.Z, 0.0f }},
    };
    sky_fs_params_t sky_fs = {
        .sun_direction = (HMM_Vec4){{
            r->sun_direction.X, r->sun_direction.Y, r->sun_direction.Z, 0.0f
        }},
    };
    sg_apply_pipeline(r->sky_pip);
    sg_apply_bindings(&r->sky_bind);
    sg_apply_uniforms(UB_sky_vs_params, &SG_RANGE(sky_vs));
    sg_apply_uniforms(UB_sky_fs_params, &SG_RANGE(sky_fs));
    sg_draw(0, 3, 1);

    // ---- 2. Draw atmosphere (fullscreen, additive blend on top of sky) ----
    atmosphere_render(&r->atmosphere, &r->sky_bind,
                      cam->position, r->sun_direction,
                      cam_S, cam_U, cam_F, tan_half_fov, aspect,
                      r->visual_config.scale_height);

    // ---- 3. Draw planet terrain ----
    // Log depth: backend-aware Fcoef and z_bias
    float far_plane = 10000000.0f;  // 10,000 km
    float Fcoef, z_bias;
    sg_backend backend = sg_query_backend();
    if (backend == SG_BACKEND_GLCORE || backend == SG_BACKEND_GLES3) {
        // OpenGL NDC z in [-1, 1]
        Fcoef = 2.0f / log2f(far_plane + 1.0f);
        z_bias = -1.0f;
    } else {
        // D3D11, Metal, WebGPU: NDC z in [0, 1]
        Fcoef = 1.0f / log2f(far_plane + 1.0f);
        z_bias = 0.0f;
    }

    // Camera-relative rendering: reuse rotation-only VP (computed earlier for sky/atmosphere)
    HMM_Mat4 vp_terrain = vp_rot;

    // Double-float camera offset relative to floating origin.
    // Mesh vertices are stored as (world_pos - origin), so the camera offset
    // must also be (cam_pos - origin) for correct camera-relative subtraction.
    double cam_rel_origin[3] = {
        cam->pos_d[0] - r->lod_tree.world_origin[0],
        cam->pos_d[1] - r->lod_tree.world_origin[1],
        cam->pos_d[2] - r->lod_tree.world_origin[2],
    };
    float pos_hi_x = (float)cam_rel_origin[0];
    float pos_hi_y = (float)cam_rel_origin[1];
    float pos_hi_z = (float)cam_rel_origin[2];
    float pos_lo_x = (float)(cam_rel_origin[0] - (double)pos_hi_x);
    float pos_lo_y = (float)(cam_rel_origin[1] - (double)pos_hi_y);
    float pos_lo_z = (float)(cam_rel_origin[2] - (double)pos_hi_z);

    planet_vs_params_t vs_params = {
        .mvp = vp_terrain,
        .camera_offset = (HMM_Vec4){{pos_hi_x, pos_hi_y, pos_hi_z, 0.0f}},
        .camera_offset_low = (HMM_Vec4){{pos_lo_x, pos_lo_y, pos_lo_z, 0.0f}},
        .log_depth = (HMM_Vec4){{Fcoef, far_plane, z_bias, 0.0f}},
    };
    planet_fs_params_t fs_params = {
        .sun_direction = (HMM_Vec4){{
            r->sun_direction.X,
            r->sun_direction.Y,
            r->sun_direction.Z,
            r->visual_config.fog_scale_height,
        }},
        .camera_pos = (HMM_Vec4){{
            cam->position.X,
            cam->position.Y,
            cam->position.Z,
            0.0f
        }},
        .atmos_params = (HMM_Vec4){{
            r->atmosphere.config.planet_radius,
            r->atmosphere.config.atmosphere_radius,
            r->atmosphere.config.rayleigh_scale,
            r->atmosphere.config.sun_intensity,
        }},
        .lod_debug = (HMM_Vec4){{0.0f, 0.0f, 0.0f, 0.0f}},
        .dusk_sun_color = (HMM_Vec4){{
            r->visual_config.dusk_sun_color.X,
            r->visual_config.dusk_sun_color.Y,
            r->visual_config.dusk_sun_color.Z,
            0.0f,
        }},
        .day_sun_color = (HMM_Vec4){{
            r->visual_config.day_sun_color.X,
            r->visual_config.day_sun_color.Y,
            r->visual_config.day_sun_color.Z,
            0.0f,
        }},
    };

    sg_apply_pipeline(r->pip);
    sg_apply_uniforms(UB_planet_vs_params, &SG_RANGE(vs_params));
    sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(fs_params));

    if (r->show_lod_debug) {
        // LOD debug mode: per-patch callback sets depth-based color
        LodDebugState debug_state = {
            .fs_params = fs_params,
            .max_depth = (float)LOD_MAX_DEPTH,
        };
        lod_tree_render(&r->lod_tree, r->pip, vp,
                        lod_debug_pre_draw, &debug_state);
    } else {
        // Normal mode: uniforms already set, no callback needed
        lod_tree_render(&r->lod_tree, r->pip, vp, NULL, NULL);
    }

    // ---- 3b. Draw hex terrain (close-range voxel grid, textured) ----
    sg_apply_pipeline(r->hex_pip);
    {
        hex_terrain_vs_params_t hex_vs = {
            .mvp = vp_terrain,
            .camera_offset = vs_params.camera_offset,
            .camera_offset_low = vs_params.camera_offset_low,
            .log_depth = vs_params.log_depth,
        };
        hex_terrain_fs_params_t hex_fs = {
            .sun_direction = fs_params.sun_direction,
            .camera_pos = fs_params.camera_pos,
            .atmos_params = fs_params.atmos_params,
            .lod_debug = fs_params.lod_debug,
            .dusk_sun_color = fs_params.dusk_sun_color,
            .day_sun_color = fs_params.day_sun_color,
        };
        sg_apply_uniforms(UB_hex_terrain_vs_params, &SG_RANGE(hex_vs));
        sg_apply_uniforms(UB_hex_terrain_fs_params, &SG_RANGE(hex_fs));
    }
    hex_terrain_render(&r->hex_terrain, r->hex_pip,
                       r->hex_atlas_view, r->hex_atlas_smp);

    // ---- 3c. Draw hex selection wireframe ----
    if (r->hex_selection.valid) {
        sg_apply_pipeline(r->highlight_pip);

        highlight_vs_params_t hl_vs = {
            .mvp = vp_terrain,
            .camera_offset = vs_params.camera_offset,
            .camera_offset_low = vs_params.camera_offset_low,
            .log_depth = vs_params.log_depth,
        };
        highlight_fs_params_t hl_fs = {
            .color = (HMM_Vec4){{1.0f, 1.0f, 1.0f, 1.0f}},
        };
        sg_apply_uniforms(UB_highlight_vs_params, &SG_RANGE(hl_vs));
        sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(hl_fs));

        sg_bindings hl_bind = { .vertex_buffers[0] = r->highlight_buf };
        sg_apply_bindings(&hl_bind);
        sg_draw(0, 12, 1);
    }

    // ---- Profiler timing ----
    uint64_t render_end = stm_now();
    r->profile.accum_render_ms += (float)stm_ms(stm_diff(render_end, render_start));

    r->profile.accum_frames++;
    r->profile.accum_time += dt;
    if (r->profile.accum_time >= 0.5f && r->profile.accum_frames > 0) {
        float n = (float)r->profile.accum_frames;
        r->profile.lod_update_ms = r->profile.accum_lod_update_ms / n;
        r->profile.lod_upload_ms = r->profile.accum_lod_upload_ms / n;
        r->profile.hex_update_ms = r->profile.accum_hex_update_ms / n;
        r->profile.hex_upload_ms = r->profile.accum_hex_upload_ms / n;
        r->profile.camera_ms = r->profile.accum_camera_ms / n;
        r->profile.render_ms = r->profile.accum_render_ms / n;
        r->profile.frame_total_ms = r->profile.camera_ms +
            r->profile.lod_update_ms + r->profile.lod_upload_ms +
            r->profile.hex_update_ms + r->profile.hex_upload_ms +
            r->profile.render_ms;

        // Snapshot resource counts
        r->profile.lod_patches = lod_tree_active_leaves(&r->lod_tree);
        r->profile.lod_gpu_verts = r->lod_tree.total_vertex_count;
        r->profile.lod_nodes = r->lod_tree.node_count;
        r->profile.hex_chunks_active = r->hex_terrain.active_count;
        r->profile.hex_chunks_meshed = r->hex_terrain.chunks_rendered;
        r->profile.hex_gpu_verts = r->hex_terrain.total_vertex_count;
        r->profile.total_triangles = (r->profile.lod_gpu_verts + r->profile.hex_gpu_verts) / 3;
        r->profile.draw_calls = r->profile.lod_patches + r->profile.hex_chunks_meshed + 2; // +sky+atmos
        r->profile.job_pending = job_system_pending(r->lod_tree.jobs);

        // Reset accumulators
        r->profile.accum_time = 0;
        r->profile.accum_frames = 0;
        r->profile.accum_lod_update_ms = 0;
        r->profile.accum_lod_upload_ms = 0;
        r->profile.accum_hex_update_ms = 0;
        r->profile.accum_hex_upload_ms = 0;
        r->profile.accum_camera_ms = 0;
        r->profile.accum_render_ms = 0;
        r->profile.accum_frame_ms = 0;
    }

    // ---- 4. HUD overlay ----
    sdtx_canvas(sapp_widthf() * 0.5f, sapp_heightf() * 0.5f);
    sdtx_origin(0.5f, 0.5f);
    sdtx_font(0);

    sdtx_color3f(1.0f, 1.0f, 0.0f);
    sdtx_printf("FPS: %.0f\n", r->display_fps);

    // Player position + altitude
    {
        float cam_r = sqrtf(cam->position.X * cam->position.X +
                            cam->position.Y * cam->position.Y +
                            cam->position.Z * cam->position.Z);
        float altitude = cam_r - r->lod_tree.planet_radius;
        sdtx_color3f(0.6f, 0.8f, 1.0f);
        sdtx_printf("alt: %.0fm  r: %.0fkm\n", altitude, cam_r / 1000.0f);
    }

    // Jetpack status
    if (cam->jetpack_active) {
        sdtx_color3f(0.3f, 1.0f, 0.5f);
        sdtx_printf("JETPACK ON  speed: %.0fx\n", cam->jetpack_speed_mult);
    } else {
        sdtx_color3f(0.5f, 0.5f, 0.5f);
        sdtx_puts("JETPACK OFF\n");
    }

    // ---- F3 Profiler overlay ----
    if (r->show_profiler) {
        ProfileStats* p = &r->profile;
        float total = p->frame_total_ms;
        if (total < 0.01f) total = 0.01f;

        sdtx_puts("\n");
        sdtx_color3f(1.0f, 0.8f, 0.3f);
        sdtx_puts("=== PROFILER (F3) ===\n");

        // Global coordinates (double precision)
        sdtx_color3f(0.9f, 0.7f, 1.0f);
        sdtx_printf("Global: %.1f, %.1f, %.1f\n",
            cam->pos_d[0], cam->pos_d[1], cam->pos_d[2]);

        // Local coordinates (relative to floating origin)
        double lx = cam->pos_d[0] - r->lod_tree.world_origin[0];
        double ly = cam->pos_d[1] - r->lod_tree.world_origin[1];
        double lz = cam->pos_d[2] - r->lod_tree.world_origin[2];
        sdtx_printf("Local:  %.3f, %.3f, %.3f\n", lx, ly, lz);
        sdtx_printf("Origin: %.1f, %.1f, %.1f\n",
            r->lod_tree.world_origin[0],
            r->lod_tree.world_origin[1],
            r->lod_tree.world_origin[2]);

        // Frame timing breakdown
        sdtx_color3f(0.9f, 0.9f, 0.9f);
        sdtx_printf("Frame:       %5.2f ms\n", total);
        sdtx_color3f(0.7f, 0.9f, 0.7f);
        sdtx_printf("  Camera:    %5.2f ms  %4.1f%%\n", p->camera_ms, p->camera_ms / total * 100.0f);
        sdtx_printf("  LOD upd:   %5.2f ms  %4.1f%%\n", p->lod_update_ms, p->lod_update_ms / total * 100.0f);
        sdtx_printf("  LOD upl:   %5.2f ms  %4.1f%%\n", p->lod_upload_ms, p->lod_upload_ms / total * 100.0f);
        sdtx_printf("  Hex upd:   %5.2f ms  %4.1f%%\n", p->hex_update_ms, p->hex_update_ms / total * 100.0f);
        sdtx_printf("  Hex upl:   %5.2f ms  %4.1f%%\n", p->hex_upload_ms, p->hex_upload_ms / total * 100.0f);
        sdtx_printf("  Render:    %5.2f ms  %4.1f%%\n", p->render_ms, p->render_ms / total * 100.0f);

        // Draw calls and triangles
        sdtx_puts("\n");
        sdtx_color3f(0.8f, 0.8f, 1.0f);
        sdtx_printf("Draw calls:  %d\n", p->draw_calls);
        sdtx_printf("Triangles:   %dk\n", p->total_triangles / 1000);

        // Memory / resource breakdown
        sdtx_puts("\n");
        sdtx_color3f(0.6f, 1.0f, 0.8f);
        sdtx_printf("LOD patches: %d  verts: %dk  nodes: %d/%d\n",
            p->lod_patches, p->lod_gpu_verts / 1000,
            p->lod_nodes, r->lod_tree.node_capacity);
        sdtx_printf("Hex chunks:  %d/%d  verts: %dk\n",
            p->hex_chunks_meshed, p->hex_chunks_active,
            p->hex_gpu_verts / 1000);
        sdtx_printf("Jobs pending: %d\n", p->job_pending);

        // Estimated GPU memory (vertex buffers only)
        int lod_vb_mb = (p->lod_gpu_verts * (int)sizeof(LodVertex)) / (1024 * 1024);
        int hex_vb_mb = (p->hex_gpu_verts * (int)sizeof(LodVertex)) / (1024 * 1024);
        sdtx_printf("GPU vbuf:    %dMB lod + %dMB hex = %dMB\n",
            lod_vb_mb, hex_vb_mb, lod_vb_mb + hex_vb_mb);
    }

    if (r->show_lod_debug) {
        sdtx_puts("\n");
        // Per-depth geometry breakdown (shown in LOD debug mode)
        sdtx_color3f(0.5f, 0.8f, 0.8f);
        sdtx_puts("Depth  Patches    Verts   MinDist   MaxDist\n");
        sdtx_color3f(0.9f, 0.9f, 0.9f);
        for (int d = 0; d <= LOD_MAX_DEPTH; d++) {
            if (r->lod_tree.level_stats[d].patch_count == 0) continue;
            sdtx_printf(" %3d   %5d   %7d   %6.1fkm  %6.1fkm\n",
                d,
                r->lod_tree.level_stats[d].patch_count,
                r->lod_tree.level_stats[d].vertex_count,
                r->lod_tree.level_stats[d].min_distance / 1000.0f,
                r->lod_tree.level_stats[d].max_distance / 1000.0f);
        }
    }
    sdtx_draw();

    sg_end_pass();
    sg_commit();
}

void render_shutdown(Renderer* r) {
    hex_terrain_destroy(&r->hex_terrain);
    atmosphere_destroy(&r->atmosphere);
    lod_tree_destroy(&r->lod_tree);
    if (r->highlight_buf.id != SG_INVALID_ID) {
        sg_destroy_buffer(r->highlight_buf);
    }
    if (r->highlight_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->highlight_pip);
    }
    if (r->hex_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->hex_pip);
    }
    if (r->hex_atlas_view.id != SG_INVALID_ID) {
        sg_destroy_view(r->hex_atlas_view);
    }
    if (r->hex_atlas_img.id != SG_INVALID_ID) {
        sg_destroy_image(r->hex_atlas_img);
    }
    if (r->hex_atlas_smp.id != SG_INVALID_ID) {
        sg_destroy_sampler(r->hex_atlas_smp);
    }
}

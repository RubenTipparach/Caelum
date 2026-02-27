#include "render.h"
#include "sokol_gfx.h"
#include "sokol_app.h"
#include "sokol_glue.h"
#include "util/sokol_debugtext.h"
#include "planet.glsl.h"
#include "sky.glsl.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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
    r->show_lod_debug = false;

    sg_shader planet_shd = sg_make_shader(planet_shader_desc(sg_query_backend()));
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
        .label = "planet-pipeline",
    });

    // ---- Sky pipeline (fullscreen triangle) ----
    float sky_verts[] = {
        -1.0f, -1.0f,
         3.0f, -1.0f,
        -1.0f,  3.0f,
    };
    r->sky_bind.vertex_buffers[0] = sg_make_buffer(&(sg_buffer_desc){
        .data = SG_RANGE(sky_verts),
        .label = "sky-vertices",
    });

    sg_shader sky_shd = sg_make_shader(sky_shader_desc(sg_query_backend()));
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

    // ---- Atmosphere ----
    float surface_r = planet->radius + planet->sea_level * planet->layer_thickness;
    AtmosphereConfig atmos_config = atmosphere_default_config(surface_r);
    atmosphere_init(&r->atmosphere, atmos_config);

    // ---- LOD tree ----
    lod_tree_init(&r->lod_tree, planet->radius, planet->layer_thickness,
                  planet->sea_level, 42);

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

    printf("[RENDER] All pipelines created\n");
    fflush(stdout);
}

void render_update_mesh(Renderer* r, Planet* planet, const Camera* cam) {
    (void)planet;
    HMM_Mat4 vp = HMM_MulM4(cam->proj, cam->view);
    lod_tree_update(&r->lod_tree, cam->position, vp);
    lod_tree_upload_meshes(&r->lod_tree);
}

void render_frame(Renderer* r, const Camera* cam, float dt) {
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
    HMM_Mat4 inv_vp = HMM_InvGeneralM4(vp);

    // Begin render pass
    sg_begin_pass(&(sg_pass){
        .action = r->pass_action,
        .swapchain = sglue_swapchain()
    });

    // ---- 1. Draw sky (fullscreen, no depth) ----
    sky_vs_params_t sky_vs = {
        .inv_vp = inv_vp,
        .camera_pos = (HMM_Vec4){{
            cam->position.X, cam->position.Y, cam->position.Z, 0.0f
        }},
    };
    sg_apply_pipeline(r->sky_pip);
    sg_apply_bindings(&r->sky_bind);
    sg_apply_uniforms(UB_sky_vs_params, &SG_RANGE(sky_vs));
    sg_draw(0, 3, 1);

    // ---- 2. Draw atmosphere (fullscreen, additive blend on top of sky) ----
    atmosphere_render(&r->atmosphere, &r->sky_bind,
                      cam->position, r->sun_direction, inv_vp);

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

    // Camera-relative rendering: rotation-only view matrix for precision at large distances
    HMM_Mat4 view_rot = cam->view;
    view_rot.Elements[3][0] = 0.0f;
    view_rot.Elements[3][1] = 0.0f;
    view_rot.Elements[3][2] = 0.0f;
    HMM_Mat4 vp_terrain = HMM_MulM4(cam->proj, view_rot);

    planet_vs_params_t vs_params = {
        .mvp = vp_terrain,
        .camera_offset = (HMM_Vec4){{
            cam->position.X, cam->position.Y, cam->position.Z, 0.0f
        }},
        .log_depth = (HMM_Vec4){{Fcoef, far_plane, z_bias, 0.0f}},
    };
    planet_fs_params_t fs_params = {
        .sun_direction = (HMM_Vec4){{
            r->sun_direction.X,
            r->sun_direction.Y,
            r->sun_direction.Z,
            0.0f
        }},
        .camera_pos = (HMM_Vec4){{
            cam->position.X,
            cam->position.Y,
            cam->position.Z,
            0.0f
        }},
        .fog_params = (HMM_Vec4){{
            0.00005f, // fog_density (scaled for 796km planet)
            50000.0f, // fog_max_distance (50km)
            0.5f,     // inscatter_strength
            0.0f
        }},
        .lod_debug = (HMM_Vec4){{0.0f, 0.0f, 0.0f, 0.0f}},
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

    // ---- 4. HUD overlay ----
    sdtx_canvas(sapp_widthf() * 0.5f, sapp_heightf() * 0.5f);
    sdtx_origin(0.5f, 0.5f);
    sdtx_font(0);

    sdtx_color3f(1.0f, 1.0f, 0.0f);
    sdtx_printf("FPS: %.0f  patches: %d  verts: %d  nodes: %d/%d\n",
        r->display_fps,
        lod_tree_active_leaves(&r->lod_tree),
        r->lod_tree.total_vertex_count,
        r->lod_tree.node_count,
        r->lod_tree.node_capacity);

    if (r->show_lod_debug) {
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
    atmosphere_destroy(&r->atmosphere);
    lod_tree_destroy(&r->lod_tree);
}

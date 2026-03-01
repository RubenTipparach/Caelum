#include "atmosphere.h"
#include "atmosphere.glsl.h"
#include <stdio.h>

AtmosphereConfig atmosphere_default_config(float effective_surface_radius) {
    return (AtmosphereConfig){
        .planet_radius = effective_surface_radius,
        .atmosphere_radius = effective_surface_radius + 50000.0f,  // 50km atmosphere
        .rayleigh_scale = 0.015f,    // Tenebris default (scale-independent with normalized steps)
        .mie_scale = 0.01f,          // Tenebris default
        .mie_g = 0.85f,
        .sun_intensity = 5.0f,       // Tenebris default
    };
}

void atmosphere_init(Atmosphere* atmos, AtmosphereConfig config) {
    atmos->config = config;

    sg_shader shd = sg_make_shader(atmosphere_shader_desc(sg_query_backend()));

    atmos->pipeline = sg_make_pipeline(&(sg_pipeline_desc){
        .shader = shd,
        .layout = {
            .attrs = {
                [ATTR_atmosphere_a_pos] = { .format = SG_VERTEXFORMAT_FLOAT2 },
            }
        },
        .depth = {
            .compare = SG_COMPAREFUNC_ALWAYS,
            .write_enabled = false,
        },
        .cull_mode = SG_CULLMODE_NONE,
        .colors[0] = {
            .blend = {
                .enabled = true,
                .src_factor_rgb = SG_BLENDFACTOR_SRC_ALPHA,
                .dst_factor_rgb = SG_BLENDFACTOR_ONE,
                .src_factor_alpha = SG_BLENDFACTOR_ZERO,
                .dst_factor_alpha = SG_BLENDFACTOR_ONE,
            }
        },
        .label = "atmosphere-pipeline",
    });

    printf("[ATMOS] Initialized: planet_r=%.1f atmos_r=%.1f rayleigh=%.3f mie=%.3f\n",
        config.planet_radius, config.atmosphere_radius,
        config.rayleigh_scale, config.mie_scale);
    fflush(stdout);
}

void atmosphere_render(const Atmosphere* atmos, const sg_bindings* fullscreen_bind,
                       HMM_Vec3 camera_pos, HMM_Vec3 sun_dir,
                       HMM_Vec3 cam_right, HMM_Vec3 cam_up, HMM_Vec3 cam_forward,
                       float tan_half_fov, float aspect,
                       float scale_height) {
    atmosphere_vs_params_t vs = {
        .cam_right   = (HMM_Vec4){{ cam_right.X, cam_right.Y, cam_right.Z, tan_half_fov }},
        .cam_up      = (HMM_Vec4){{ cam_up.X, cam_up.Y, cam_up.Z, aspect }},
        .cam_forward = (HMM_Vec4){{ cam_forward.X, cam_forward.Y, cam_forward.Z, 0.0f }},
        .camera_pos  = (HMM_Vec4){{ camera_pos.X, camera_pos.Y, camera_pos.Z, 0.0f }},
    };
    atmosphere_fs_params_t fs = {
        .sun_direction = (HMM_Vec4){{
            sun_dir.X, sun_dir.Y, sun_dir.Z, 0.0f
        }},
        .radii = (HMM_Vec4){{
            atmos->config.planet_radius,
            atmos->config.atmosphere_radius,
            scale_height, 0.0f
        }},
        .scatter_coeffs = (HMM_Vec4){{
            atmos->config.rayleigh_scale,
            atmos->config.mie_scale,
            atmos->config.mie_g,
            atmos->config.sun_intensity,
        }},
    };

    sg_apply_pipeline(atmos->pipeline);
    sg_apply_bindings(fullscreen_bind);
    sg_apply_uniforms(UB_atmosphere_vs_params, &SG_RANGE(vs));
    sg_apply_uniforms(UB_atmosphere_fs_params, &SG_RANGE(fs));
    sg_draw(0, 3, 1);
}

void atmosphere_destroy(Atmosphere* atmos) {
    if (atmos->pipeline.id != SG_INVALID_ID) {
        sg_destroy_pipeline(atmos->pipeline);
        atmos->pipeline.id = SG_INVALID_ID;
    }
}

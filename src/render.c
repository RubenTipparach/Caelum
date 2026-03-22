#include "render.h"
#include "touch_controls.h"
#include "sokol_gfx.h"
#include "sokol_app.h"
#include "sokol_glue.h"
#include "sokol_time.h"
#include "util/sokol_debugtext.h"
#include "planet.glsl.h"
#include "sky.glsl.h"
#include "highlight.glsl.h"
#include "hex_terrain.glsl.h"
#include "hotbar.glsl.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "terrain_noise.h"
#include "ai_agent.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---- Planet texture generation ----
// Generates an equirectangular texture showing continents and oceans.
// Uses the same noise + biome coloring as LOD mesh generation.
static void generate_planet_texture(Renderer* r, float planet_radius, int seed, int width, int height) {
    printf("[RENDER] Generating %dx%d planet texture...\n", width, height); fflush(stdout);

    unsigned char* pixels = (unsigned char*)malloc(width * height * 4);

    fnl_state continental = ht_create_continental_noise(seed);
    fnl_state mountain = ht_create_mountain_noise(seed);
    fnl_state warp = ht_create_warp_noise(seed);
    fnl_state detail = ht_create_detail_noise(seed);

    // Biome reference colors (same as lod.c terrain_color_m)
    const float water_deep[3]    = {0.06f, 0.10f, 0.25f};
    const float water_shallow[3] = {0.12f, 0.21f, 0.39f};
    const float sand_col[3]      = {0.94f, 0.84f, 0.77f};
    const float grass_col[3]     = {0.07f, 0.58f, 0.35f};
    const float rock_col[3]      = {0.46f, 0.45f, 0.45f};
    const float high_stone[3]    = {0.50f, 0.50f, 0.50f};
    const float ice_col[3]       = {0.44f, 0.77f, 0.97f};

    for (int y = 0; y < height; y++) {
        float lat = ((float)y / (float)(height - 1) - 0.5f) * (float)M_PI;  // [-pi/2, pi/2]
        float cos_lat = cosf(lat);
        float sin_lat = sinf(lat);

        for (int x = 0; x < width; x++) {
            float lon = ((float)x / (float)(width - 1) - 0.5f) * 2.0f * (float)M_PI;  // [-pi, pi]

            // Unit sphere direction
            HMM_Vec3 dir = {{
                cos_lat * cosf(lon),
                sin_lat,
                cos_lat * sinf(lon),
            }};

            // Sample terrain height
            float h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, dir);
            float rel = h_m - TERRAIN_SEA_LEVEL_M;

            // Biome coloring (same logic as lod.c terrain_color_m)
            float r_c, g_c, b_c;
            if (h_m < TERRAIN_SEA_LEVEL_M) {
                float depth_t = -rel / 2000.0f;
                if (depth_t > 1.0f) depth_t = 1.0f;
                r_c = water_shallow[0] + (water_deep[0] - water_shallow[0]) * depth_t;
                g_c = water_shallow[1] + (water_deep[1] - water_shallow[1]) * depth_t;
                b_c = water_shallow[2] + (water_deep[2] - water_shallow[2]) * depth_t;
            } else {
                float blend = 100.0f;
                float t1 = ht_smoothstepf(200.0f - blend, 200.0f + blend, rel);
                float t2 = ht_smoothstepf(1500.0f - blend, 1500.0f + blend, rel);
                float t3 = ht_smoothstepf(3000.0f - blend, 3000.0f + blend, rel);
                float t4 = ht_smoothstepf(4500.0f - blend, 4500.0f + blend, rel);

                r_c = sand_col[0]; g_c = sand_col[1]; b_c = sand_col[2];
                r_c += (grass_col[0] - r_c) * t1; g_c += (grass_col[1] - g_c) * t1; b_c += (grass_col[2] - b_c) * t1;
                r_c += (rock_col[0] - r_c) * t2;  g_c += (rock_col[1] - g_c) * t2;  b_c += (rock_col[2] - b_c) * t2;
                r_c += (high_stone[0] - r_c) * t3; g_c += (high_stone[1] - g_c) * t3; b_c += (high_stone[2] - b_c) * t3;
                r_c += (ice_col[0] - r_c) * t4;   g_c += (ice_col[1] - g_c) * t4;   b_c += (ice_col[2] - b_c) * t4;
            }

            int idx = (y * width + x) * 4;
            pixels[idx + 0] = (unsigned char)(fminf(1.0f, fmaxf(0.0f, r_c)) * 255.0f);
            pixels[idx + 1] = (unsigned char)(fminf(1.0f, fmaxf(0.0f, g_c)) * 255.0f);
            pixels[idx + 2] = (unsigned char)(fminf(1.0f, fmaxf(0.0f, b_c)) * 255.0f);
            pixels[idx + 3] = 255;
        }
    }

    r->planet_tex_img = sg_make_image(&(sg_image_desc){
        .width = width,
        .height = height,
        .pixel_format = SG_PIXELFORMAT_RGBA8,
        .data.mip_levels[0] = (sg_range){ pixels, (size_t)(width * height * 4) },
        .label = "planet-color-texture",
    });
    r->planet_tex_view = sg_make_view(&(sg_view_desc){
        .texture.image = r->planet_tex_img,
        .label = "planet-color-view",
    });
    r->planet_tex_smp = sg_make_sampler(&(sg_sampler_desc){
        .min_filter = SG_FILTER_LINEAR,
        .mag_filter = SG_FILTER_LINEAR,
        .wrap_u = SG_WRAP_REPEAT,
        .wrap_v = SG_WRAP_CLAMP_TO_EDGE,
        .label = "planet-color-sampler",
    });

    free(pixels);
    printf("[RENDER] Planet texture generated.\n"); fflush(stdout);
}

// LOD debug pre-draw callback: sets per-patch lod_debug uniform
typedef struct {
    planet_fs_params_t fs_params;
    float max_depth;
} LodDebugState;

static void lod_debug_pre_draw(int depth, void* user_data) {
    LodDebugState* state = (LodDebugState*)user_data;
    state->fs_params.lod_debug.X = (float)depth;
    state->fs_params.lod_debug.Y = state->max_depth;
    // .Z and .W preserved (hex_fade_start/end)
    sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(state->fs_params));
}

void render_init(Renderer* r, Planet* planet, const Camera* cam, const char* edits_dir) {
    // ---- Planet pipeline ----
    printf("[RENDER] render_init: creating planet shader...\n"); fflush(stdout);
    r->show_lod_debug = false;
    r->lod_current_body = -1;   // Start targeting Tenebris

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
    r->lod_tree.split_factor = r->visual_config.lod_split_factor;
    printf("[RENDER] render_init: LOD tree done\n"); fflush(stdout);

    // ---- Hex terrain (close-range voxel grid) ----
    printf("[RENDER] render_init: creating hex terrain...\n"); fflush(stdout);
    hex_terrain_init(&r->hex_terrain, planet->radius, planet->layer_thickness,
                     planet->sea_level, 42, r->lod_tree.jobs, edits_dir);
    printf("[RENDER] render_init: hex terrain done\n"); fflush(stdout);

    // ---- Planet color texture (equirectangular, baked from terrain noise) ----
    generate_planet_texture(r, planet->radius, 42, 2048, 1024);

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
        .size = 36 * 3 * sizeof(float),  // 36 vertices * 3 floats (full hex prism)
        .usage.vertex_buffer = true,
        .usage.stream_update = true,
        .label = "highlight-vertices",
    });

    // ---- Wireframe line pipeline for placement face highlight (green inset) ----
    r->highlight_face_pip = sg_make_pipeline(&(sg_pipeline_desc){
        .shader = highlight_shd,
        .layout = {
            .attrs = {
                [ATTR_highlight_a_position] = { .format = SG_VERTEXFORMAT_FLOAT3 },
            }
        },
        .primitive_type = SG_PRIMITIVETYPE_LINES,
        .depth = {
            .compare = SG_COMPAREFUNC_ALWAYS,
            .write_enabled = false,
        },
        .colors[0] = {
            .blend = {
                .enabled = true,
                .src_factor_rgb = SG_BLENDFACTOR_SRC_ALPHA,
                .dst_factor_rgb = SG_BLENDFACTOR_ONE_MINUS_SRC_ALPHA,
                .src_factor_alpha = SG_BLENDFACTOR_ONE,
                .dst_factor_alpha = SG_BLENDFACTOR_ZERO,
            }
        },
        .label = "highlight-face-pipeline",
    });
    r->highlight_face_buf = sg_make_buffer(&(sg_buffer_desc){
        .size = 12 * 3 * sizeof(float),  // max 12 line verts * 3 floats (hex = 6 segments)
        .usage.vertex_buffer = true,
        .usage.stream_update = true,
        .label = "highlight-face-vertices",
    });

    r->hex_selection.valid = false;

    // ---- Normal debug arrows buffer (G key) ----
    // 2 arrows: each = shaft(2 verts) + arrowhead(6 verts) = 8 verts. Total = 16 verts * 3 floats.
    r->normal_debug_buf = sg_make_buffer(&(sg_buffer_desc){
        .size = 16 * 3 * sizeof(float),
        .usage.vertex_buffer = true,
        .usage.stream_update = true,
        .label = "normal-debug-arrows",
    });

    // ---- Hex terrain textured pipeline ----
    printf("[RENDER] render_init: loading hex terrain textures...\n"); fflush(stdout);
    {
        // Load 11 terrain textures and build atlas
        // Layout: water, sand, dirt, grass, stone, ice, snow, dirt_grass, dirt_snow, torch, moon
        const char* tex_files[11] = {
            "assets/textures/water.png",
            "assets/textures/sand.png",
            "assets/textures/dirt.png",
            "assets/textures/grass.png",
            "assets/textures/rocks.png",
            "assets/textures/ice.png",
            "assets/textures/snow.png",
            "assets/textures/dirt_grass.png",
            NULL,  // dirt_snow: generated from dirt+snow blend
            "assets/textures/torch.png",
            "assets/textures/moon.png",
        };
        const int TILE_SIZE = 16;
        const int NUM_TILES = 11;
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
                    [ATTR_hex_terrain_a_position]    = { .format = SG_VERTEXFORMAT_FLOAT3 },
                    [ATTR_hex_terrain_a_normal]      = { .format = SG_VERTEXFORMAT_FLOAT3 },
                    [ATTR_hex_terrain_a_uv]          = { .format = SG_VERTEXFORMAT_FLOAT2 },
                    [ATTR_hex_terrain_a_color]       = { .format = SG_VERTEXFORMAT_FLOAT3 },
                    [ATTR_hex_terrain_a_sky_light]   = { .format = SG_VERTEXFORMAT_FLOAT },
                    [ATTR_hex_terrain_a_torch_light] = { .format = SG_VERTEXFORMAT_FLOAT },
                }
            },
            .depth = {
                .compare = SG_COMPAREFUNC_LESS_EQUAL,
                .write_enabled = true,
                .bias = -4.0f,             // Small bias: hex wins depth test at LOD overlap boundary
                .bias_slope_scale = -1.0f,
            },
            .cull_mode = SG_CULLMODE_BACK,
            .face_winding = SG_FACEWINDING_CCW,
            .label = "hex-terrain-pipeline",
        });
    }

    // ---- Wireframe overlay pipelines (LOD debug mode) ----
    {
        sg_shader wire_shd = sg_make_shader(highlight_shader_desc(sg_query_backend()));

        // Pipeline for LodVertex (36-byte stride)
        r->lod_wireframe_pip = sg_make_pipeline(&(sg_pipeline_desc){
            .shader = wire_shd,
            .layout = {
                .buffers[0] = { .stride = 36 },
                .attrs = {
                    [ATTR_highlight_a_position] = { .format = SG_VERTEXFORMAT_FLOAT3 },
                }
            },
            .index_type = SG_INDEXTYPE_UINT32,
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
            .label = "lod-wireframe-pipeline",
        });

        // Pipeline for HexVertex (52-byte stride)
        r->hex_wireframe_pip = sg_make_pipeline(&(sg_pipeline_desc){
            .shader = wire_shd,
            .layout = {
                .buffers[0] = { .stride = 52 },
                .attrs = {
                    [ATTR_highlight_a_position] = { .format = SG_VERTEXFORMAT_FLOAT3 },
                }
            },
            .index_type = SG_INDEXTYPE_UINT32,
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
            .label = "hex-wireframe-pipeline",
        });

        // Shared wireframe index buffer: maps triangle lists to line edges.
        // For triangle (v0,v1,v2): edges (v0,v1), (v1,v2), (v2,v0).
        #define WIREFRAME_MAX_VERTS 199998  // must be multiple of 3 for triangle wireframe
        #define WIREFRAME_NUM_TRIS  (WIREFRAME_MAX_VERTS / 3)
        #define WIREFRAME_NUM_IDX   (WIREFRAME_NUM_TRIS * 6)
        uint32_t* wire_indices = (uint32_t*)malloc(WIREFRAME_NUM_IDX * sizeof(uint32_t));
        for (int i = 0; i < WIREFRAME_MAX_VERTS; i += 3) {
            int base = (i / 3) * 6;
            wire_indices[base + 0] = i;     wire_indices[base + 1] = i + 1;
            wire_indices[base + 2] = i + 1; wire_indices[base + 3] = i + 2;
            wire_indices[base + 4] = i + 2; wire_indices[base + 5] = i;
        }
        r->wireframe_idx = sg_make_buffer(&(sg_buffer_desc){
            .usage.index_buffer = true,
            .data = { wire_indices, WIREFRAME_NUM_IDX * sizeof(uint32_t) },
            .label = "wireframe-indices",
        });
        free(wire_indices);
    }

    // ---- Hotbar UI pipeline ----
    {
        sg_shader hotbar_shd = sg_make_shader(hotbar_shader_desc(sg_query_backend()));
        r->hotbar_pip = sg_make_pipeline(&(sg_pipeline_desc){
            .shader = hotbar_shd,
            .layout = {
                .attrs = {
                    [ATTR_hotbar_a_position] = { .format = SG_VERTEXFORMAT_FLOAT2 },
                    [ATTR_hotbar_a_uv]       = { .format = SG_VERTEXFORMAT_FLOAT2 },
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
                    .dst_factor_rgb = SG_BLENDFACTOR_ONE_MINUS_SRC_ALPHA,
                    .src_factor_alpha = SG_BLENDFACTOR_ONE,
                    .dst_factor_alpha = SG_BLENDFACTOR_ZERO,
                }
            },
            .label = "hotbar-pipeline",
        });

        // Dynamic vertex buffer for hotbar quads (8 slots * 6 verts * 4 floats)
        r->hotbar_buf = sg_make_buffer(&(sg_buffer_desc){
            .usage.vertex_buffer = true,
            .usage.stream_update = true,
            .size = 8 * 6 * 4 * sizeof(float),
            .label = "hotbar-vbuf",
        });
    }

    // ---- Torch rendering system ----
    torch_init(&r->torch_system);

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

    // Hex terrain: 3D voxel system (close-range, planet + moons)
    uint64_t t3 = stm_now();
    hex_terrain_update(&r->hex_terrain, cam->position, r->lod_tree.world_origin);
    uint64_t t4 = stm_now();
    hex_terrain_upload_meshes(&r->hex_terrain);
    uint64_t t5 = stm_now();
    double hex_update_ms = stm_ms(stm_diff(t4, t3));
    double hex_upload_ms = stm_ms(stm_diff(t5, t4));

    // Raycast for block selection (crosshair → hex terrain)
    r->hex_selection = hex_terrain_raycast(&r->hex_terrain,
        cam->position, cam->forward, 20.0f);

    // Placement target: same block, ctrl inverts to far side face
    r->hex_placement = r->hex_selection;
    if (r->ctrl_mode && r->hex_placement.valid) {
        hex_terrain_ctrl_placement(&r->hex_terrain, &r->hex_placement, cam->position);
    }

    // No LOD suppression — both LOD and hex terrain render, hex wins via
    // HEX_SURFACE_BIAS (0.5m) in the log depth buffer. Zero flickering/gaps.
    // suppress_range stays 0.0 (set in render_init)

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
    // Always render — Tenebris fallback mesh provides solid geometry when on a moon.
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
            r->lod_tree.suppress_range,  // hex terrain suppress range (0 = disabled)
        }},
        .atmos_params = (HMM_Vec4){{
            r->atmosphere.config.planet_radius,
            r->atmosphere.config.atmosphere_radius,
            r->atmosphere.config.rayleigh_scale * r->visual_config.terrain_fog_density,
            r->atmosphere.config.sun_intensity,
        }},
        .lod_debug = (HMM_Vec4){{0.0f, 0.0f, r->lod_fade_start, r->lod_fade_end}},
        .dusk_sun_color = (HMM_Vec4){{
            r->visual_config.dusk_sun_color.X,
            r->visual_config.dusk_sun_color.Y,
            r->visual_config.dusk_sun_color.Z,
            r->visual_config.coltab_factor,
        }},
        .day_sun_color = (HMM_Vec4){{
            r->visual_config.day_sun_color.X,
            r->visual_config.day_sun_color.Y,
            r->visual_config.day_sun_color.Z,
            r->visual_config.coltab_max_dist,
        }},
        .planet_tex_params = (HMM_Vec4){{
            (r->planet_tex_img.id != SG_INVALID_ID) ? 1.0f : 0.0f,  // enabled
            50000.0f,   // blend_start: vertex colors dominate below 50km (LOD 6+)
            200000.0f,  // blend_end: fully texture above 200km (LOD 1-5, orbital)
            0.0f,
        }},
    };

    // Snapshot Tenebris-mode fs_params BEFORE applying moon overrides.
    // This is passed to solar_system_render for the Tenebris fallback mesh.
    planet_fs_params_t planet_fs_for_tenebris = fs_params;

    // When LOD tree targets a moon, adjust uniforms for correct lighting/fog
    if (r->lod_current_body >= 0) {
        // camera_pos relative to frozen moon center (body_center_d, not live pos_d)
        float cx = (float)(cam->pos_d[0] - r->lod_tree.body_center_d[0]);
        float cy = (float)(cam->pos_d[1] - r->lod_tree.body_center_d[1]);
        float cz = (float)(cam->pos_d[2] - r->lod_tree.body_center_d[2]);
        fs_params.camera_pos = (HMM_Vec4){{cx, cy, cz, 0.0f}};
        // Disable aerial perspective fog on moons (no atmosphere)
        fs_params.atmos_params = (HMM_Vec4){{1.0f, 2.0f, 0.0f, 0.0f}};
        // Disable planet texture on moons (they use vertex colors only)
        fs_params.planet_tex_params = (HMM_Vec4){{0.0f, 0.0f, 0.0f, 0.0f}};
    }

    sg_apply_pipeline(r->pip);
    sg_apply_uniforms(UB_planet_vs_params, &SG_RANGE(vs_params));
    sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(fs_params));

    // Set up hex terrain render info for depth-13 hex mesh nodes
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
        .hex_fade = (HMM_Vec4){{
            r->visual_config.hex_fade_start,
            r->visual_config.hex_fade_end,
            0.0f, 0.0f,
        }},
    };
    LodHexRenderInfo hex_info = {
        .pip = r->hex_pip,
        .atlas_view = r->hex_atlas_view,
        .atlas_smp = r->hex_atlas_smp,
        .vs_ub = { .slot = UB_hex_terrain_vs_params, .data = &hex_vs, .size = sizeof(hex_vs) },
        .fs_ub = { .slot = UB_hex_terrain_fs_params, .data = &hex_fs, .size = sizeof(hex_fs) },
        .planet_atlas_view = r->planet_tex_view,
        .planet_atlas_smp = r->planet_tex_smp,
    };
    LodUniformBlock planet_ub[2] = {
        { .slot = UB_planet_vs_params, .data = &vs_params, .size = sizeof(vs_params) },
        { .slot = UB_planet_fs_params, .data = &fs_params, .size = sizeof(fs_params) },
    };

    // suppress_range was already computed in render_update via hex_terrain_effective_range.
    // Only suppress where hex terrain has actually meshed chunks (not the full HEX_RANGE).

    if (r->show_lod_debug) {
        LodDebugState debug_state = {
            .fs_params = fs_params,
            .max_depth = (float)LOD_MAX_DEPTH,
        };
        lod_tree_render(&r->lod_tree, r->pip, vp,
                        lod_debug_pre_draw, &debug_state, planet_ub, &hex_info);

        // Wireframe overlay: draw triangle edges on top of LOD patches
        highlight_vs_params_t wire_vs = {
            .mvp = vp_terrain,
            .camera_offset = vs_params.camera_offset,
            .camera_offset_low = vs_params.camera_offset_low,
            .log_depth = vs_params.log_depth,
        };
        highlight_fs_params_t wire_fs = {
            .color = (HMM_Vec4){{0.0f, 0.0f, 0.0f, 0.4f}},
        };
        LodUniformBlock wire_ub[2] = {
            { .slot = UB_highlight_vs_params, .data = &wire_vs, .size = sizeof(wire_vs) },
            { .slot = UB_highlight_fs_params, .data = &wire_fs, .size = sizeof(wire_fs) },
        };
        lod_tree_render_wireframe(&r->lod_tree,
            r->lod_wireframe_pip, r->hex_wireframe_pip,
            r->wireframe_idx, wire_ub);
    } else {
        lod_tree_render(&r->lod_tree, r->pip, vp, NULL, NULL, planet_ub, &hex_info);
    }

    // ---- 3a+. Draw moons ----
    // Pass the original (pre-moon-override) fs_params for the Tenebris fallback mesh
    // so it looks identical to the LOD-rendered Tenebris from space.
    solar_system_render(&r->solar_system, cam, r->sun_direction,
                        vp_terrain, r->pip,
                        Fcoef, far_plane, z_bias,
                        r->lod_tree.world_origin,
                        r->lod_current_body,
                        r->planet_tex_view, r->planet_tex_smp,
                        &planet_fs_for_tenebris);

    // ---- 3a++. Draw moon orbit lines (space mode) ----
    solar_system_draw_orbits(&r->solar_system, cam, vp_terrain);

    // ---- 3b. Draw hex terrain (close-range 3D voxel chunks) ----
    if (r->hex_terrain.active_count > 0) {
        sg_apply_pipeline(r->hex_pip);
        sg_apply_uniforms(UB_hex_terrain_vs_params, &SG_RANGE(hex_vs));
        sg_apply_uniforms(UB_hex_terrain_fs_params, &SG_RANGE(hex_fs));
        hex_terrain_render(&r->hex_terrain, r->hex_pip,
                           r->hex_atlas_view, r->hex_atlas_smp);
    }

    // ---- 3b+. Draw torch models ----
    torch_update(&r->torch_system, dt);
    if (r->torch_system.instance_count > 0) {
        double cam_off[3] = {
            (double)vs_params.camera_offset.X,
            (double)vs_params.camera_offset.Y,
            (double)vs_params.camera_offset.Z
        };
        double cam_off_low[3] = {
            (double)vs_params.camera_offset_low.X,
            (double)vs_params.camera_offset_low.Y,
            (double)vs_params.camera_offset_low.Z
        };
        torch_render(&r->torch_system, vp_terrain,
                     cam_off, cam_off_low,
                     r->lod_tree.world_origin,
                     vs_params.log_depth.X, vs_params.log_depth.Y, vs_params.log_depth.Z,
                     r->sun_direction);
    }

    // ---- 3b++. Physics wireframe overlay (P key): hex prism outlines ----
    if (r->show_wireframe && r->hex_terrain.active_count > 0) {
        sg_apply_pipeline(r->highlight_pip);
        highlight_vs_params_t wire_vs = {
            .mvp = vp_terrain,
            .camera_offset = vs_params.camera_offset,
            .camera_offset_low = vs_params.camera_offset_low,
            .log_depth = vs_params.log_depth,
        };
        highlight_fs_params_t wire_fs = {
            .color = (HMM_Vec4){{0.0f, 1.0f, 0.0f, 0.5f}},
        };
        sg_apply_uniforms(UB_highlight_vs_params, &SG_RANGE(wire_vs));
        sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(wire_fs));

        for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
            HexChunk* chunk = &r->hex_terrain.chunks[i];
            if (!chunk->active || chunk->wire_vertex_count == 0 ||
                chunk->wire_buf.id == SG_INVALID_ID) continue;
            sg_bindings wb = {0};
            wb.vertex_buffers[0] = chunk->wire_buf;
            sg_apply_bindings(&wb);
            sg_draw(0, chunk->wire_vertex_count, 1);
        }
    }

    // ---- 3c. Draw hex selection highlight (full prism outline + green placement face) ----
    if (r->hex_selection.valid) {
        float highlight_verts[108];  // 36 verts * 3 floats
        if (hex_terrain_build_highlight(&r->hex_terrain, &r->hex_selection,
                                         r->lod_tree.world_origin, highlight_verts)) {
            // White wireframe prism outline
            sg_update_buffer(r->highlight_buf, &(sg_range){
                highlight_verts, sizeof(highlight_verts)
            });
            sg_apply_pipeline(r->highlight_pip);
            highlight_vs_params_t hl_vs = {
                .mvp = vp_terrain,
                .camera_offset = vs_params.camera_offset,
                .camera_offset_low = vs_params.camera_offset_low,
                .log_depth = vs_params.log_depth,
            };
            highlight_fs_params_t hl_fs = {
                .color = (HMM_Vec4){{1.0f, 1.0f, 1.0f, 0.8f}},
            };
            sg_apply_uniforms(UB_highlight_vs_params, &SG_RANGE(hl_vs));
            sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(hl_fs));
            sg_bindings hl_bind = {0};
            hl_bind.vertex_buffers[0] = r->highlight_buf;
            sg_apply_bindings(&hl_bind);
            sg_draw(0, 36, 1);

            // Green wireframe placement face (uses placement target, not selection)
            float face_verts[36];  // max 12 verts * 3 floats
            int face_vert_count = r->hex_placement.valid ?
                hex_terrain_build_placement_face(
                    &r->hex_terrain, &r->hex_placement,
                    r->lod_tree.world_origin, face_verts) : 0;
            if (face_vert_count > 0) {
                sg_update_buffer(r->highlight_face_buf, &(sg_range){
                    face_verts, (size_t)(face_vert_count * 3) * sizeof(float)
                });
                sg_apply_pipeline(r->highlight_face_pip);
                sg_apply_uniforms(UB_highlight_vs_params, &SG_RANGE(hl_vs));
                highlight_fs_params_t face_fs = {
                    .color = (HMM_Vec4){{0.0f, 1.0f, 0.0f, 1.0f}},
                };
                sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(face_fs));
                sg_bindings face_bind = {0};
                face_bind.vertex_buffers[0] = r->highlight_face_buf;
                sg_apply_bindings(&face_bind);
                sg_draw(0, face_vert_count, 1);
            }
        }
    }

    // ---- G key: debug arrows (red=raycast origin+up, green=camera forward 50m) ----
    if (cam->show_normal_debug && cam->gravity_body >= 0 && r->hex_terrain.frame_valid) {
        float arrow_len = 50.0f;
        float head_len = 8.0f;
        float head_spread = 3.0f;
        float verts[16 * 3]; // 16 vertices * 3 floats (2 arrows x 6 verts each + spare)
        int vi = 0;

        // Helper: write vertex relative to world_origin
        #define ARROW_V(px, py, pz) do { \
            verts[vi++] = (px) - (float)r->lod_tree.world_origin[0]; \
            verts[vi++] = (py) - (float)r->lod_tree.world_origin[1]; \
            verts[vi++] = (pz) - (float)r->lod_tree.world_origin[2]; \
        } while(0)

        // RED arrow: raycast origin (cam->position) + local_up direction
        {
            HMM_Vec3 perp;
            HMM_Vec3 ref = {{1, 0, 0}};
            if (fabsf(HMM_DotV3(cam->local_up, ref)) > 0.9f) ref = (HMM_Vec3){{0, 1, 0}};
            perp = HMM_NormV3(HMM_Cross(cam->local_up, ref));

            float bx = cam->position.X, by = cam->position.Y, bz = cam->position.Z;
            float tx = bx + cam->local_up.X * arrow_len;
            float ty = by + cam->local_up.Y * arrow_len;
            float tz = bz + cam->local_up.Z * arrow_len;
            // Shaft
            ARROW_V(bx, by, bz); ARROW_V(tx, ty, tz);
            // Arrowhead
            float hx = tx - cam->local_up.X * head_len;
            float hy = ty - cam->local_up.Y * head_len;
            float hz = tz - cam->local_up.Z * head_len;
            ARROW_V(tx, ty, tz);
            ARROW_V(hx + perp.X * head_spread, hy + perp.Y * head_spread, hz + perp.Z * head_spread);
            ARROW_V(tx, ty, tz);
            ARROW_V(hx - perp.X * head_spread, hy - perp.Y * head_spread, hz - perp.Z * head_spread);
        }

        // GREEN arrow: camera center + forward direction (50m)
        {
            HMM_Vec3 perp2;
            HMM_Vec3 ref2 = {{1, 0, 0}};
            if (fabsf(HMM_DotV3(cam->forward, ref2)) > 0.9f) ref2 = (HMM_Vec3){{0, 1, 0}};
            perp2 = HMM_NormV3(HMM_Cross(cam->forward, ref2));

            float bx = cam->position.X, by = cam->position.Y, bz = cam->position.Z;
            float tx = bx + cam->forward.X * arrow_len;
            float ty = by + cam->forward.Y * arrow_len;
            float tz = bz + cam->forward.Z * arrow_len;
            // Shaft
            ARROW_V(bx, by, bz); ARROW_V(tx, ty, tz);
            // Arrowhead
            float hx = tx - cam->forward.X * head_len;
            float hy = ty - cam->forward.Y * head_len;
            float hz = tz - cam->forward.Z * head_len;
            ARROW_V(tx, ty, tz);
            ARROW_V(hx + perp2.X * head_spread, hy + perp2.Y * head_spread, hz + perp2.Z * head_spread);
            ARROW_V(tx, ty, tz);
            ARROW_V(hx - perp2.X * head_spread, hy - perp2.Y * head_spread, hz - perp2.Z * head_spread);
        }
        #undef ARROW_V

        sg_update_buffer(r->normal_debug_buf, &(sg_range){ verts, sizeof(verts) });

        highlight_vs_params_t dbg_vs = {
            .mvp = vp_terrain,
            .camera_offset = vs_params.camera_offset,
            .camera_offset_low = vs_params.camera_offset_low,
            .log_depth = vs_params.log_depth,
        };
        sg_bindings dbg_bind = {0};
        dbg_bind.vertex_buffers[0] = r->normal_debug_buf;

        // Draw RED arrow (raycast origin + up)
        sg_apply_pipeline(r->highlight_face_pip);
        sg_apply_uniforms(UB_highlight_vs_params, &SG_RANGE(dbg_vs));
        highlight_fs_params_t red_fs = { .color = (HMM_Vec4){{1.0f, 0.0f, 0.0f, 1.0f}} };
        sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(red_fs));
        sg_apply_bindings(&dbg_bind);
        sg_draw(0, 6, 1);

        // Draw GREEN arrow (camera forward)
        highlight_fs_params_t green_fs = { .color = (HMM_Vec4){{0.0f, 1.0f, 0.0f, 1.0f}} };
        sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(green_fs));
        sg_draw(6, 6, 1);
    }

    // ---- Draw AI agent solid mesh (planet pipeline, TRIANGLES) ----
    if (r->agent_system) {
        AiAgentSystem* asys = (AiAgentSystem*)r->agent_system;

        static sg_buffer s_mesh_buf = {0};
        if (s_mesh_buf.id == SG_INVALID_ID) {
            s_mesh_buf = sg_make_buffer(&(sg_buffer_desc){
                .size = AI_AGENT_MAX * AI_AGENT_MESH_MAX_VERTS * 36,
                .usage.vertex_buffer = true,
                .usage.stream_update = true,
                .label = "agent-mesh",
            });
        }

        sg_apply_pipeline(r->pip);
        sg_apply_uniforms(UB_planet_vs_params, &SG_RANGE(vs_params));
        sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(fs_params));

        for (int ai = 0; ai < asys->count; ai++) {
            AiAgent* agent = &asys->agents[ai];
            if (!agent->active || !agent->mesh_valid || !agent->mesh_verts ||
                agent->state == AGENT_STATE_SLEEPING) continue;

            float px = (float)agent->pos_d[0];
            float py = (float)agent->pos_d[1];
            float pz = (float)agent->pos_d[2];
            float wo0 = (float)r->lod_tree.world_origin[0];
            float wo1 = (float)r->lod_tree.world_origin[1];
            float wo2 = (float)r->lod_tree.world_origin[2];

            HMM_Vec3 aup = agent->local_up;
            HMM_Vec3 afwd = agent->forward;
            HMM_Vec3 art = HMM_NormV3(HMM_Cross(afwd, aup));
            afwd = HMM_NormV3(HMM_Cross(aup, art));

            // Transform model-space verts to floating-origin space
            // Same coordinate space as terrain: (world_pos - world_origin)
            int nv = agent->vert_count_total;
            static float xf[AI_AGENT_MESH_MAX_VERTS * 9];
            float* src = agent->mesh_verts;
            for (int v = 0; v < nv; v++) {
                float mx = src[v*9], my = src[v*9+1], mz = src[v*9+2];
                float mnx = src[v*9+3], mny = src[v*9+4], mnz = src[v*9+5];
                // Rotate + translate to world, subtract world_origin
                xf[v*9+0] = px + art.X*mx + aup.X*my + afwd.X*mz - wo0;
                xf[v*9+1] = py + art.Y*mx + aup.Y*my + afwd.Y*mz - wo1;
                xf[v*9+2] = pz + art.Z*mx + aup.Z*my + afwd.Z*mz - wo2;
                xf[v*9+3] = art.X*mnx + aup.X*mny + afwd.X*mnz;
                xf[v*9+4] = art.Y*mnx + aup.Y*mny + afwd.Y*mnz;
                xf[v*9+5] = art.Z*mnx + aup.Z*mny + afwd.Z*mnz;
                xf[v*9+6] = src[v*9+6];
                xf[v*9+7] = src[v*9+7];
                xf[v*9+8] = src[v*9+8];
            }

            int off = sg_append_buffer(s_mesh_buf, &(sg_range){xf, nv * 36});
            if (!sg_query_buffer_overflow(s_mesh_buf)) {
                sg_bindings mb = {0};
                mb.vertex_buffers[0] = s_mesh_buf;
                mb.vertex_buffer_offsets[0] = off;
                sg_apply_bindings(&mb);
                sg_draw(0, nv, 1);
            }
        }
    }

    // ---- Draw AI agents wireframe ----
    // 3 parts per agent: (1) magenta up-arrow, (2) bounding box, (3) wireframe body
    // All use highlight_face_pip (LINES with depth test + log depth).
    // Verts in floating-origin space (world_pos - world_origin), same as debug arrows.
    if (r->agent_system) {
        AiAgentSystem* asys = (AiAgentSystem*)r->agent_system;

        // Dedicated agent line buffer — 256 verts max (created once)
        static sg_buffer s_abuf = {0};
        if (s_abuf.id == SG_INVALID_ID) {
            s_abuf = sg_make_buffer(&(sg_buffer_desc){
                .size = 256 * 3 * sizeof(float),
                .usage.vertex_buffer = true,
                .usage.stream_update = true,
                .label = "agent-lines",
            });
        }

        highlight_vs_params_t avs = {
            .mvp = vp_terrain,
            .camera_offset = vs_params.camera_offset,
            .camera_offset_low = vs_params.camera_offset_low,
            .log_depth = vs_params.log_depth,
        };

        for (int ai = 0; ai < asys->count; ai++) {
            AiAgent* agent = &asys->agents[ai];
            if (!agent->active || agent->state == AGENT_STATE_SLEEPING) continue;

            float px = (float)agent->pos_d[0];
            float py = (float)agent->pos_d[1];
            float pz = (float)agent->pos_d[2];
            float wo0 = (float)r->lod_tree.world_origin[0];
            float wo1 = (float)r->lod_tree.world_origin[1];
            float wo2 = (float)r->lod_tree.world_origin[2];

            HMM_Vec3 up = agent->local_up;
            HMM_Vec3 fwd = agent->forward;
            HMM_Vec3 rt = HMM_NormV3(HMM_Cross(fwd, up));
            fwd = HMM_NormV3(HMM_Cross(up, rt));

            // Vertex helper: local (sx,sy,sz) -> floating-origin space
            #define V(buf,i, sx,sy,sz) do { \
                (buf)[(i)*3  ] = px+rt.X*(sx)+up.X*(sy)+fwd.X*(sz)-wo0; \
                (buf)[(i)*3+1] = py+rt.Y*(sx)+up.Y*(sy)+fwd.Y*(sz)-wo1; \
                (buf)[(i)*3+2] = pz+rt.Z*(sx)+up.Z*(sy)+fwd.Z*(sz)-wo2; \
            } while(0)

            float verts[256*3];
            int n = 0;
            float h = 1.5f;

            // ---- (1) Magenta up-arrow (20m tall) ----
            V(verts,n, 0,0,0); n++; V(verts,n, 0,20,0); n++; // shaft
            V(verts,n, 0,20,0); n++; V(verts,n, 2,15,0); n++; // head right
            V(verts,n, 0,20,0); n++; V(verts,n, -2,15,0); n++; // head left
            int arrow_count = n;

            // ---- (2) Bounding box (orange) ----
            float bw=0.3f, bd=0.3f;
            #define BOX_EDGE(x1,y1,z1,x2,y2,z2) V(verts,n,x1,y1,z1);n++;V(verts,n,x2,y2,z2);n++
            // bottom
            BOX_EDGE(-bw,0,-bd, bw,0,-bd); BOX_EDGE(-bw,0,bd, bw,0,bd);
            BOX_EDGE(-bw,0,-bd, -bw,0,bd); BOX_EDGE(bw,0,-bd, bw,0,bd);
            // top
            BOX_EDGE(-bw,h,-bd, bw,h,-bd); BOX_EDGE(-bw,h,bd, bw,h,bd);
            BOX_EDGE(-bw,h,-bd, -bw,h,bd); BOX_EDGE(bw,h,-bd, bw,h,bd);
            // verticals
            BOX_EDGE(-bw,0,-bd, -bw,h,-bd); BOX_EDGE(bw,0,-bd, bw,h,-bd);
            BOX_EDGE(-bw,0,bd, -bw,h,bd); BOX_EDGE(bw,0,bd, bw,h,bd);
            #undef BOX_EDGE
            int box_count = n - arrow_count;

            // ---- (3) Wireframe body ----
            float tw=0.18f, td2=0.12f;
            float tb=h*0.3f, tt=h*0.7f;
            float hb=h*0.78f;
            float hdw=0.13f, hdd=0.1f;
            float ls=0.07f, as=tw+0.04f;
            int body_start = n;
            // torso verticals
            V(verts,n,-tw,tb,-td2);n++;V(verts,n,-tw,tt,-td2);n++;
            V(verts,n, tw,tb,-td2);n++;V(verts,n, tw,tt,-td2);n++;
            V(verts,n,-tw,tb, td2);n++;V(verts,n,-tw,tt, td2);n++;
            V(verts,n, tw,tb, td2);n++;V(verts,n, tw,tt, td2);n++;
            // torso top ring
            V(verts,n,-tw,tt,-td2);n++;V(verts,n, tw,tt,-td2);n++;
            V(verts,n,-tw,tt, td2);n++;V(verts,n, tw,tt, td2);n++;
            V(verts,n,-tw,tt,-td2);n++;V(verts,n,-tw,tt, td2);n++;
            V(verts,n, tw,tt,-td2);n++;V(verts,n, tw,tt, td2);n++;
            // torso bottom ring
            V(verts,n,-tw,tb,-td2);n++;V(verts,n, tw,tb,-td2);n++;
            V(verts,n,-tw,tb, td2);n++;V(verts,n, tw,tb, td2);n++;
            V(verts,n,-tw,tb,-td2);n++;V(verts,n,-tw,tb, td2);n++;
            V(verts,n, tw,tb,-td2);n++;V(verts,n, tw,tb, td2);n++;
            // head box
            V(verts,n,-hdw,hb,-hdd);n++;V(verts,n, hdw,hb,-hdd);n++;
            V(verts,n,-hdw,hb, hdd);n++;V(verts,n, hdw,hb, hdd);n++;
            V(verts,n,-hdw,hb,-hdd);n++;V(verts,n,-hdw,hb, hdd);n++;
            V(verts,n, hdw,hb,-hdd);n++;V(verts,n, hdw,hb, hdd);n++;
            V(verts,n,-hdw,h,-hdd);n++;V(verts,n, hdw,h,-hdd);n++;
            V(verts,n,-hdw,h, hdd);n++;V(verts,n, hdw,h, hdd);n++;
            V(verts,n,-hdw,h,-hdd);n++;V(verts,n,-hdw,h, hdd);n++;
            V(verts,n, hdw,h,-hdd);n++;V(verts,n, hdw,h, hdd);n++;
            V(verts,n,-hdw,hb,-hdd);n++;V(verts,n,-hdw,h,-hdd);n++;
            V(verts,n, hdw,hb,-hdd);n++;V(verts,n, hdw,h,-hdd);n++;
            V(verts,n,-hdw,hb, hdd);n++;V(verts,n,-hdw,h, hdd);n++;
            V(verts,n, hdw,hb, hdd);n++;V(verts,n, hdw,h, hdd);n++;
            // legs
            V(verts,n,-ls,0,0);n++;V(verts,n,-ls,tb,0);n++;
            V(verts,n, ls,0,0);n++;V(verts,n, ls,tb,0);n++;
            // arms
            V(verts,n,-as,h*0.5f,0);n++;V(verts,n,-as,tt,0);n++;
            V(verts,n, as,h*0.5f,0);n++;V(verts,n, as,tt,0);n++;
            int body_count = n - body_start;

            #undef V

            // Upload all verts at once
            sg_update_buffer(s_abuf, &(sg_range){ verts, n*3*sizeof(float) });

            sg_apply_pipeline(r->highlight_face_pip);
            sg_apply_uniforms(UB_highlight_vs_params, &SG_RANGE(avs));

            sg_bindings ab = {0};
            ab.vertex_buffers[0] = s_abuf;
            sg_apply_bindings(&ab);

            // Draw (1) arrow in magenta
            highlight_fs_params_t c1 = {.color=(HMM_Vec4){{1,0,1,1}}};
            sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(c1));
            sg_draw(0, arrow_count, 1);

            // Draw (2) bounding box in orange
            highlight_fs_params_t c2 = {.color=(HMM_Vec4){{1,0.5f,0,0.8f}}};
            sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(c2));
            sg_draw(arrow_count, box_count, 1);

            // Draw (3) body in white (or green if talking)
            HMM_Vec4 bc = (agent->state==AGENT_STATE_TALKING)
                ? (HMM_Vec4){{0,1,0,1}} : (HMM_Vec4){{1,1,1,1}};
            highlight_fs_params_t c3 = {.color=bc};
            sg_apply_uniforms(UB_highlight_fs_params, &SG_RANGE(c3));
            sg_draw(body_start, body_count, 1);
        }
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
        r->profile.draw_calls = r->profile.lod_patches + 2; // +sky+atmos
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

    // Crosshair: small white filled circle via sgl
    touch_render_crosshair();

    sdtx_canvas(sapp_widthf() * 0.5f, sapp_heightf() * 0.5f);
    sdtx_origin(0.5f, 0.5f);
    sdtx_font(0);

    sdtx_color3f(1.0f, 1.0f, 0.0f);
    sdtx_printf("FPS: %.0f\n", r->display_fps);

    // Player altitude (relative to gravity body surface)
    {
        double bdx = cam->pos_d[0] - r->lod_tree.body_center_d[0];
        double bdy = cam->pos_d[1] - r->lod_tree.body_center_d[1];
        double bdz = cam->pos_d[2] - r->lod_tree.body_center_d[2];
        float cam_r = (float)sqrt(bdx*bdx + bdy*bdy + bdz*bdz);
        float surface_r = r->lod_tree.planet_radius;
        if (r->lod_tree.body_type == LOD_BODY_MOON && cam_r > 1.0f) {
            HMM_Vec3 cam_dir = {{ (float)(bdx/cam_r), (float)(bdy/cam_r), (float)(bdz/cam_r) }};
            surface_r = moon_surface_radius(&r->lod_tree.moon_shape, cam_dir);
        }
        float altitude = cam_r - surface_r;
        float body_r_km = r->lod_tree.planet_radius / 1000.0f;
        sdtx_color3f(0.6f, 0.8f, 1.0f);
        sdtx_printf("alt: %.0fm  R: %.0fkm\n", altitude, body_r_km);
    }

    // Jetpack status
    if (cam->jetpack_active) {
        sdtx_color3f(0.3f, 1.0f, 0.5f);
        sdtx_printf("JETPACK ON  speed: %.0fx\n", cam->jetpack_speed_mult);
    } else {
        sdtx_color3f(0.5f, 0.5f, 0.5f);
        sdtx_puts("JETPACK OFF\n");
    }

    // Selected block type
    {
        static const char* block_names[] = {
            "AIR", "WATER", "SAND", "DIRT", "GRASS", "STONE", "ICE", "BEDROCK", "TORCH"
        };
        // Hotbar types: STONE, DIRT, GRASS, SAND, WATER, ICE, TORCH
        static const int hotbar_voxel[] = { 5, 3, 4, 2, 1, 6, 8 };
        int vtype = hotbar_voxel[r->hotbar_selected_slot % 7];
        const char* name = (vtype < 9) ? block_names[vtype] : "???";
        sdtx_color3f(1.0f, 0.8f, 0.4f);
        sdtx_printf("[%d] %s\n", r->hotbar_selected_slot + 1, name);
    }

    // Reference frame indicator (top-right corner)
    {
        float canvas_w = sapp_widthf() * 0.5f;
        const char* body_name;
        if (cam->gravity_body < 0) {
            body_name = "Tenebris";
        } else if (cam->gravity_body < r->solar_system.moon_count) {
            body_name = r->solar_system.moons[cam->gravity_body].name;
        } else {
            body_name = "???";
        }
        int name_len = 0;
        for (const char* p = body_name; *p; p++) name_len++;
        // Position: right edge minus name length, top row
        float col = canvas_w / 8.0f - (float)name_len - 1.0f;
        sdtx_color3f(0.6f, 0.8f, 1.0f);
        sdtx_pos(col, 0.0f);
        sdtx_puts(body_name);
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

    // ---- 4b. Moon name labels (space mode only) ----
    solar_system_draw_labels(&r->solar_system, cam, vp_rot);

    // Space mode indicator
    if (cam->space_mode) {
        sdtx_color3f(0.3f, 0.8f, 1.0f);
        sdtx_pos(0.5f, 0.5f + 4.0f);
        sdtx_puts("SPACE FLIGHT  Q/E roll  Scroll=speed");
    } else if (cam->gravity_body >= 0 && cam->gravity_body < r->solar_system.moon_count) {
        sdtx_color3f(0.4f, 1.0f, 0.6f);
        sdtx_pos(0.5f, 0.5f + 4.0f);
        sdtx_printf("Moon: %s", r->solar_system.moons[cam->gravity_body].name);
    }

    sdtx_draw();

    // ---- 5. Hotbar UI overlay ----
    {
        // Atlas tile indices for each hotbar slot (matches hotbar_types[] in main.c)
        // STONE=4, DIRT=2, GRASS=3, SAND=1, WATER=0, ICE=5, TORCH=9
        static const int hotbar_atlas[] = { 4, 2, 3, 1, 0, 5, 9 };
        static const int HOTBAR_SLOTS = 7;
        static const float ATLAS_TILES = 11.0f;

        float sw = sapp_widthf();
        float sh = sapp_heightf();
        float slot_size = 40.0f;   // pixels per slot
        float padding = 4.0f;      // gap between slots
        float total_w = HOTBAR_SLOTS * slot_size + (HOTBAR_SLOTS - 1) * padding;
        float x0 = (sw - total_w) * 0.5f;
        float y0 = sh - slot_size - 12.0f;  // 12px from bottom

        // Build 8 quads (6 verts each): pos.xy + uv.xy = 4 floats per vert
        float verts[8 * 6 * 4];
        int vi = 0;
        for (int i = 0; i < HOTBAR_SLOTS; i++) {
            float sx = x0 + i * (slot_size + padding);
            float sy = y0;
            float u0 = (float)hotbar_atlas[i] / ATLAS_TILES;
            float u1 = ((float)hotbar_atlas[i] + 1.0f) / ATLAS_TILES;

            // Triangle 1: top-left, top-right, bottom-right
            verts[vi++] = sx;             verts[vi++] = sy;               verts[vi++] = u0; verts[vi++] = 0.0f;
            verts[vi++] = sx + slot_size; verts[vi++] = sy;               verts[vi++] = u1; verts[vi++] = 0.0f;
            verts[vi++] = sx + slot_size; verts[vi++] = sy + slot_size;   verts[vi++] = u1; verts[vi++] = 1.0f;
            // Triangle 2: top-left, bottom-right, bottom-left
            verts[vi++] = sx;             verts[vi++] = sy;               verts[vi++] = u0; verts[vi++] = 0.0f;
            verts[vi++] = sx + slot_size; verts[vi++] = sy + slot_size;   verts[vi++] = u1; verts[vi++] = 1.0f;
            verts[vi++] = sx;             verts[vi++] = sy + slot_size;   verts[vi++] = u0; verts[vi++] = 1.0f;
        }
        sg_update_buffer(r->hotbar_buf, &(sg_range){ verts, sizeof(verts) });

        // Orthographic projection for screen-space rendering
        HMM_Mat4 ortho = HMM_Orthographic_LH_ZO(0.0f, sw, sh, 0.0f, -1.0f, 1.0f);

        // Draw each slot (all as one batch, then highlight selected)
        sg_apply_pipeline(r->hotbar_pip);

        // Draw all slots with dim tint first
        sg_bindings hotbar_bind = {
            .vertex_buffers[0] = r->hotbar_buf,
            .views[VIEW_hotbar_tex] = r->hex_atlas_view,
            .samplers[SMP_hotbar_smp] = r->hex_atlas_smp,
        };
        sg_apply_bindings(&hotbar_bind);

        hotbar_vs_params_t vs_p = { .ortho = ortho };
        sg_apply_uniforms(UB_hotbar_vs_params, &SG_RANGE(vs_p));

        // Draw unselected slots with semi-transparent tint
        for (int i = 0; i < HOTBAR_SLOTS; i++) {
            float brightness = (i == r->hotbar_selected_slot) ? 1.0f : 0.5f;
            float alpha = (i == r->hotbar_selected_slot) ? 1.0f : 0.7f;
            hotbar_fs_params_t fs_p = { .tint = {{ brightness, brightness, brightness, alpha }} };
            sg_apply_uniforms(UB_hotbar_fs_params, &SG_RANGE(fs_p));
            sg_draw(i * 6, 6, 1);
        }
    }

    // Touch controls overlay (joysticks + buttons)
    if (r->touch) {
        touch_render(r->touch, r->hotbar_selected_slot, 7);
    }

    sg_end_pass();
    sg_commit();
}

void render_shutdown(Renderer* r) {
    solar_system_shutdown(&r->solar_system);
    hex_terrain_destroy(&r->hex_terrain);
    atmosphere_destroy(&r->atmosphere);
    lod_tree_destroy(&r->lod_tree);
    torch_destroy(&r->torch_system);
    if (r->wireframe_idx.id != SG_INVALID_ID) {
        sg_destroy_buffer(r->wireframe_idx);
    }
    if (r->lod_wireframe_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->lod_wireframe_pip);
    }
    if (r->hex_wireframe_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->hex_wireframe_pip);
    }
    if (r->highlight_buf.id != SG_INVALID_ID) {
        sg_destroy_buffer(r->highlight_buf);
    }
    if (r->highlight_face_buf.id != SG_INVALID_ID) {
        sg_destroy_buffer(r->highlight_face_buf);
    }
    if (r->highlight_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->highlight_pip);
    }
    if (r->highlight_face_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->highlight_face_pip);
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
    if (r->planet_tex_view.id != SG_INVALID_ID) {
        sg_destroy_view(r->planet_tex_view);
    }
    if (r->planet_tex_img.id != SG_INVALID_ID) {
        sg_destroy_image(r->planet_tex_img);
    }
    if (r->planet_tex_smp.id != SG_INVALID_ID) {
        sg_destroy_sampler(r->planet_tex_smp);
    }
    if (r->hotbar_pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(r->hotbar_pip);
    }
    if (r->hotbar_buf.id != SG_INVALID_ID) {
        sg_destroy_buffer(r->hotbar_buf);
    }
}

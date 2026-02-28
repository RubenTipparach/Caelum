#ifndef RENDER_H
#define RENDER_H

#include "HandmadeMath.h"
#include "sokol_gfx.h"
#include "planet.h"
#include "camera.h"
#include "atmosphere.h"
#include "lod.h"
#include "hex_terrain.h"

// ---- Profiler stats (F3 overlay) ----
typedef struct ProfileStats {
    // Per-frame CPU timings (ms), smoothed over ~0.5s
    float camera_ms;
    float lod_update_ms;
    float lod_upload_ms;
    float hex_update_ms;
    float hex_upload_ms;
    float render_ms;
    float frame_total_ms;

    // GPU resource counts
    int draw_calls;
    int total_triangles;

    // Memory estimates
    int hex_chunks_active;
    int hex_chunks_meshed;
    int hex_gpu_verts;
    int lod_gpu_verts;
    int lod_patches;
    int lod_nodes;
    int job_pending;

    // Smoothing accumulator
    float accum_time;
    int accum_frames;
    float accum_camera_ms;
    float accum_lod_update_ms;
    float accum_lod_upload_ms;
    float accum_hex_update_ms;
    float accum_hex_upload_ms;
    float accum_render_ms;
    float accum_frame_ms;
} ProfileStats;

typedef struct Renderer {
    // Planet rendering
    sg_pipeline pip;
    sg_bindings bind;
    int vertex_count;

    // Sky rendering
    sg_pipeline sky_pip;
    sg_bindings sky_bind;

    // Atmosphere
    Atmosphere atmosphere;

    // LOD tree (large-scale planet rendering)
    LodTree lod_tree;
    bool show_lod_debug;    // L key: color patches by LOD depth + show depth stats

    // Hex terrain (close-range voxel grid)
    HexTerrain hex_terrain;

    // Common
    sg_pass_action pass_action;
    HMM_Vec3 sun_direction;
    float sun_angle;

    // FPS counter
    float fps_accumulator;
    int fps_frame_count;
    float display_fps;

    // Profiler (F3)
    bool show_profiler;
    ProfileStats profile;
} Renderer;

void render_init(Renderer* r, Planet* planet, const Camera* cam);
void render_update_mesh(Renderer* r, Planet* planet, const Camera* cam);
void render_frame(Renderer* r, const Camera* cam, float dt);
void render_shutdown(Renderer* r);

#endif

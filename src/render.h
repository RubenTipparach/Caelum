#ifndef RENDER_H
#define RENDER_H

#include "HandmadeMath.h"
#include "sokol_gfx.h"
#include "planet.h"
#include "camera.h"
#include "atmosphere.h"
#include "lod.h"
#include "hex_terrain.h"

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
} Renderer;

void render_init(Renderer* r, Planet* planet, const Camera* cam);
void render_update_mesh(Renderer* r, Planet* planet, const Camera* cam);
void render_frame(Renderer* r, const Camera* cam, float dt);
void render_shutdown(Renderer* r);

#endif

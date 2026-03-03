// torch.h -- Torch rendering system API
// Manages torch model instances on the planet surface with animated flames.

#ifndef TORCH_H
#define TORCH_H

#include "HandmadeMath.h"
#include "sokol_gfx.h"
#include <stdbool.h>
#include <stdint.h>

#define TORCH_MAX_INSTANCES 256

typedef struct TorchInstance {
    bool active;
    float world_pos[3];  // Double-precision world position stored as doubles
    float up_dir[3];     // Surface normal (local up direction)
    int chunk_cx, chunk_cz;  // Owning chunk coords (for chunk-based removal)
    int gcol, grow, layer;   // Grid coords (for individual removal)
    float time;              // Animation time accumulator
} TorchInstance;

typedef struct TorchSystem {
    sg_pipeline pip;
    sg_buffer model_buf;      // Static torch geometry VBO
    int model_vertex_count;   // Number of vertices in torch model

    TorchInstance instances[TORCH_MAX_INSTANCES];
    int instance_count;
} TorchSystem;

void torch_init(TorchSystem* ts);
void torch_destroy(TorchSystem* ts);
void torch_update(TorchSystem* ts, float dt);

// Add a torch instance. Returns index or -1 if full.
int torch_add(TorchSystem* ts, float world_x, float world_y, float world_z,
              float up_x, float up_y, float up_z,
              int chunk_cx, int chunk_cz, int gcol, int grow, int layer);

// Remove torch at specific grid position. Returns true if found.
bool torch_remove(TorchSystem* ts, int gcol, int grow, int layer);

// Remove all torches belonging to a chunk.
void torch_remove_chunk(TorchSystem* ts, int chunk_cx, int chunk_cz);

// Render all active torch instances.
void torch_render(TorchSystem* ts, HMM_Mat4 vp,
                  const double camera_offset[3], const double camera_offset_low[3],
                  float fcoef, float far_plane, float z_bias,
                  HMM_Vec3 sun_dir);

#endif // TORCH_H

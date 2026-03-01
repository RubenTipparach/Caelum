#ifndef HEX_TERRAIN_H
#define HEX_TERRAIN_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "sokol_gfx.h"
#include "hex_vertex.h"
#include "job_system.h"

// ---- Configuration ----
#define HEX_RADIUS          1.0f        // Circumradius of each hex (center to vertex) in meters
#define HEX_HEIGHT          1.0f        // Height of one voxel layer in meters
#define HEX_RANGE           500.0f      // Activation range from camera in meters
#define HEX_INNER_RANGE     350.0f      // Full voxel prisms inside this range
#define HEX_TRANSITION_ON   375.0f      // Switch chunk TO transition mode (hysteresis high)
#define HEX_TRANSITION_OFF  325.0f      // Revert chunk FROM transition mode (hysteresis low)
#define HEX_CHUNK_SIZE      32          // Hex columns per chunk side (32x32 grid)
#define HEX_MAX_CHUNKS      512         // Maximum active chunks
#define HEX_MAX_COLUMN_H    16384       // Maximum terrain height in layers (int16_t, covers up to 16km)
#define HEX_MAX_UPLOADS      32         // GPU uploads per frame
#define HEX_MAX_ACTIVATIONS  16         // Max new chunks activated per frame (spread load)
#define HEX_MAX_DRAW_ALT    2000.0f     // Hex terrain disabled above this altitude (m above ground)
#define HEX_SURFACE_BIAS     0.05f      // Radial offset (m) to prevent Z-fighting with smooth LOD
#define HEX_COL_SPACING     (1.5f * HEX_RADIUS)               // Horizontal center-to-center
#define HEX_ROW_SPACING     (1.7320508f * HEX_RADIUS)         // sqrt(3) * HEX_RADIUS

// Chunk world dimensions
#define HEX_CHUNK_WIDTH     (HEX_CHUNK_SIZE * HEX_COL_SPACING)    // ~48m
#define HEX_CHUNK_DEPTH     (HEX_CHUNK_SIZE * HEX_ROW_SPACING)    // ~55m

// ---- Hex chunk ----
typedef struct HexChunk {
    int cx, cz;                 // Chunk grid coordinates
    bool active;                // Currently in use
    bool dirty;                 // Needs mesh regeneration
    bool generating;            // Mesh is being generated on worker thread
    bool is_transition;         // True = smooth surface-only mesh (transition zone)

    // Per-column terrain data (offset grid within chunk)
    int16_t heights[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE];   // Terrain height in layers
    uint8_t types[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE];     // VoxelType of surface block

    // GPU mesh
    sg_buffer gpu_buffer;
    int gpu_vertex_count;

    // CPU mesh (generated, then uploaded)
    HexVertex* cpu_vertices;
    int cpu_vertex_count;

    void* pending_job;          // HexMeshJob* (opaque)
} HexChunk;

// ---- Hex terrain system ----
typedef struct HexTerrain {
    HexChunk chunks[HEX_MAX_CHUNKS];
    int active_count;

    // Planet parameters
    float planet_radius;
    float layer_thickness;
    int sea_level;
    int seed;

    // Tangent frame (computed from camera ground-projected position)
    HMM_Vec3 tangent_origin;   // World position of grid origin on surface
    HMM_Vec3 tangent_up;       // Local up (radial outward)
    HMM_Vec3 tangent_east;     // Tangent plane X axis
    HMM_Vec3 tangent_north;    // Tangent plane Z axis

    // Stable frame: tangent frame is locked and only re-anchored when camera
    // drifts far from the frame origin. This prevents constant chunk regen.
    bool frame_valid;
    HMM_Vec3 frame_anchor;     // Camera pos when frame was last anchored

    // Floating origin (shared with LOD tree)
    double world_origin[3];

    // Camera position (float, relative to floating origin)
    HMM_Vec3 camera_pos;

    // Job system (shared with LOD tree)
    JobSystem* jobs;

    // Stats
    int total_vertex_count;
    int chunks_rendered;
} HexTerrain;

// ---- Hex selection / interaction ----
typedef struct HexHitResult {
    bool valid;              // True if a hex was hit
    int chunk_index;         // Index into ht->chunks[]
    int col, row;            // Local column/row within chunk
    int gcol, grow;          // Global grid coordinates
    int height;              // Terrain height at hit hex (in layers)
    uint8_t type;            // VoxelType of surface block
} HexHitResult;

// ---- API ----

void hex_terrain_init(HexTerrain* ht, float planet_radius, float layer_thickness,
                      int sea_level, int seed, JobSystem* jobs);

void hex_terrain_destroy(HexTerrain* ht);

// Update: load/unload chunks based on camera position, generate meshes
void hex_terrain_update(HexTerrain* ht, HMM_Vec3 camera_pos,
                        const double world_origin[3]);

// Upload CPU meshes to GPU
void hex_terrain_upload_meshes(HexTerrain* ht);

// Render all active hex chunks (pipeline must be applied by caller)
void hex_terrain_render(HexTerrain* ht, sg_pipeline pip,
                        sg_view atlas_view, sg_sampler atlas_smp);

// Check if a world position is within hex terrain range (for LOD suppression)
bool hex_terrain_covers_position(const HexTerrain* ht, HMM_Vec3 world_pos);

// Get the activation range (for LOD suppression)
float hex_terrain_get_range(void);

// Get effective suppress range based on current mesh coverage.
// Returns HEX_RANGE when >=80% of active chunks have GPU meshes, 0 otherwise.
// Use to prevent black holes during chunk generation.
float hex_terrain_effective_range(const HexTerrain* ht);

// Raycast from camera into hex terrain. Returns hit info.
HexHitResult hex_terrain_raycast(const HexTerrain* ht, HMM_Vec3 ray_origin,
                                  HMM_Vec3 ray_dir, float max_dist);

// Remove top layer at the given hex. Marks chunk dirty for remeshing.
bool hex_terrain_break(HexTerrain* ht, const HexHitResult* hit);

// Place a block on top of the given hex. Marks chunk dirty for remeshing.
bool hex_terrain_place(HexTerrain* ht, const HexHitResult* hit, uint8_t voxel_type);

// Generate wireframe vertices for a hex selection highlight.
// Writes 12 float3 vertices (6 line segments) into out_verts (must hold 36 floats).
// Returns true if successful.
bool hex_terrain_build_highlight(const HexTerrain* ht, const HexHitResult* hit,
                                  const double world_origin[3], float* out_verts);

#endif

#ifndef LOD_H
#define LOD_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "sokol_gfx.h"
#include "job_system.h"

// ---- Configuration ----
#define LOD_MAX_DEPTH       13      // Max subdivision depth
#define LOD_MAX_NODES       16384   // Node pool (intermediate + leaf nodes)
#define LOD_ROOT_COUNT      20      // 20 icosahedron faces
#define LOD_CHILDREN        4       // Aperture-4 subdivision
#define LOD_MAX_UPLOADS_PER_FRAME 64 // Max GPU uploads per frame
#define LOD_MAX_SPLITS_PER_FRAME 256 // Tree builds fast; capped at voxel depth anyway
#define LOD_NUM_WORKERS     4       // Worker threads for mesh generation
#define LOD_SPLIT_FACTOR    8.0f    // Split when distance < patch_arc * factor (higher = more detail)

// Vertex format for LOD mesh patches (same as existing Vertex)
typedef struct LodVertex {
    float pos[3];
    float normal[3];
    float color[3];
} LodVertex;

// ---- Spherical triangle on the unit sphere ----
typedef struct SphericalTriangle {
    HMM_Vec3 v0, v1, v2;       // Vertices on unit sphere
    HMM_Vec3 center;            // Centroid (normalized to unit sphere)
    float angular_radius;        // Half-angle subtended from center
} SphericalTriangle;

// ---- LOD node states ----
typedef enum {
    LOD_UNLOADED,       // No mesh data
    LOD_GENERATING,     // Mesh is being generated (future: on worker thread)
    LOD_READY,          // CPU mesh ready, needs GPU upload
    LOD_ACTIVE,         // GPU buffer uploaded, ready to render
} LodNodeState;

// ---- Single LOD tree node ----
typedef struct LodNode {
    SphericalTriangle tri;
    int depth;
    int parent;                 // Index in nodes array, -1 for roots
    int children[LOD_CHILDREN]; // -1 = not subdivided
    LodNodeState state;

    // GPU mesh data
    sg_buffer gpu_buffer;
    int gpu_vertex_count;

    // CPU mesh data (generated, then uploaded to GPU)
    LodVertex* cpu_vertices;
    int cpu_vertex_count;

    bool is_leaf;               // True if this node has no children
    bool wants_split;           // Marked for subdivision
    bool wants_merge;           // Marked for merge
    void* pending_job;          // MeshGenJob* (opaque, defined in lod.c)
} LodNode;

// ---- The LOD tree ----
typedef struct LodTree {
    LodNode* nodes;
    int node_count;
    int node_capacity;

    int root_nodes[LOD_ROOT_COUNT];

    // Planet parameters
    float planet_radius;        // Base sphere radius (e.g., 796000 for 1/8 Earth)
    float layer_thickness;
    int sea_level;
    int seed;

    // Floating origin (double precision) — all mesh vertex positions are stored
    // relative to this origin to keep float32 vertex coords near zero.
    // When camera moves far from origin, recenter and regenerate meshes.
    double world_origin[3];

    // Camera state (updated each frame)
    HMM_Vec3 camera_pos;

    // Job system for threaded mesh generation
    JobSystem* jobs;

    // Per-level stats (collected each frame in lod_tree_update)
    struct {
        int patch_count;        // Number of visible leaf patches at this depth
        int vertex_count;       // Total vertices rendered at this depth
        float min_distance;     // Closest patch distance from camera
        float max_distance;     // Farthest patch distance from camera
    } level_stats[LOD_MAX_DEPTH + 1];
    int stats_frame_counter;    // For periodic console printing

    // Depth-based split thresholds (precomputed for consistent LOD ordering)
    float depth_arc[LOD_MAX_DEPTH + 1];  // Canonical arc length at each depth

    // Per-frame budget
    int splits_this_frame;

    // Hex terrain suppression: don't render LOD patches closer than this distance.
    // Set to 0 to disable suppression. When > 0, leaf patches with
    // patch_distance < suppress_range are skipped during rendering.
    float suppress_range;

    // Stats
    int active_leaf_count;
    int total_vertex_count;
} LodTree;

// ---- API ----

// Initialize the LOD tree with 20 root icosahedron faces
void lod_tree_init(LodTree* tree, float planet_radius, float layer_thickness,
                   int sea_level, int seed);

// Free all resources
void lod_tree_destroy(LodTree* tree);

// Update the tree: split/merge nodes based on camera distance
// Call once per frame before rendering
void lod_tree_update(LodTree* tree, HMM_Vec3 camera_pos, HMM_Mat4 view_proj);

// Upload CPU meshes to GPU (limited per frame to avoid stalls)
void lod_tree_upload_meshes(LodTree* tree);

// Callback called before each patch draw, receives node depth and user data.
// Use to set per-patch uniforms (e.g., LOD debug color).
// NULL = no callback.
typedef void (*LodPreDrawFunc)(int depth, void* user_data);

// Render all active leaf nodes
void lod_tree_render(LodTree* tree, sg_pipeline pip, HMM_Mat4 vp,
                     LodPreDrawFunc pre_draw, void* user_data);

// Find the finest leaf node containing a world position

// Get terrain height at a world position (for collision)
// Returns double to avoid float quantization at large radii (6cm steps at 800km)
double lod_tree_terrain_height(const LodTree* tree, HMM_Vec3 world_pos);

// Check and recenter floating origin if camera has moved far from current origin.
// Call once per frame with the camera's double-precision position.
// Returns true if origin was recentered (all meshes will regenerate).
bool lod_tree_update_origin(LodTree* tree, const double cam_pos_d[3]);

// Get stats
int lod_tree_active_leaves(const LodTree* tree);

#endif

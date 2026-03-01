#include "lod.h"
#include "math_utils.h"
#include "log_config.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

// FastNoiseLite: implementation is in planet.c (FNL_IMPL defined there).
// We only need the header for fnl_state/fnlCreateState/fnlGetNoise3D.
#include "FastNoiseLite.h"
#include "planet.h"  // VoxelType, MAX_VOXEL_HEIGHT

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---- Forward declarations ----
static void free_node_mesh(LodNode* node);
static int allocate_node(LodTree* tree);
static void init_spherical_triangle(SphericalTriangle* tri,
                                     HMM_Vec3 v0, HMM_Vec3 v1, HMM_Vec3 v2);

// ---- Mesh generation job (for threaded generation) ----
typedef struct MeshGenJob {
    // Inputs (copied by value - safe from node array realloc)
    SphericalTriangle tri;
    int depth;
    float planet_radius;
    float layer_thickness;
    int sea_level;
    int seed;
    float origin[3];  // Floating origin — subtract from world positions

    // Shared hex tangent frame (copied from LodTree for thread safety)
    HMM_Vec3 hex_frame_origin;
    HMM_Vec3 hex_frame_east;
    HMM_Vec3 hex_frame_north;

    // Output (written by worker thread)
    void* vertices;
    int vertex_count;
    int vertex_stride;      // sizeof(LodVertex) or sizeof(HexVertex)
    bool is_hex_mesh;       // True if vertices are HexVertex[]

    // Linkage
    int node_index;
    volatile int completed;  // 0 = pending, 1 = done
} MeshGenJob;

// ---- Helpers ----

// Smooth hermite interpolation between edge0 and edge1
static float smoothstepf(float edge0, float edge1, float x) {
    float t = (x - edge0) / (edge1 - edge0);
    t = fminf(1.0f, fmaxf(0.0f, t));
    return t * t * (3.0f - 2.0f * t);
}

// ---- Terrain noise ----

// Continental landmass noise: very low frequency for large land/ocean shapes
static fnl_state create_continental_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 3;
    noise.frequency = 0.6f;   // Large-scale (effective ~1.8 on unit sphere w/ scale=3)
    noise.seed = seed;
    return noise;
}

// Mountain chain noise: ridged multifractal for sharp peaks + V-shaped valleys
static fnl_state create_mountain_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_RIDGED;
    noise.octaves = 5;
    noise.frequency = 1.5f;   // Mountain-chain scale
    noise.seed = seed + 4000;
    return noise;
}

// Domain warping noise: displaces sample coordinates for organic coastlines
static fnl_state create_warp_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 3;
    noise.frequency = 4.0f;
    noise.seed = seed + 1000;
    return noise;
}

// Detail noise: high-frequency local variation (ridges, bumps)
static fnl_state create_detail_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_RIDGED;
    noise.octaves = 3;
    noise.frequency = 16.0f;
    noise.seed = seed + 2000;
    return noise;
}

// Color variation noise: subtle per-vertex color offsets for visual texture
static fnl_state create_color_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 2;
    noise.frequency = 40.0f;
    noise.seed = seed + 3000;
    return noise;
}

// ---- Terrain height constants ----
// For a 796km radius planet, ~8km total relief (~1% of radius, KSP-like exaggeration)
#define TERRAIN_SEA_LEVEL_M   4000.0f   // Sea level = 4km above planet_radius
#define TERRAIN_AMPLITUDE_M   8000.0f   // Total relief: 0 to 8km above planet_radius
#define TERRAIN_MIN_M         500.0f    // Minimum terrain (deep ocean floor)

// Sample terrain noise value at a point on the unit sphere.
// Multi-layer approach: continental base + ridged mountains + local detail.
// Returns combined noise in roughly [-1, 1] range.
static float sample_terrain_noise(fnl_state* continental, fnl_state* mountain,
                                   fnl_state* warp_noise, fnl_state* detail_noise,
                                   HMM_Vec3 unit_pos) {
    float scale = 3.0f;
    float px = unit_pos.X * scale;
    float py = unit_pos.Y * scale;
    float pz = unit_pos.Z * scale;

    // Domain warping: organic coastlines and terrain shapes
    float warp_strength = 0.5f;
    float wx = fnlGetNoise3D(warp_noise, px + 5.2f, py + 1.3f, pz + 3.7f);
    float wy = fnlGetNoise3D(warp_noise, px + 9.1f, py + 4.8f, pz + 7.2f);
    float wz = fnlGetNoise3D(warp_noise, px + 2.6f, py + 8.4f, pz + 0.9f);
    float wpx = px + wx * warp_strength;
    float wpy = py + wy * warp_strength;
    float wpz = pz + wz * warp_strength;

    // Layer 1: Continental structure — large landmasses vs ocean basins
    float continent = fnlGetNoise3D(continental, wpx, wpy, wpz);

    // Layer 2: Mountain chains — ridged multifractal, masked to land only
    float mountain_raw = fnlGetNoise3D(mountain, wpx, wpy, wpz);
    float mountain_val = (mountain_raw + 1.0f) * 0.5f;  // Remap [-1,1] to [0,1]
    mountain_val *= mountain_val;  // Sharpen: emphasize peaks, flatten valleys
    // Mountains only on solid land (ramp up towards continental interior)
    float land_factor = smoothstepf(-0.05f, 0.35f, continent);
    float mountain_height = mountain_val * land_factor;

    // Layer 3: Local detail — high-freq texture on both land and ocean floor
    float detail = fnlGetNoise3D(detail_noise, px, py, pz);
    float detail_weight = 0.05f + land_factor * 0.10f;  // More detail on land

    // Combine layers (ocean bias shifts land/ocean boundary for ~70% ocean coverage)
    float height = continent * 0.55f - 0.22f;  // Continental base shape
    height += mountain_height * 0.45f;         // Mountain chains on land
    height += detail * detail_weight;           // Local variation

    // Clamp to [-1, 1]
    if (height > 1.0f) height = 1.0f;
    if (height < -1.0f) height = -1.0f;

    return height;
}

// Sample terrain height in METERS above planet_radius.
// Returns continuous floating-point height for coarse LOD rendering and collision.
static float sample_terrain_height_m(fnl_state* continental, fnl_state* mountain,
                                      fnl_state* warp_noise, fnl_state* detail_noise,
                                      HMM_Vec3 unit_pos) {
    float n = sample_terrain_noise(continental, mountain, warp_noise, detail_noise, unit_pos);
    // Map noise [-1,1] to height [TERRAIN_MIN_M, TERRAIN_MIN_M + TERRAIN_AMPLITUDE_M]
    float height = TERRAIN_MIN_M + (n + 1.0f) * 0.5f * TERRAIN_AMPLITUDE_M;
    if (height < 0.0f) height = 0.0f;
    return height;
}

// Legacy: sample terrain height in discrete layers (for voxel mesh at depth 16).
// Returns integer layer index compatible with planet.c voxel system.
static int sample_terrain_height(fnl_state* continental, fnl_state* mountain,
                                  fnl_state* warp_noise, fnl_state* detail_noise,
                                  HMM_Vec3 unit_pos, int sea_level) {
    float n = sample_terrain_noise(continental, mountain, warp_noise, detail_noise, unit_pos);
    int height = 8 + (int)((n + 1.0f) * 0.5f * 40.0f);
    if (height < 1) height = 1;
    if (height >= 64) height = 63;
    (void)sea_level;
    return height;
}

// Perturb a base color with noise for visual variety
static HMM_Vec3 perturb_color(HMM_Vec3 base, fnl_state* color_noise, HMM_Vec3 unit_pos) {
    float cn = fnlGetNoise3D(color_noise,
        unit_pos.X * 3.0f, unit_pos.Y * 3.0f, unit_pos.Z * 3.0f);
    float variation = cn * 0.12f;  // +/- 12% brightness variation
    return (HMM_Vec3){{
        fminf(1.0f, fmaxf(0.0f, base.X + variation)),
        fminf(1.0f, fmaxf(0.0f, base.Y + variation * 0.8f)),
        fminf(1.0f, fmaxf(0.0f, base.Z + variation * 0.6f))
    }};
}

// Get terrain color based on height in meters above planet_radius.
// Colors matched to hex terrain texture atlas average colors for seamless LOD transitions.
// Uses smooth blending across biome boundaries (±100m blend zones).
static HMM_Vec3 terrain_color_m(float height_m) {
    float rel = height_m - TERRAIN_SEA_LEVEL_M;

    // Biome reference colors (matched to texture atlas averages)
    const HMM_Vec3 water_deep    = {{0.06f, 0.10f, 0.25f}};
    const HMM_Vec3 water_shallow = {{0.12f, 0.21f, 0.39f}};
    const HMM_Vec3 sand          = {{0.94f, 0.84f, 0.77f}};
    const HMM_Vec3 grass         = {{0.07f, 0.58f, 0.35f}};
    const HMM_Vec3 rock          = {{0.46f, 0.45f, 0.45f}};
    const HMM_Vec3 high_stone    = {{0.50f, 0.50f, 0.50f}};
    const HMM_Vec3 ice           = {{0.44f, 0.77f, 0.97f}};

    // Ocean: smooth deep→shallow gradient
    if (height_m < TERRAIN_SEA_LEVEL_M) {
        float depth_t = fminf(1.0f, -rel / 2000.0f);
        return vec3_lerp(water_shallow, water_deep, depth_t);
    }

    // Land biome transitions — smooth blend across ±100m zones
    float blend = 100.0f;
    float t1 = smoothstepf(200.0f  - blend, 200.0f  + blend, rel);  // Sand → Grass
    float t2 = smoothstepf(1500.0f - blend, 1500.0f + blend, rel);  // Grass → Rock
    float t3 = smoothstepf(3000.0f - blend, 3000.0f + blend, rel);  // Rock → High Stone
    float t4 = smoothstepf(4500.0f - blend, 4500.0f + blend, rel);  // High Stone → Ice

    HMM_Vec3 color = sand;
    color = vec3_lerp(color, grass,      t1);
    color = vec3_lerp(color, rock,       t2);
    color = vec3_lerp(color, high_stone, t3);
    color = vec3_lerp(color, ice,        t4);
    return color;
}

// Legacy: get terrain color from layer index (for voxel mesh)
static HMM_Vec3 terrain_color(int height, int sea_level) {
    int rel = height - sea_level;
    if (height < sea_level) {
        return (HMM_Vec3){{0.12f, 0.21f, 0.39f}};  // Water
    } else if (rel <= 1) {
        return (HMM_Vec3){{0.94f, 0.84f, 0.77f}};  // Sand
    } else if (rel <= 6) {
        return (HMM_Vec3){{0.07f, 0.58f, 0.35f}};  // Grass
    } else if (rel <= 12) {
        return (HMM_Vec3){{0.46f, 0.45f, 0.45f}};  // Stone
    } else {
        return (HMM_Vec3){{0.44f, 0.77f, 0.97f}};  // Ice
    }
}

// ---- Voxel type helpers (matches planet.c generate_terrain logic) ----

// Compute voxel type at a given layer for a column with given terrain height
static VoxelType compute_voxel_type(int layer, int terrain_height, int sea_level) {
    if (layer > terrain_height) {
        if (layer <= sea_level) return VOXEL_WATER;
        return VOXEL_AIR;
    }
    int rel = layer - sea_level;
    if (rel < -3) return VOXEL_STONE;
    if (rel < 0) return VOXEL_DIRT;
    if (rel == 0) {
        if (terrain_height <= sea_level + 1) return VOXEL_SAND;
        return VOXEL_DIRT;
    }
    if (rel == 1 && terrain_height <= sea_level + 2) return VOXEL_SAND;
    if (rel <= 6) return VOXEL_GRASS;
    if (rel <= 12) return VOXEL_STONE;
    return VOXEL_ICE;
}

// Get color for a voxel type (thread-safe, matched to hex terrain texture atlas)
static HMM_Vec3 voxel_type_color(VoxelType type) {
    switch (type) {
        case VOXEL_WATER: return (HMM_Vec3){{0.12f, 0.21f, 0.39f}};
        case VOXEL_SAND:  return (HMM_Vec3){{0.94f, 0.84f, 0.77f}};
        case VOXEL_DIRT:  return (HMM_Vec3){{0.59f, 0.29f, 0.24f}};
        case VOXEL_GRASS: return (HMM_Vec3){{0.07f, 0.58f, 0.35f}};
        case VOXEL_STONE: return (HMM_Vec3){{0.46f, 0.45f, 0.45f}};
        case VOXEL_ICE:   return (HMM_Vec3){{0.44f, 0.77f, 0.97f}};
        default:          return (HMM_Vec3){{1.0f, 0.0f, 1.0f}};  // Magenta = error
    }
}

// ---- Spherical triangle utilities ----

static void init_spherical_triangle(SphericalTriangle* tri,
                                     HMM_Vec3 v0, HMM_Vec3 v1, HMM_Vec3 v2) {
    tri->v0 = vec3_normalize(v0);
    tri->v1 = vec3_normalize(v1);
    tri->v2 = vec3_normalize(v2);

    // Centroid on unit sphere
    HMM_Vec3 c = vec3_add(vec3_add(tri->v0, tri->v1), tri->v2);
    tri->center = vec3_normalize(c);

    // Angular radius: max angle from center to any vertex
    float d0 = acosf(fminf(1.0f, fmaxf(-1.0f, vec3_dot(tri->center, tri->v0))));
    float d1 = acosf(fminf(1.0f, fmaxf(-1.0f, vec3_dot(tri->center, tri->v1))));
    float d2 = acosf(fminf(1.0f, fmaxf(-1.0f, vec3_dot(tri->center, tri->v2))));
    tri->angular_radius = fmaxf(d0, fmaxf(d1, d2));
}

// Aperture-4 subdivision: split a spherical triangle into 4 children
// by taking midpoints of edges and projecting to the unit sphere
static void subdivide_triangle(const SphericalTriangle* parent,
                                SphericalTriangle children[4]) {
    HMM_Vec3 m01 = vec3_normalize(vec3_scale(
        vec3_add(parent->v0, parent->v1), 0.5f));
    HMM_Vec3 m12 = vec3_normalize(vec3_scale(
        vec3_add(parent->v1, parent->v2), 0.5f));
    HMM_Vec3 m20 = vec3_normalize(vec3_scale(
        vec3_add(parent->v2, parent->v0), 0.5f));

    // Child 0: v0, m01, m20 (corner at v0)
    init_spherical_triangle(&children[0], parent->v0, m01, m20);
    // Child 1: m01, v1, m12 (corner at v1)
    init_spherical_triangle(&children[1], m01, parent->v1, m12);
    // Child 2: m20, m12, v2 (corner at v2)
    init_spherical_triangle(&children[2], m20, m12, parent->v2);
    // Child 3: m01, m12, m20 (center triangle)
    init_spherical_triangle(&children[3], m01, m12, m20);
}

// Check if a point on the unit sphere is inside a spherical triangle
// Uses the winding-number approach (all cross products same sign)
static bool point_in_spherical_triangle(const SphericalTriangle* tri, HMM_Vec3 p) {
    HMM_Vec3 c0 = vec3_cross(vec3_sub(tri->v1, tri->v0), vec3_sub(p, tri->v0));
    HMM_Vec3 c1 = vec3_cross(vec3_sub(tri->v2, tri->v1), vec3_sub(p, tri->v1));
    HMM_Vec3 c2 = vec3_cross(vec3_sub(tri->v0, tri->v2), vec3_sub(p, tri->v2));

    float d0 = vec3_dot(c0, p);
    float d1 = vec3_dot(c1, p);
    float d2 = vec3_dot(c2, p);

    return (d0 >= 0 && d1 >= 0 && d2 >= 0) ||
           (d0 <= 0 && d1 <= 0 && d2 <= 0);
}

// ---- Node allocation ----

static int allocate_node(LodTree* tree) {
    // First, try to reuse a dead node (freed by merge_node).
    // Dead nodes have: state=LOD_UNLOADED, parent=-1, no children, is_leaf=true.
    for (int i = 0; i < tree->node_count; i++) {
        LodNode* n = &tree->nodes[i];
        if (n->state == LOD_UNLOADED && n->parent == -1 && n->is_leaf &&
            n->children[0] == -1 && n->children[1] == -1 &&
            n->children[2] == -1 && n->children[3] == -1 &&
            n->gpu_buffer.id == SG_INVALID_ID && n->pending_job == NULL) {
            // Skip root nodes (they're always alive)
            bool is_root = false;
            for (int r = 0; r < LOD_ROOT_COUNT; r++) {
                if (tree->root_nodes[r] == i) { is_root = true; break; }
            }
            if (is_root) continue;

            // Reuse this slot
            memset(n, 0, sizeof(LodNode));
            n->parent = -1;
            for (int j = 0; j < LOD_CHILDREN; j++) n->children[j] = -1;
            n->gpu_buffer.id = SG_INVALID_ID;
            n->is_leaf = true;
            n->state = LOD_UNLOADED;
            n->pending_job = NULL;
            return i;
        }
    }

    // No dead node available — grow the array
    if (tree->node_count >= tree->node_capacity) {
        int new_cap = tree->node_capacity * 2;
        if (new_cap > LOD_MAX_NODES) new_cap = LOD_MAX_NODES;
        if (tree->node_count >= new_cap) {
            return -1;  // At capacity
        }
        tree->nodes = (LodNode*)realloc(tree->nodes, new_cap * sizeof(LodNode));
        tree->node_capacity = new_cap;
    }
    int idx = tree->node_count++;
    memset(&tree->nodes[idx], 0, sizeof(LodNode));
    tree->nodes[idx].parent = -1;
    for (int i = 0; i < LOD_CHILDREN; i++) {
        tree->nodes[idx].children[i] = -1;
    }
    tree->nodes[idx].gpu_buffer.id = SG_INVALID_ID;
    tree->nodes[idx].is_leaf = true;
    tree->nodes[idx].state = LOD_UNLOADED;
    tree->nodes[idx].pending_job = NULL;
    return idx;
}

// ---- Tree initialization ----

void lod_tree_init(LodTree* tree, float planet_radius, float layer_thickness,
                   int sea_level, int seed) {
    tree->planet_radius = planet_radius;
    tree->layer_thickness = layer_thickness;
    tree->sea_level = sea_level;
    tree->seed = seed;
    tree->camera_pos = (HMM_Vec3){{0, 0, 0}};
    tree->world_origin[0] = 0.0;
    tree->world_origin[1] = 0.0;
    tree->world_origin[2] = 0.0;
    tree->active_leaf_count = 0;
    tree->total_vertex_count = 0;
    tree->stats_frame_counter = 0;
    tree->suppress_range = 0.0f;
    tree->split_factor = LOD_SPLIT_FACTOR;

    // Initial capacity
    tree->node_capacity = 256;
    tree->nodes = (LodNode*)calloc(tree->node_capacity, sizeof(LodNode));
    tree->node_count = 0;

    // Create job system for threaded mesh generation
    tree->jobs = job_system_create(LOD_NUM_WORKERS);

    // Create 20 root nodes from icosahedron faces
    for (int i = 0; i < ICO_FACE_COUNT; i++) {
        int idx = allocate_node(tree);
        if (idx < 0) break;
        tree->root_nodes[i] = idx;

        LodNode* node = &tree->nodes[idx];
        node->depth = 0;
        node->parent = -1;

        HMM_Vec3 v0 = vec3_normalize(ICO_VERTICES[ICO_FACES[i][0]]);
        HMM_Vec3 v1 = vec3_normalize(ICO_VERTICES[ICO_FACES[i][1]]);
        HMM_Vec3 v2 = vec3_normalize(ICO_VERTICES[ICO_FACES[i][2]]);

        init_spherical_triangle(&node->tri, v0, v1, v2);
    }

    // Precompute depth-based arc lengths for consistent split/merge thresholds.
    // Use the average root angular_radius as baseline, halving each depth.
    // This ensures all patches at the same depth have the same split criterion,
    // regardless of aperture-4 subdivision asymmetry (corner vs center children).
    {
        float avg_root_ar = 0.0f;
        for (int i = 0; i < LOD_ROOT_COUNT; i++) {
            avg_root_ar += tree->nodes[tree->root_nodes[i]].tri.angular_radius;
        }
        avg_root_ar /= (float)LOD_ROOT_COUNT;

        for (int d = 0; d <= LOD_MAX_DEPTH; d++) {
            float ar = avg_root_ar;
            for (int i = 0; i < d; i++) ar *= 0.5f;
            tree->depth_arc[d] = ar * planet_radius;
        }
    }

    printf("[LOD] Tree initialized: radius=%.0f, %d root nodes, max_depth=%d, %d workers\n",
        planet_radius, LOD_ROOT_COUNT, LOD_MAX_DEPTH, LOD_NUM_WORKERS);
    fflush(stdout);
}

void lod_tree_destroy(LodTree* tree) {
    // Destroy job system first (waits for pending jobs to finish)
    if (tree->jobs) {
        job_system_destroy(tree->jobs);
        tree->jobs = NULL;
    }

    // Clean up all nodes (including any completed pending jobs)
    for (int i = 0; i < tree->node_count; i++) {
        LodNode* node = &tree->nodes[i];

        // Free pending job data if any
        if (node->pending_job) {
            MeshGenJob* job = (MeshGenJob*)node->pending_job;
            if (job->vertices) {
                free(job->vertices);
            }
            free(job);
            node->pending_job = NULL;
        }

        free_node_mesh(node);
        if (node->gpu_buffer.id != SG_INVALID_ID) {
            sg_destroy_buffer(node->gpu_buffer);
        }
    }
    free(tree->nodes);
    tree->nodes = NULL;
    tree->node_count = 0;
    tree->node_capacity = 0;
}

// ---- Hex grid constants (for LOD depth-13 hex mesh generation) ----
#define LOD_HEX_RADIUS      1.0f
#define LOD_HEX_COL_SPACING (1.5f * LOD_HEX_RADIUS)
#define LOD_HEX_ROW_SPACING (1.7320508f * LOD_HEX_RADIUS)  // sqrt(3)

static const float LOD_HEX_COS[6] = {
    1.0f, 0.5f, -0.5f, -1.0f, -0.5f, 0.5f
};
static const float LOD_HEX_SIN[6] = {
    0.0f, 0.8660254f, 0.8660254f, 0.0f, -0.8660254f, -0.8660254f
};

// ---- Mesh generation ----

// Helper: emit a triangle into the vertex array
// Automatically ensures winding matches the given normal direction (CCW when
// viewed from the normal side), so that SG_FACEWINDING_CCW renders front faces.
static void lod_emit_tri(LodVertex** verts, int* count, int* cap,
                          HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                          HMM_Vec3 normal, HMM_Vec3 color) {
    // Ensure winding matches normal (CCW from the normal's side = front face)
    HMM_Vec3 e1 = vec3_sub(p1, p0);
    HMM_Vec3 e2 = vec3_sub(p2, p0);
    HMM_Vec3 face_n = vec3_cross(e1, e2);
    if (vec3_dot(face_n, normal) < 0.0f) {
        HMM_Vec3 tmp = p1; p1 = p2; p2 = tmp;
    }

    if (*count + 3 > *cap) {
        *cap *= 2;
        *verts = (LodVertex*)realloc(*verts, *cap * sizeof(LodVertex));
    }
    LodVertex* v = &(*verts)[*count];

    v[0].pos[0] = p0.X; v[0].pos[1] = p0.Y; v[0].pos[2] = p0.Z;
    v[0].normal[0] = normal.X; v[0].normal[1] = normal.Y; v[0].normal[2] = normal.Z;
    v[0].color[0] = color.X; v[0].color[1] = color.Y; v[0].color[2] = color.Z;

    v[1].pos[0] = p1.X; v[1].pos[1] = p1.Y; v[1].pos[2] = p1.Z;
    v[1].normal[0] = normal.X; v[1].normal[1] = normal.Y; v[1].normal[2] = normal.Z;
    v[1].color[0] = color.X; v[1].color[1] = color.Y; v[1].color[2] = color.Z;

    v[2].pos[0] = p2.X; v[2].pos[1] = p2.Y; v[2].pos[2] = p2.Z;
    v[2].normal[0] = normal.X; v[2].normal[1] = normal.Y; v[2].normal[2] = normal.Z;
    v[2].color[0] = color.X; v[2].color[1] = color.Y; v[2].color[2] = color.Z;

    *count += 3;
}

// ---- Hex grid mesh generation for depth-13 LOD nodes ----

// Texture atlas mapping (must match hex_terrain.c)
#define LOD_HEX_ATLAS_TILES 9
#define LOD_ATLAS_WATER  0
#define LOD_ATLAS_SAND   1
#define LOD_ATLAS_GRASS  3
#define LOD_ATLAS_STONE  4
#define LOD_ATLAS_ICE    5

static int lod_hex_atlas(float height_m) {
    float rel = height_m - TERRAIN_SEA_LEVEL_M;
    if (height_m < TERRAIN_SEA_LEVEL_M) return LOD_ATLAS_WATER;
    if (rel < 200.0f) return LOD_ATLAS_SAND;
    if (rel < 1500.0f) return LOD_ATLAS_GRASS;
    if (rel < 3000.0f) return LOD_ATLAS_STONE;
    if (rel < 4500.0f) return LOD_ATLAS_STONE;
    return LOD_ATLAS_ICE;
}

typedef struct { float u, v; } LodHexUV;

static LodHexUV lod_hex_top_uv(float vx, float vz, float cx, float cz, int atlas_idx) {
    float local_u = (vx - cx) / (2.0f * LOD_HEX_RADIUS) + 0.5f;
    float local_v = (vz - cz) / (LOD_HEX_ROW_SPACING) + 0.5f;
    LodHexUV uv;
    uv.u = ((float)atlas_idx + local_u) / (float)LOD_HEX_ATLAS_TILES;
    uv.v = local_v;
    return uv;
}

// ---- Hex neighbor lookup (axial coordinate system) ----

static const int LOD_AX_NEIGHBORS[6][2] = {
    { +1,  0 },   // dir 0: NE
    { +1, -1 },   // dir 1: SE
    {  0, -1 },   // dir 2: S
    { -1,  0 },   // dir 3: SW
    { -1, +1 },   // dir 4: NW
    {  0, +1 },   // dir 5: N
};

// Map neighbor direction to hex edge index (edge i = vertex[i] to vertex[(i+1)%6])
static const int LOD_DIR_TO_EDGE[6] = {
    0,  // dir 0 (NE) → edge 0 (v0-v1)
    5,  // dir 1 (SE) → edge 5 (v5-v0)
    4,  // dir 2 (S)  → edge 4 (v4-v5)
    3,  // dir 3 (SW) → edge 3 (v3-v4)
    2,  // dir 4 (NW) → edge 2 (v2-v3)
    1,  // dir 5 (N)  → edge 1 (v1-v2)
};

static void lod_hex_neighbor(int col, int row, int dir, int* ncol, int* nrow) {
    int q = col;
    int r = row - (col - (col & 1)) / 2;
    int nq = q + LOD_AX_NEIGHBORS[dir][0];
    int nr = r + LOD_AX_NEIGHBORS[dir][1];
    *ncol = nq;
    *nrow = nr + (nq - (nq & 1)) / 2;
}

// ---- Wall atlas helpers ----

#define LOD_ATLAS_DIRT       2
#define LOD_ATLAS_DIRT_GRASS 7
#define LOD_ATLAS_DIRT_SNOW  8

#define LOD_HEX_SURFACE_BIAS 0.05f
#define LOD_HEX_MAX_WALL     32

static int lod_side_atlas(int cap_atlas) {
    if (cap_atlas == LOD_ATLAS_GRASS) return LOD_ATLAS_DIRT_GRASS;
    if (cap_atlas == LOD_ATLAS_ICE)   return LOD_ATLAS_DIRT_SNOW;
    return cap_atlas;
}

static int lod_wall_atlas(int cap_atlas, float surface_h, float wall_y) {
    float depth = surface_h - wall_y;
    if (depth < 1.0f) return lod_side_atlas(cap_atlas);
    if (depth < 3.0f) return LOD_ATLAS_DIRT;
    return LOD_ATLAS_STONE;
}

// Emit a hex triangle into HexVertex array.
// winding_normal: used for CCW winding check (e.g., wall_outward for walls).
// shading_normal: stored in vertex for lighting (e.g., local_up for voxel-style).
static void lod_hex_emit_tri(HexVertex** verts, int* count, int* cap,
                              HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                              HMM_Vec3 winding_normal, HMM_Vec3 shading_normal,
                              LodHexUV uv0, LodHexUV uv1, LodHexUV uv2,
                              HMM_Vec3 color) {
    HMM_Vec3 e1 = vec3_sub(p1, p0);
    HMM_Vec3 e2 = vec3_sub(p2, p0);
    HMM_Vec3 face_n = vec3_cross(e1, e2);
    if (vec3_dot(face_n, winding_normal) < 0.0f) {
        HMM_Vec3 tmp = p1; p1 = p2; p2 = tmp;
        LodHexUV utmp = uv1; uv1 = uv2; uv2 = utmp;
    }

    if (*count + 3 > *cap) {
        *cap = (*cap) * 2;
        *verts = (HexVertex*)realloc(*verts, *cap * sizeof(HexVertex));
    }
    HexVertex* v = &(*verts)[*count];

    v[0].pos[0] = p0.X; v[0].pos[1] = p0.Y; v[0].pos[2] = p0.Z;
    v[0].normal[0] = shading_normal.X; v[0].normal[1] = shading_normal.Y; v[0].normal[2] = shading_normal.Z;
    v[0].uv[0] = uv0.u; v[0].uv[1] = uv0.v;
    v[0].color[0] = color.X; v[0].color[1] = color.Y; v[0].color[2] = color.Z;

    v[1].pos[0] = p1.X; v[1].pos[1] = p1.Y; v[1].pos[2] = p1.Z;
    v[1].normal[0] = shading_normal.X; v[1].normal[1] = shading_normal.Y; v[1].normal[2] = shading_normal.Z;
    v[1].uv[0] = uv1.u; v[1].uv[1] = uv1.v;
    v[1].color[0] = color.X; v[1].color[1] = color.Y; v[1].color[2] = color.Z;

    v[2].pos[0] = p2.X; v[2].pos[1] = p2.Y; v[2].pos[2] = p2.Z;
    v[2].normal[0] = shading_normal.X; v[2].normal[1] = shading_normal.Y; v[2].normal[2] = shading_normal.Z;
    v[2].uv[0] = uv2.u; v[2].uv[1] = uv2.v;
    v[2].color[0] = color.X; v[2].color[1] = color.Y; v[2].color[2] = color.Z;

    *count += 3;
}

// Point-in-triangle test (2D) with margin.
static bool point_in_triangle_margin(float px, float pz,
                                      const float tx[3], const float tz[3],
                                      float margin) {
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        float ex = tx[j] - tx[i];
        float ez = tz[j] - tz[i];
        float nx = -ez;
        float nz = ex;
        float len = sqrtf(nx * nx + nz * nz);
        if (len < 1e-8f) continue;
        nx /= len;
        nz /= len;
        float d = (px - tx[i]) * nx + (pz - tz[i]) * nz;
        if (d < -margin) return false;
    }
    return true;
}

// Generate hex prism mesh filling a spherical triangle.
// Each hex gets a quantized height (1m steps), flat cap perpendicular to planet
// radial direction, and vertical walls where adjacent hexes differ in height.
// Uses a shared tangent frame (hex_origin/east/north) so all depth-13 patches
// generate hexes from the same global grid — no seams between patches.
// Outputs HexVertex[] (pos + normal + uv) for rendering with hex terrain pipeline.
static void generate_hex_mesh_for_triangle(
    const SphericalTriangle* tri, float planet_radius, int seed,
    HMM_Vec3 hex_origin, HMM_Vec3 hex_east, HMM_Vec3 hex_north,
    HexVertex** out_vertices, int* out_count)
{
    // 1. Use shared tangent frame (same for ALL depth-13 patches → seamless tiling)
    HMM_Vec3 center = hex_origin;  // frame center, NOT per-triangle center
    HMM_Vec3 east = hex_east;
    HMM_Vec3 north = hex_north;

    // 2. Project triangle vertices to 2D tangent plane
    float tvx[3], tvz[3];
    HMM_Vec3 verts[3] = { tri->v0, tri->v1, tri->v2 };
    for (int i = 0; i < 3; i++) {
        float dot_vc = vec3_dot(verts[i], center);
        HMM_Vec3 offset = vec3_sub(verts[i], vec3_scale(center, dot_vc));
        tvx[i] = vec3_dot(offset, east) * planet_radius;
        tvz[i] = vec3_dot(offset, north) * planet_radius;
    }

    // 3. Compute bounding box with margin
    float margin = LOD_HEX_RADIUS * 1.5f;
    float min_x = fminf(tvx[0], fminf(tvx[1], tvx[2])) - margin;
    float max_x = fmaxf(tvx[0], fmaxf(tvx[1], tvx[2])) + margin;
    float min_z = fminf(tvz[0], fminf(tvz[1], tvz[2])) - margin;
    float max_z = fmaxf(tvz[0], fmaxf(tvz[1], tvz[2])) + margin;

    // 4. Create noise states
    fnl_state continental = create_continental_noise(seed);
    fnl_state mountain_n = create_mountain_noise(seed);
    fnl_state warp = create_warp_noise(seed);
    fnl_state detail = create_detail_noise(seed);
    fnl_state cnoise = create_color_noise(seed);

    HMM_Vec3 center_scaled = vec3_scale(center, planet_radius);

    // 5. Column/row range from bounding box.
    // All patches share the same tangent frame, so col=0/row=0 maps to the
    // frame center for every patch.  Simple direct grid enumeration works.
    int col_min = (int)floorf(min_x / LOD_HEX_COL_SPACING);
    int col_max = (int)ceilf(max_x / LOD_HEX_COL_SPACING);
    int row_min = (int)floorf(min_z / LOD_HEX_ROW_SPACING) - 1;
    int row_max = (int)ceilf(max_z / LOD_HEX_ROW_SPACING) + 1;

    // Helper: convert (col, row) to tangent-plane position (hx, hz).
    #define HEX_POS(col, row, hx, hz) do { \
        (hx) = (float)(col) * LOD_HEX_COL_SPACING; \
        (hz) = (((col) & 1) ? ((row) + 0.5f) : (float)(row)) * LOD_HEX_ROW_SPACING; \
    } while(0)

    // 6. Phase 1: Build height map — sample terrain at hex centers, quantize
    int grid_cols = col_max - col_min + 1;
    int grid_rows = row_max - row_min + 1;
    int grid_size = grid_cols * grid_rows;
    int16_t* hm_heights = (int16_t*)malloc(grid_size * sizeof(int16_t));
    uint8_t* hm_atlas = (uint8_t*)malloc(grid_size * sizeof(uint8_t));
    HMM_Vec3* hm_colors = (HMM_Vec3*)malloc(grid_size * sizeof(HMM_Vec3));
    for (int i = 0; i < grid_size; i++) {
        hm_heights[i] = INT16_MIN;
    }

    for (int col = col_min; col <= col_max; col++) {
        for (int row = row_min; row <= row_max; row++) {
            float hx, hz;
            HEX_POS(col, row, hx, hz);
            if (!point_in_triangle_margin(hx, hz, tvx, tvz, margin))
                continue;

            HMM_Vec3 cdir = vec3_normalize(
                vec3_add(center_scaled,
                    vec3_add(vec3_scale(east, hx), vec3_scale(north, hz))));
            float h_m = sample_terrain_height_m(&continental, &mountain_n, &warp, &detail, cdir);
            float eff = fmaxf(h_m, TERRAIN_SEA_LEVEL_M);
            int h = (int)ceilf(eff);
            if (h < 0) h = 0;

            int gi = (col - col_min) * grid_rows + (row - row_min);

            hm_heights[gi] = (int16_t)h;
            hm_atlas[gi] = (uint8_t)lod_hex_atlas(h_m);
            hm_colors[gi] = perturb_color(terrain_color_m(h_m), &cnoise, cdir);
        }
    }

    // 7. Phase 1.5: Compute terrain slope normals from height map.
    //    The slope normal encodes the macro terrain orientation so that hexes
    //    on the shadow side of a mountain are shaded darker (mountain shadow illusion).
    HMM_Vec3* hm_slopes = (HMM_Vec3*)malloc(grid_size * sizeof(HMM_Vec3));
    for (int col = col_min; col <= col_max; col++) {
        for (int row = row_min; row <= row_max; row++) {
            int gi = (col - col_min) * grid_rows + (row - row_min);
            if (hm_heights[gi] == INT16_MIN) {
                hm_slopes[gi] = center;  // fallback
                continue;
            }

            float hx, hz;
            HEX_POS(col, row, hx, hz);
            int h = hm_heights[gi];

            // Accumulate gradient from neighbor height differences
            float grad_e = 0.0f, grad_n = 0.0f;
            int ncount = 0;
            for (int dir = 0; dir < 6; dir++) {
                int ncol, nrow;
                lod_hex_neighbor(col, row, dir, &ncol, &nrow);
                if (ncol < col_min || ncol > col_max ||
                    nrow < row_min || nrow > row_max) continue;
                int ni = (ncol - col_min) * grid_rows + (nrow - row_min);
                if (hm_heights[ni] == INT16_MIN) continue;

                float nhx, nhz;
                HEX_POS(ncol, nrow, nhx, nhz);
                float dx = nhx - hx;
                float dz = nhz - hz;
                float dh = (float)(hm_heights[ni] - h);
                float dist_sq = dx * dx + dz * dz;
                if (dist_sq < 0.01f) continue;
                grad_e += dh * dx / dist_sq;
                grad_n += dh * dz / dist_sq;
                ncount++;
            }
            if (ncount > 0) {
                grad_e /= ncount;
                grad_n /= ncount;
            }

            // Slope normal = local_up tilted by the terrain gradient
            HMM_Vec3 cdir = vec3_normalize(
                vec3_add(center_scaled,
                    vec3_add(vec3_scale(east, hx), vec3_scale(north, hz))));
            HMM_Vec3 slope = vec3_normalize(
                vec3_sub(cdir,
                    vec3_add(vec3_scale(east, grad_e), vec3_scale(north, grad_n))));
            // Ensure slope still points outward
            if (vec3_dot(slope, cdir) < 0.1f) slope = cdir;
            hm_slopes[gi] = slope;
        }
    }

    // 7.5. Apply slope-based brightness to vertex colors.
    //       Steeper slopes get darker, matching how the planet shader's AO + NdotL
    //       look on smooth LOD terrain. This makes the distance fade seamless.
    for (int col = col_min; col <= col_max; col++) {
        for (int row = row_min; row <= row_max; row++) {
            int gi = (col - col_min) * grid_rows + (row - row_min);
            if (hm_heights[gi] == INT16_MIN) continue;

            float hx, hz;
            HEX_POS(col, row, hx, hz);
            HMM_Vec3 cdir = vec3_normalize(
                vec3_add(center_scaled,
                    vec3_add(vec3_scale(east, hx), vec3_scale(north, hz))));

            // Slope vs vertical: steeper = darker (mimics AO + diffuse falloff)
            float slope_dot = vec3_dot(hm_slopes[gi], cdir);
            float brightness = 0.45f + 0.55f * fmaxf(0.0f, slope_dot);
            hm_colors[gi] = vec3_scale(hm_colors[gi], brightness);
        }
    }

    // 8. Allocate output (cap=12 verts + walls ~36 avg per hex)
    int vert_cap = grid_size * 48;
    if (vert_cap < 64) vert_cap = 64;
    *out_vertices = (HexVertex*)malloc(vert_cap * sizeof(HexVertex));
    *out_count = 0;

    // 9. Phase 2: Generate cap + wall geometry
    for (int col = col_min; col <= col_max; col++) {
        for (int row = row_min; row <= row_max; row++) {
            int gi = (col - col_min) * grid_rows + (row - row_min);
            if (hm_heights[gi] == INT16_MIN) continue;

            int h = hm_heights[gi];
            int cap_atlas_idx = hm_atlas[gi];
            HMM_Vec3 hex_color = hm_colors[gi];
            HMM_Vec3 slope_normal = hm_slopes[gi];

            float hx, hz;
            HEX_POS(col, row, hx, hz);

            // Hex center direction = planet radial up
            HMM_Vec3 cdir = vec3_normalize(
                vec3_add(center_scaled,
                    vec3_add(vec3_scale(east, hx), vec3_scale(north, hz))));
            HMM_Vec3 local_up = cdir;

            // Cap radius: planet_radius + quantized_height + surface bias
            float cap_r = planet_radius + (float)h + LOD_HEX_SURFACE_BIAS;

            // Compute 6 corner positions — all at the same quantized height
            HMM_Vec3 hex_verts[6];
            HMM_Vec3 hex_dirs[6];
            LodHexUV hex_uvs[6];
            for (int i = 0; i < 6; i++) {
                float vx = hx + LOD_HEX_RADIUS * LOD_HEX_COS[i];
                float vz = hz + LOD_HEX_RADIUS * LOD_HEX_SIN[i];
                hex_dirs[i] = vec3_normalize(
                    vec3_add(center_scaled,
                        vec3_add(vec3_scale(east, vx), vec3_scale(north, vz))));
                hex_verts[i] = vec3_scale(hex_dirs[i], cap_r);
                hex_uvs[i] = lod_hex_top_uv(vx, vz, hx, hz, cap_atlas_idx);
            }

            // Emit 4 cap triangles (fan from vertex 0)
            // Winding = local_up, shading = slope_normal (terrain shadow illusion)
            for (int i = 0; i < 4; i++) {
                lod_hex_emit_tri(out_vertices, out_count, &vert_cap,
                    hex_verts[0], hex_verts[i + 1], hex_verts[i + 2],
                    local_up, slope_normal,
                    hex_uvs[0], hex_uvs[i + 1], hex_uvs[i + 2],
                    hex_color);
            }

            // ---- Walls: check 6 neighbors ----
            for (int dir = 0; dir < 6; dir++) {
                int ncol, nrow;
                lod_hex_neighbor(col, row, dir, &ncol, &nrow);

                // Look up neighbor height in grid
                if (ncol < col_min || ncol > col_max ||
                    nrow < row_min || nrow > row_max)
                    continue;  // outside grid — skip (shared frame means no seams)
                int ni = (ncol - col_min) * grid_rows + (nrow - row_min);
                if (hm_heights[ni] == INT16_MIN) continue;  // invalid neighbor
                int nh = hm_heights[ni];
                if (h <= nh) continue;  // neighbor same or taller — no wall

                // Edge vertices for this wall direction
                int edge = LOD_DIR_TO_EDGE[dir];
                int vi0 = edge;
                int vi1 = (edge + 1) % 6;

                // Wall outward direction (hex center → edge midpoint)
                float mx = (LOD_HEX_COS[vi0] + LOD_HEX_COS[vi1]) * 0.5f * LOD_HEX_RADIUS;
                float mz = (LOD_HEX_SIN[vi0] + LOD_HEX_SIN[vi1]) * 0.5f * LOD_HEX_RADIUS;
                HMM_Vec3 wall_outward = vec3_normalize(
                    vec3_add(vec3_scale(east, mx), vec3_scale(north, mz)));

                // Cap wall height to limit geometry
                int wall_bottom = nh;
                int wall_top = h;
                if (wall_top - wall_bottom > LOD_HEX_MAX_WALL)
                    wall_bottom = wall_top - LOD_HEX_MAX_WALL;

                // Emit per-layer wall quads
                // Winding = wall_outward, shading = slope_normal (shadow continuity)
                for (int layer = wall_bottom; layer < wall_top; layer++) {
                    float bot_r = planet_radius + (float)layer + LOD_HEX_SURFACE_BIAS;
                    float top_r = planet_radius + (float)(layer + 1) + LOD_HEX_SURFACE_BIAS;

                    HMM_Vec3 p0_bot = vec3_scale(hex_dirs[vi0], bot_r);
                    HMM_Vec3 p1_bot = vec3_scale(hex_dirs[vi1], bot_r);
                    HMM_Vec3 p0_top = vec3_scale(hex_dirs[vi0], top_r);
                    HMM_Vec3 p1_top = vec3_scale(hex_dirs[vi1], top_r);

                    // Wall texture by depth below surface
                    float wall_y = (float)layer + 0.5f;
                    float surface_h = (float)h;
                    int wall_atlas_idx = lod_wall_atlas(cap_atlas_idx, surface_h, wall_y);

                    LodHexUV wuv_bl = { (float)wall_atlas_idx / (float)LOD_HEX_ATLAS_TILES, 1.0f };
                    LodHexUV wuv_br = { ((float)wall_atlas_idx + 1.0f) / (float)LOD_HEX_ATLAS_TILES, 1.0f };
                    LodHexUV wuv_tl = { (float)wall_atlas_idx / (float)LOD_HEX_ATLAS_TILES, 0.0f };
                    LodHexUV wuv_tr = { ((float)wall_atlas_idx + 1.0f) / (float)LOD_HEX_ATLAS_TILES, 0.0f };

                    lod_hex_emit_tri(out_vertices, out_count, &vert_cap,
                        p0_bot, p1_bot, p1_top,
                        wall_outward, slope_normal,
                        wuv_bl, wuv_br, wuv_tr,
                        hex_color);
                    lod_hex_emit_tri(out_vertices, out_count, &vert_cap,
                        p0_bot, p1_top, p0_top,
                        wall_outward, slope_normal,
                        wuv_bl, wuv_tr, wuv_tl,
                        hex_color);
                }
            }
        }
    }

    // 10. Cleanup
    free(hm_heights);
    free(hm_atlas);
    free(hm_colors);
    free(hm_slopes);
    #undef HEX_POS
}

// Generate mesh for a coarse node (tessellated spherical triangle).
// Thread-safe: takes all inputs by value, writes to output pointers.
static void generate_coarse_mesh_params(
    const SphericalTriangle* tri, int depth,
    float planet_radius, float layer_thickness, int sea_level, int seed,
    LodVertex** out_vertices, int* out_count)
{
    // Tessellation resolution depends on depth.
    // Coarse depths (0-3) are fallback-only, kept light.
    // Mid depths (4-9) are visible at the horizon, need smoother tessellation.
    // Depths 10-12 are close-range, need high tessellation.
    int tess = 4;
    if (depth >= 4) tess = 6;
    if (depth >= 8) tess = 8;
    if (depth >= 10) tess = 16;
    if (depth >= 11) tess = 24;
    if (depth >= 12) tess = 32;

    fnl_state continental = create_continental_noise(seed);
    fnl_state mountain_n = create_mountain_noise(seed);
    fnl_state warp = create_warp_noise(seed);
    fnl_state detail = create_detail_noise(seed);
    fnl_state cnoise = create_color_noise(seed);

    int max_verts = tess * tess * 6 * 3;
    *out_vertices = (LodVertex*)malloc(max_verts * sizeof(LodVertex));
    *out_count = 0;
    int cap = max_verts;

    // Barycentric tessellation of the spherical triangle
    int num_points = (tess + 1) * (tess + 2) / 2;
    HMM_Vec3* points = (HMM_Vec3*)malloc(num_points * sizeof(HMM_Vec3));
    float* heights = (float*)malloc(num_points * sizeof(float));
    HMM_Vec3* colors = (HMM_Vec3*)malloc(num_points * sizeof(HMM_Vec3));

    int idx = 0;
    for (int row = 0; row <= tess; row++) {
        for (int col = 0; col <= tess - row; col++) {
            float u = (float)col / (float)tess;
            float v = (float)row / (float)tess;
            float w = 1.0f - u - v;

            HMM_Vec3 p = vec3_add(
                vec3_add(
                    vec3_scale(tri->v0, w),
                    vec3_scale(tri->v1, u)
                ),
                vec3_scale(tri->v2, v)
            );
            p = vec3_normalize(p);

            float h_m = sample_terrain_height_m(&continental, &mountain_n, &warp, &detail, p);
            // Clamp to sea level for ocean surface
            float effective_h_m = fmaxf(h_m, TERRAIN_SEA_LEVEL_M);

            float radius = planet_radius + effective_h_m;

            points[idx] = vec3_scale(p, radius);
            heights[idx] = h_m;
            colors[idx] = perturb_color(terrain_color_m(h_m), &cnoise, p);
            idx++;
        }
    }

    // Generate triangles from the tessellated grid
    int* row_offsets = (int*)malloc((tess + 1) * sizeof(int));
    row_offsets[0] = 0;
    for (int r = 1; r <= tess; r++) {
        row_offsets[r] = row_offsets[r - 1] + (tess - (r - 1) + 1);
    }

    for (int row = 0; row < tess; row++) {
        int cols_this_row = tess - row + 1;
        for (int col = 0; col < cols_this_row - 1; col++) {
            // Upward-pointing triangle
            int i0 = row_offsets[row] + col;
            int i1 = row_offsets[row] + col + 1;
            int i2 = row_offsets[row + 1] + col;

            HMM_Vec3 p0 = points[i0];
            HMM_Vec3 p1 = points[i1];
            HMM_Vec3 p2 = points[i2];

            HMM_Vec3 edge1 = vec3_sub(p1, p0);
            HMM_Vec3 edge2 = vec3_sub(p2, p0);
            HMM_Vec3 normal = vec3_normalize(vec3_cross(edge1, edge2));

            HMM_Vec3 face_center = vec3_scale(vec3_add(vec3_add(p0, p1), p2), 1.0f / 3.0f);
            if (vec3_dot(normal, face_center) < 0) {
                normal = vec3_scale(normal, -1.0f);
                HMM_Vec3 tmp = p1; p1 = p2; p2 = tmp;
            }

            HMM_Vec3 c = vec3_scale(
                vec3_add(vec3_add(colors[i0], colors[i1]), colors[i2]),
                1.0f / 3.0f);

            lod_emit_tri(out_vertices, out_count, &cap,
                         p0, p1, p2, normal, c);

            // Downward-pointing triangle
            if (col < cols_this_row - 2) {
                int i3 = row_offsets[row + 1] + col;
                int i4 = row_offsets[row] + col + 1;
                int i5 = row_offsets[row + 1] + col + 1;

                HMM_Vec3 q0 = points[i3];
                HMM_Vec3 q1 = points[i4];
                HMM_Vec3 q2 = points[i5];

                HMM_Vec3 e1 = vec3_sub(q1, q0);
                HMM_Vec3 e2 = vec3_sub(q2, q0);
                HMM_Vec3 n = vec3_normalize(vec3_cross(e1, e2));

                HMM_Vec3 fc = vec3_scale(vec3_add(vec3_add(q0, q1), q2), 1.0f / 3.0f);
                if (vec3_dot(n, fc) < 0) {
                    n = vec3_scale(n, -1.0f);
                    HMM_Vec3 tmp = q1; q1 = q2; q2 = tmp;
                }

                HMM_Vec3 c2 = vec3_scale(
                    vec3_add(vec3_add(colors[i3], colors[i4]), colors[i5]),
                    1.0f / 3.0f);

                lod_emit_tri(out_vertices, out_count, &cap,
                             q0, q1, q2, n, c2);
            }
        }
    }

    // ---- Boundary skirts: hide gaps between patches at different LOD levels ----
    // Only needed at depths where neighbors may have different LOD levels.
    // At depth >= LOD_MAX_DEPTH-1, all neighbors are the same depth — no seams.
    if (depth < LOD_MAX_DEPTH - 1) {
        float patch_world_size = tri->angular_radius * planet_radius;
        float skirt_drop = fmaxf(3.0f * layer_thickness, patch_world_size * 0.015f);

        #define EMIT_SKIRT(ia, ib) do { \
            HMM_Vec3 pa = points[(ia)], pb = points[(ib)]; \
            float ra = sqrtf(vec3_dot(pa, pa)); \
            float rb = sqrtf(vec3_dot(pb, pb)); \
            HMM_Vec3 ua = vec3_scale(pa, 1.0f / ra); \
            HMM_Vec3 ub = vec3_scale(pb, 1.0f / rb); \
            HMM_Vec3 pa_bot = vec3_scale(ua, ra - skirt_drop); \
            HMM_Vec3 pb_bot = vec3_scale(ub, rb - skirt_drop); \
            HMM_Vec3 smid = vec3_normalize(vec3_scale(vec3_add(ua, ub), 0.5f)); \
            HMM_Vec3 outward = vec3_normalize(vec3_sub(smid, tri->center)); \
            HMM_Vec3 skc = vec3_scale(vec3_add(colors[(ia)], colors[(ib)]), 0.5f); \
            lod_emit_tri(out_vertices, out_count, &cap, pa, pb_bot, pb, outward, skc); \
            lod_emit_tri(out_vertices, out_count, &cap, pa, pa_bot, pb_bot, outward, skc); \
        } while(0)

        for (int col = 0; col < tess; col++) {
            EMIT_SKIRT(row_offsets[0] + col, row_offsets[0] + col + 1);
        }
        for (int row = 0; row < tess; row++) {
            EMIT_SKIRT(row_offsets[row + 1], row_offsets[row]);
        }
        for (int row = 0; row < tess; row++) {
            EMIT_SKIRT(row_offsets[row] + (tess - row),
                       row_offsets[row + 1] + (tess - row - 1));
        }

        #undef EMIT_SKIRT
    }

    free(points);
    free(heights);
    free(colors);
    free(row_offsets);
}



// Apply floating origin offset to all vertex positions in a mesh.
// Called after mesh generation (normals/winding already computed from world coords).
// Vertex positions are translated so they're relative to origin, keeping float32 precise.
// Works with any vertex format where pos[3] is the first field (LodVertex, HexVertex).
static void apply_origin_offset(void* verts, int count, int stride, const float origin[3]) {
    char* base = (char*)verts;
    for (int i = 0; i < count; i++) {
        float* pos = (float*)(base + i * stride);
        pos[0] -= origin[0];
        pos[1] -= origin[1];
        pos[2] -= origin[2];
    }
}

// Worker thread callback for async mesh generation
static void mesh_gen_worker(void* data) {
    MeshGenJob* job = (MeshGenJob*)data;

    if (job->depth >= LOD_MAX_DEPTH) {
        // Hex mesh: output is HexVertex[]
        HexVertex* hex_verts = NULL;
        int hex_count = 0;
        generate_hex_mesh_for_triangle(
            &job->tri, job->planet_radius, job->seed,
            job->hex_frame_origin, job->hex_frame_east, job->hex_frame_north,
            &hex_verts, &hex_count);
        job->vertices = hex_verts;
        job->vertex_count = hex_count;
        job->vertex_stride = sizeof(HexVertex);
        job->is_hex_mesh = true;
    } else {
        // Tessellated mesh: smooth terrain (all non-hex depths)
        LodVertex* lod_verts = NULL;
        int lod_count = 0;
        generate_coarse_mesh_params(
            &job->tri, job->depth,
            job->planet_radius, job->layer_thickness, job->sea_level, job->seed,
            &lod_verts, &lod_count);
        job->vertices = lod_verts;
        job->vertex_count = lod_count;
        job->vertex_stride = sizeof(LodVertex);
        job->is_hex_mesh = false;
    }

    // Apply floating origin: shift vertex positions so they're near camera
    apply_origin_offset(job->vertices, job->vertex_count, job->vertex_stride, job->origin);

    // Signal completion (volatile write - safe on x86 with store ordering)
    job->completed = 1;
}

// Synchronous mesh generation (fallback when no job system)
static void generate_node_mesh(LodTree* tree, LodNode* node) {
    if (node->depth >= LOD_MAX_DEPTH) {
        HexVertex* hex_verts = NULL;
        int hex_count = 0;
        generate_hex_mesh_for_triangle(
            &node->tri, tree->planet_radius, tree->seed,
            tree->hex_frame_origin, tree->hex_frame_east, tree->hex_frame_north,
            &hex_verts, &hex_count);
        node->cpu_vertices = hex_verts;
        node->cpu_vertex_count = hex_count;
        node->cpu_vertex_stride = sizeof(HexVertex);
        node->is_hex_mesh = true;
    } else {
        LodVertex* lod_verts = NULL;
        int lod_count = 0;
        generate_coarse_mesh_params(
            &node->tri, node->depth,
            tree->planet_radius, tree->layer_thickness, tree->sea_level, tree->seed,
            &lod_verts, &lod_count);
        node->cpu_vertices = lod_verts;
        node->cpu_vertex_count = lod_count;
        node->cpu_vertex_stride = sizeof(LodVertex);
        node->is_hex_mesh = false;
    }
    float origin[3] = {(float)tree->world_origin[0],
                       (float)tree->world_origin[1],
                       (float)tree->world_origin[2]};
    apply_origin_offset(node->cpu_vertices, node->cpu_vertex_count,
                        node->cpu_vertex_stride, origin);
    node->state = LOD_READY;
}

// Submit async mesh generation job for a node.
// Returns true if job was submitted, false if queue is full (node stays UNLOADED).
static bool submit_mesh_job(LodTree* tree, int node_idx) {
    LodNode* node = &tree->nodes[node_idx];

    MeshGenJob* job = (MeshGenJob*)calloc(1, sizeof(MeshGenJob));
    job->tri = node->tri;
    job->depth = node->depth;
    job->planet_radius = tree->planet_radius;
    job->layer_thickness = tree->layer_thickness;
    job->sea_level = tree->sea_level;
    job->seed = tree->seed;
    job->origin[0] = (float)tree->world_origin[0];
    job->origin[1] = (float)tree->world_origin[1];
    job->origin[2] = (float)tree->world_origin[2];
    job->hex_frame_origin = tree->hex_frame_origin;
    job->hex_frame_east = tree->hex_frame_east;
    job->hex_frame_north = tree->hex_frame_north;
    job->node_index = node_idx;
    job->vertices = NULL;
    job->vertex_count = 0;
    job->completed = 0;

    // Non-blocking submit: if queue is full, free the job and let caller retry next frame
    if (!job_system_try_submit(tree->jobs, mesh_gen_worker, job)) {
        free(job);
        return false;  // Node stays in current state, will retry next frame
    }

    node->pending_job = job;
    node->state = LOD_GENERATING;
    return true;
}

static void free_node_mesh(LodNode* node) {
    if (node->cpu_vertices) {
        free(node->cpu_vertices);
        node->cpu_vertices = NULL;
    }
    node->cpu_vertex_count = 0;
}

// ---- Process completed async jobs ----

// Called from main thread: check for completed jobs and transfer results to nodes
static void process_completed_jobs(LodTree* tree) {
    for (int i = 0; i < tree->node_count; i++) {
        LodNode* node = &tree->nodes[i];
        if (node->state != LOD_GENERATING) continue;
        if (!node->pending_job) continue;

        MeshGenJob* job = (MeshGenJob*)node->pending_job;
        if (!job->completed) continue;

        // Check if job's origin matches current tree origin (may be stale after recenter)
        float cur_origin[3] = {(float)tree->world_origin[0],
                               (float)tree->world_origin[1],
                               (float)tree->world_origin[2]};
        bool stale = (job->origin[0] != cur_origin[0] ||
                      job->origin[1] != cur_origin[1] ||
                      job->origin[2] != cur_origin[2]);

        if (stale) {
            // Discard stale mesh and reset node for regeneration
            if (job->vertices) free(job->vertices);
            free(job);
            node->pending_job = NULL;
            node->state = LOD_UNLOADED;
            continue;
        }

        // Transfer mesh data from job to node
        node->cpu_vertices = job->vertices;
        node->cpu_vertex_count = job->vertex_count;
        node->cpu_vertex_stride = job->vertex_stride;
        node->is_hex_mesh = job->is_hex_mesh;
        job->vertices = NULL;  // Transfer ownership

        node->state = LOD_READY;
        free(job);
        node->pending_job = NULL;
    }
}

// ---- Split and merge ----

// Split a leaf node into 4 children
static bool split_node(LodTree* tree, int node_idx) {
    LodNode* node = &tree->nodes[node_idx];
    if (!node->is_leaf) return false;
    if (node->depth >= LOD_MAX_DEPTH) return false;
    if (tree->node_count + LOD_CHILDREN > tree->node_capacity &&
        tree->node_capacity >= LOD_MAX_NODES) {
        return false;  // At max capacity
    }

    SphericalTriangle child_tris[LOD_CHILDREN];
    subdivide_triangle(&node->tri, child_tris);

    for (int i = 0; i < LOD_CHILDREN; i++) {
        int child_idx = allocate_node(tree);
        if (child_idx < 0) {
            // Rollback: can't allocate all children
            for (int j = 0; j < i; j++) {
                tree->nodes[node->children[j]].state = LOD_UNLOADED;
                node->children[j] = -1;
            }
            return false;
        }

        // IMPORTANT: re-fetch node pointer after potential realloc in allocate_node
        node = &tree->nodes[node_idx];

        node->children[i] = child_idx;

        LodNode* child = &tree->nodes[child_idx];
        child->tri = child_tris[i];
        child->depth = node->depth + 1;
        child->parent = node_idx;
        child->is_leaf = true;
        child->state = LOD_UNLOADED;
    }

    node->is_leaf = false;
    return true;
}

// Merge a node's children back (remove children, make parent a leaf again)
static void merge_node(LodTree* tree, int node_idx) {
    LodNode* node = &tree->nodes[node_idx];
    if (node->is_leaf) return;

    // Check that all children are leaves (can't merge if grandchildren exist)
    // Also skip merge if any child is still generating its mesh
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int child_idx = node->children[i];
        if (child_idx < 0) continue;
        if (!tree->nodes[child_idx].is_leaf) return;
        if (tree->nodes[child_idx].state == LOD_GENERATING) return;
    }

    // Free children's resources and mark them as dead (reusable by allocate_node)
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int child_idx = node->children[i];
        if (child_idx < 0) continue;

        LodNode* child = &tree->nodes[child_idx];

        // Clean up pending job if any
        if (child->pending_job) {
            MeshGenJob* job = (MeshGenJob*)child->pending_job;
            if (job->vertices) free(job->vertices);
            free(job);
            child->pending_job = NULL;
        }

        free_node_mesh(child);
        if (child->gpu_buffer.id != SG_INVALID_ID) {
            sg_destroy_buffer(child->gpu_buffer);
            child->gpu_buffer.id = SG_INVALID_ID;
        }
        child->state = LOD_UNLOADED;
        child->parent = -1;  // Mark as dead — allocate_node can now reclaim this slot

        node->children[i] = -1;
    }

    node->is_leaf = true;
}

// ---- Distance-based LOD decision ----

// Returns effective distance from camera to nearest patch edge for LOD decisions.
// On the ground: pure arc distance (splits correctly near the player).
// Above max terrain: adds altitude penalty to prevent over-splitting from orbit.
static float patch_distance(const LodTree* tree, const LodNode* node) {
    float cam_r = sqrtf(vec3_dot(tree->camera_pos, tree->camera_pos));
    if (cam_r < 1.0f) cam_r = 1.0f;
    HMM_Vec3 cam_dir = vec3_scale(tree->camera_pos, 1.0f / cam_r);
    float cos_angle = vec3_dot(cam_dir, node->tri.center);
    if (cos_angle > 1.0f) cos_angle = 1.0f;
    if (cos_angle < -1.0f) cos_angle = -1.0f;
    float angle_to_center = acosf(cos_angle);
    float edge_angle = angle_to_center - node->tri.angular_radius;
    if (edge_angle < 0.0f) edge_angle = 0.0f;

    // Arc distance along the surface to nearest patch edge
    float arc_dist = edge_angle * tree->planet_radius;

    // Altitude above the highest possible terrain (prevents over-splitting from orbit)
    float max_surface_r = tree->planet_radius + TERRAIN_AMPLITUDE_M;
    float altitude = cam_r - max_surface_r;
    if (altitude < 0.0f) altitude = 0.0f;

    // On/below max terrain: pure arc distance. Above: altitude adds to distance.
    return sqrtf(arc_dist * arc_dist + altitude * altitude);
}


// Distance from camera to patch CENTER (not edge). Used for split/merge decisions
// so that all patches at the same depth get consistent treatment regardless of
// aperture-4 subdivision asymmetry (corner vs center children have different angular_radius).
static float patch_center_distance(const LodTree* tree, const LodNode* node) {
    float cam_r = sqrtf(vec3_dot(tree->camera_pos, tree->camera_pos));
    if (cam_r < 1.0f) cam_r = 1.0f;
    HMM_Vec3 cam_dir = vec3_scale(tree->camera_pos, 1.0f / cam_r);
    float cos_angle = vec3_dot(cam_dir, node->tri.center);
    if (cos_angle > 1.0f) cos_angle = 1.0f;
    if (cos_angle < -1.0f) cos_angle = -1.0f;
    float angle_to_center = acosf(cos_angle);

    float arc_dist = angle_to_center * tree->planet_radius;

    float max_surface_r = tree->planet_radius + TERRAIN_AMPLITUDE_M;
    float altitude = cam_r - max_surface_r;
    if (altitude < 0.0f) altitude = 0.0f;

    return sqrtf(arc_dist * arc_dist + altitude * altitude);
}

// Split when camera is close relative to the DEPTH-BASED arc length.
// Using depth_arc[] (precomputed) ensures all patches at the same depth
// have identical thresholds — no asymmetry from aperture-4 corner vs center children.
static bool should_split(const LodTree* tree, const LodNode* node) {
    if (node->depth >= LOD_MAX_DEPTH) return false;
    // Don't split into hex-detail depth when hex terrain covers the area
    if (node->depth + 1 >= LOD_MAX_DEPTH && tree->suppress_range > 0.0f) return false;
    float dist = patch_center_distance(tree, node);
    float patch_arc = tree->depth_arc[node->depth];
    return dist < patch_arc * tree->split_factor;
}

// Check if a patch's children should be merged back.
static bool should_merge(const LodTree* tree, const LodNode* node) {
    float dist = patch_center_distance(tree, node);
    float patch_arc = tree->depth_arc[node->depth];
    return dist > patch_arc * tree->split_factor * 2.0f;
}

// Back-hemisphere culling: skip patches entirely on the far side of the planet.
static bool patch_on_back_hemisphere(const LodTree* tree, const LodNode* node) {
    float cam_r = sqrtf(vec3_dot(tree->camera_pos, tree->camera_pos));
    if (cam_r < 1.0f) return false;
    HMM_Vec3 cam_dir = vec3_scale(tree->camera_pos, 1.0f / cam_r);
    // Horizon angle from camera
    float horizon_angle = acosf(fminf(1.0f, tree->planet_radius / cam_r));
    float cos_angle = vec3_dot(cam_dir, node->tri.center);
    if (cos_angle > 1.0f) cos_angle = 1.0f;
    if (cos_angle < -1.0f) cos_angle = -1.0f;
    float angle_to_center = acosf(cos_angle);
    // Patch edge can be closer than center
    float nearest_angle = angle_to_center - node->tri.angular_radius;
    // Include some margin past the horizon (for tall terrain + skirts)
    return nearest_angle > ((float)M_PI / 2.0f + horizon_angle + 0.1f);
}

// ---- Update: distance-based split/merge, mesh at all leaf levels ----

static void update_node(LodTree* tree, int node_idx) {
    LodNode* node = &tree->nodes[node_idx];

    // Skip patches completely behind the planet
    if (patch_on_back_hemisphere(tree, node)) {
        if (!node->is_leaf) {
            merge_node(tree, node_idx);
            // If merge failed (grandchildren exist), recurse so they can merge first
            node = &tree->nodes[node_idx];
            if (!node->is_leaf) {
                for (int i = 0; i < LOD_CHILDREN; i++) {
                    int ci = node->children[i];
                    if (ci >= 0) update_node(tree, ci);
                }
                // Retry merge now that children may have been collapsed
                merge_node(tree, node_idx);
            }
        }
        return;
    }

    if (node->is_leaf) {
        if (should_split(tree, node)) {
            // For depths >= 4, require a renderable mesh before splitting.
            // This guarantees parent-as-fallback: while children load,
            // the parent's mesh fills in seamlessly (no popping).
            // Depths 0-3 are too coarse to be individually visible, so they
            // cascade-split instantly to build the tree skeleton quickly.
            bool need_mesh_first = (node->depth >= 4);

            if (need_mesh_first) {
                if (node->state == LOD_UNLOADED) {
                    if (tree->jobs) {
                        submit_mesh_job(tree, node_idx);
                    } else {
                        generate_node_mesh(tree, node);
                    }
                    return;  // Wait for mesh before splitting
                }
                if (node->state == LOD_GENERATING || node->state == LOD_READY) {
                    return;  // Still loading/uploading, wait
                }
            }
            // Split (immediately for coarse depths, after mesh ready for fine)
            if (tree->splits_this_frame < LOD_MAX_SPLITS_PER_FRAME) {
                // Generate fallback mesh for coarse nodes before splitting
                // (start job but don't wait — best effort fallback)
                if (!need_mesh_first && node->state == LOD_UNLOADED) {
                    if (tree->jobs) {
                        submit_mesh_job(tree, node_idx);
                    }
                }
                split_node(tree, node_idx);
                tree->splits_this_frame++;
                // Recurse into new children immediately
                if (!tree->nodes[node_idx].is_leaf) {
                    for (int i = 0; i < LOD_CHILDREN; i++) {
                        int ci = tree->nodes[node_idx].children[i];
                        if (ci >= 0) {
                            update_node(tree, ci);
                        }
                    }
                }
            }
        } else {
            // This leaf is at the right detail level — generate mesh if needed
            if (node->state == LOD_UNLOADED) {
                if (tree->jobs) {
                    submit_mesh_job(tree, node_idx);
                } else {
                    generate_node_mesh(tree, node);
                }
            }
        }
    } else {
        // Internal node
        if (should_merge(tree, node)) {
            // Too far for this detail level: merge children
            merge_node(tree, node_idx);
            node = &tree->nodes[node_idx];
            if (!node->is_leaf) {
                // Merge failed (grandchildren exist) — recurse so they can merge first
                for (int i = 0; i < LOD_CHILDREN; i++) {
                    int ci = node->children[i];
                    if (ci >= 0) update_node(tree, ci);
                }
                // Retry merge now that children may have been collapsed
                merge_node(tree, node_idx);
                node = &tree->nodes[node_idx];
            }
            // After merge, this is a leaf — generate its mesh
            if (node->is_leaf && node->state == LOD_UNLOADED) {
                if (tree->jobs) {
                    submit_mesh_job(tree, node_idx);
                } else {
                    generate_node_mesh(tree, node);
                }
            }
        } else {
            // Keep children — recurse
            for (int i = 0; i < LOD_CHILDREN; i++) {
                int ci = tree->nodes[node_idx].children[i];
                if (ci >= 0) {
                    update_node(tree, ci);
                }
            }
        }
    }
}

void lod_tree_update(LodTree* tree, HMM_Vec3 camera_pos, HMM_Mat4 view_proj) {
    (void)view_proj;
    tree->camera_pos = camera_pos;
    tree->splits_this_frame = 0;

    // Update shared hex tangent frame.
    // All depth-13 patches use this single frame so their hex grids tile seamlessly.
    // Re-anchor when camera drifts ~2km from current frame center (well beyond the
    // ~1.3km depth-13 visible range, so old patches have naturally cycled out by then).
    {
        float cam_r = sqrtf(vec3_dot(camera_pos, camera_pos));
        if (cam_r < 1.0f) cam_r = 1.0f;
        HMM_Vec3 cam_dir = vec3_scale(camera_pos, 1.0f / cam_r);

        bool need_reanchor = !tree->hex_frame_valid;
        if (!need_reanchor) {
            float d = vec3_dot(cam_dir, tree->hex_frame_origin);
            // cos(angle) < 0.999997 ≈ 0.14° ≈ 2km at 800km radius
            if (d < 0.999997f) need_reanchor = true;
        }
        if (need_reanchor) {
            tree->hex_frame_origin = cam_dir;
            HMM_Vec3 wy = {{0.0f, 1.0f, 0.0f}};
            if (fabsf(vec3_dot(cam_dir, wy)) > 0.99f)
                wy = (HMM_Vec3){{1.0f, 0.0f, 0.0f}};
            tree->hex_frame_east = vec3_normalize(vec3_cross(wy, cam_dir));
            tree->hex_frame_north = vec3_normalize(vec3_cross(cam_dir, tree->hex_frame_east));
            tree->hex_frame_valid = true;
        }
    }

    // Process completed async jobs first (transfers mesh data to nodes)
    process_completed_jobs(tree);

    // Backpressure: if job queue is heavily loaded, skip tree updates this frame.
    // This prevents the main thread from queuing faster than workers can process,
    // which would eventually fill the queue and stall the frame loop.
    int pending = job_system_pending(tree->jobs);
    if (pending > 128) {  // Let workers catch up but don't starve the tree build
        // Queue heavily loaded — skip splits this frame, let workers catch up
        goto collect_stats;
    }

    // Sort root nodes by distance (nearest first) so closer faces get
    // priority for the per-frame split budget.
    int sorted_roots[LOD_ROOT_COUNT];
    float root_dists[LOD_ROOT_COUNT];
    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        sorted_roots[i] = i;
        root_dists[i] = patch_center_distance(tree, &tree->nodes[tree->root_nodes[i]]);
    }
    // Insertion sort (N=20, always fast)
    for (int i = 1; i < LOD_ROOT_COUNT; i++) {
        int key_idx = sorted_roots[i];
        float key_dist = root_dists[i];
        int j = i - 1;
        while (j >= 0 && root_dists[j] > key_dist) {
            sorted_roots[j + 1] = sorted_roots[j];
            root_dists[j + 1] = root_dists[j];
            j--;
        }
        sorted_roots[j + 1] = key_idx;
        root_dists[j + 1] = key_dist;
    }

    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        update_node(tree, tree->root_nodes[sorted_roots[i]]);
    }

collect_stats:
    // Reset per-level stats
    for (int d = 0; d <= LOD_MAX_DEPTH; d++) {
        tree->level_stats[d].patch_count = 0;
        tree->level_stats[d].vertex_count = 0;
        tree->level_stats[d].min_distance = 1e30f;
        tree->level_stats[d].max_distance = 0.0f;
    }

    // Collect stats from all nodes
    int leaves = 0;
    int total_verts = 0;
    for (int i = 0; i < tree->node_count; i++) {
        LodNode* node = &tree->nodes[i];
        if (!node->is_leaf) continue;
        if (node->state != LOD_READY && node->state != LOD_ACTIVE) continue;

        leaves++;
        int d = node->depth;
        if (d < 0) d = 0;
        if (d > LOD_MAX_DEPTH) d = LOD_MAX_DEPTH;

        tree->level_stats[d].patch_count++;

        if (node->state == LOD_ACTIVE) {
            total_verts += node->gpu_vertex_count;
            tree->level_stats[d].vertex_count += node->gpu_vertex_count;
        }

        // Compute center distance (same metric used for LOD split/merge decisions)
        float dist = patch_center_distance(tree, node);

        if (dist < tree->level_stats[d].min_distance)
            tree->level_stats[d].min_distance = dist;
        if (dist > tree->level_stats[d].max_distance)
            tree->level_stats[d].max_distance = dist;
    }
    tree->active_leaf_count = leaves;
    tree->total_vertex_count = total_verts;

    // Periodic console stats (every 120 frames)
    tree->stats_frame_counter++;
    if (tree->stats_frame_counter % 120 == 0) {
        LOG_VERBOSE(LOD, "STATS leaves=%d  total_verts=%d  nodes=%d/%d\n",
            leaves, total_verts, tree->node_count, tree->node_capacity);
        RAW_LOG_VERBOSE(LOD, "  Depth | Patches |  Vertices | Dist Min (km) | Dist Max (km)\n");
        RAW_LOG_VERBOSE(LOD, "  ------+---------+-----------+---------------+--------------\n");
        for (int d = 0; d <= LOD_MAX_DEPTH; d++) {
            if (tree->level_stats[d].patch_count == 0) continue;
            RAW_LOG_VERBOSE(LOD, "  %5d | %7d | %9d | %13.1f | %13.1f\n",
                d,
                tree->level_stats[d].patch_count,
                tree->level_stats[d].vertex_count,
                tree->level_stats[d].min_distance / 1000.0f,
                tree->level_stats[d].max_distance / 1000.0f);
        }
        if (LOG_ENABLE_LOD && log_verbose) fflush(stdout);
    }
}

// ---- GPU upload ----

void lod_tree_upload_meshes(LodTree* tree) {
    int uploads = 0;
    for (int i = 0; i < tree->node_count && uploads < LOD_MAX_UPLOADS_PER_FRAME; i++) {
        LodNode* node = &tree->nodes[i];
        if (node->state != LOD_READY) continue;
        if (!node->cpu_vertices || node->cpu_vertex_count == 0) {
            node->state = LOD_ACTIVE;
            node->gpu_vertex_count = 0;
            continue;
        }

        // Destroy old buffer if it exists
        if (node->gpu_buffer.id != SG_INVALID_ID) {
            sg_destroy_buffer(node->gpu_buffer);
        }

        int stride = node->cpu_vertex_stride > 0 ? node->cpu_vertex_stride : (int)sizeof(LodVertex);
        node->gpu_buffer = sg_make_buffer(&(sg_buffer_desc){
            .data = (sg_range){
                node->cpu_vertices,
                (size_t)node->cpu_vertex_count * stride
            },
            .label = node->is_hex_mesh ? "lod-hex-vertices" : "lod-patch-vertices",
        });

        // Check if buffer creation succeeded (pool may be exhausted)
        if (node->gpu_buffer.id == SG_INVALID_ID) {
            // Keep in READY state so we can retry next frame
            continue;
        }

        node->gpu_vertex_count = node->cpu_vertex_count;
        node->state = LOD_ACTIVE;
        uploads++;

        // Free CPU-side data after upload
        free(node->cpu_vertices);
        node->cpu_vertices = NULL;
    }
}

// ---- Rendering ----

// Internal render state passed through recursion
typedef struct {
    sg_bindings bind;
    sg_bindings hex_bind;
    sg_pipeline planet_pip;
    const LodUniformBlock* planet_ub;  // [2]: vs + fs
    const LodHexRenderInfo* hex;
    bool hex_active;         // Currently bound pipeline is hex
    int total_verts;
    LodPreDrawFunc pre_draw;
    void* user_data;
} RenderState;

// Draw a single node with the correct pipeline
static void draw_node(const LodTree* tree, const LodNode* node, RenderState* rs) {
    if (node->state != LOD_ACTIVE || node->gpu_vertex_count == 0 ||
        node->gpu_buffer.id == SG_INVALID_ID)
        return;

    // Suppress LOD hex-detail patches when hex terrain system is active
    if (node->is_hex_mesh && tree->suppress_range > 0.0f)
        return;

    if (node->is_hex_mesh && rs->hex) {
        // Switch to hex pipeline if not already active
        if (!rs->hex_active) {
            sg_apply_pipeline(rs->hex->pip);
            sg_apply_uniforms(rs->hex->vs_ub.slot,
                &(sg_range){rs->hex->vs_ub.data, rs->hex->vs_ub.size});
            sg_apply_uniforms(rs->hex->fs_ub.slot,
                &(sg_range){rs->hex->fs_ub.data, rs->hex->fs_ub.size});
            rs->hex_active = true;
        }
        // Skip pre_draw for hex nodes (pre_draw applies planet-specific uniforms)
        rs->hex_bind.vertex_buffers[0] = node->gpu_buffer;
        sg_apply_bindings(&rs->hex_bind);
    } else {
        // Switch back to planet pipeline if coming from hex
        if (rs->hex_active) {
            sg_apply_pipeline(rs->planet_pip);
            sg_apply_uniforms(rs->planet_ub[0].slot,
                &(sg_range){rs->planet_ub[0].data, rs->planet_ub[0].size});
            sg_apply_uniforms(rs->planet_ub[1].slot,
                &(sg_range){rs->planet_ub[1].data, rs->planet_ub[1].size});
            rs->hex_active = false;
        }
        if (rs->pre_draw) rs->pre_draw(node->depth, rs->user_data);
        rs->bind.vertex_buffers[0] = node->gpu_buffer;
        sg_apply_bindings(&rs->bind);
    }
    sg_draw(0, node->gpu_vertex_count, 1);
    rs->total_verts += node->gpu_vertex_count;
}

// Recursive render: draws ALL active leaf nodes at their respective LOD levels.
// Internal nodes with active meshes are used as fallback while children load.
static void render_node(const LodTree* tree, int node_idx, RenderState* rs) {
    const LodNode* node = &tree->nodes[node_idx];

    if (node->is_leaf) {
        draw_node(tree, node, rs);
        return;
    }

    // Internal node: check if all children have renderable meshes.
    bool all_children_ready = true;
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int ci = node->children[i];
        if (ci < 0) { all_children_ready = false; break; }
        const LodNode* child = &tree->nodes[ci];
        if (child->is_leaf) {
            if (child->state != LOD_ACTIVE || child->gpu_vertex_count == 0 ||
                child->gpu_buffer.id == SG_INVALID_ID) {
                all_children_ready = false;
                break;
            }
        }
    }

    if (all_children_ready) {
        for (int i = 0; i < LOD_CHILDREN; i++) {
            if (node->children[i] >= 0)
                render_node(tree, node->children[i], rs);
        }
    } else {
        if (node->state == LOD_ACTIVE && node->gpu_vertex_count > 0 &&
            node->gpu_buffer.id != SG_INVALID_ID) {
            draw_node(tree, node, rs);
        } else {
            for (int i = 0; i < LOD_CHILDREN; i++) {
                if (node->children[i] >= 0)
                    render_node(tree, node->children[i], rs);
            }
        }
    }
}

void lod_tree_render(LodTree* tree, sg_pipeline pip, HMM_Mat4 vp,
                     LodPreDrawFunc pre_draw, void* user_data,
                     const LodUniformBlock planet_ub[2],
                     const LodHexRenderInfo* hex) {
    (void)vp;

    RenderState rs = {0};
    rs.bind = (sg_bindings){0};
    rs.hex_bind = (sg_bindings){0};
    rs.planet_pip = pip;
    rs.planet_ub = planet_ub;
    if (hex) {
        rs.hex_bind.views[0] = hex->atlas_view;
        rs.hex_bind.samplers[0] = hex->atlas_smp;
    }
    rs.hex = hex;
    rs.hex_active = false;
    rs.total_verts = 0;
    rs.pre_draw = pre_draw;
    rs.user_data = user_data;

    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        render_node(tree, tree->root_nodes[i], &rs);
    }

    tree->total_vertex_count = rs.total_verts;
}

// ---- Wireframe overlay ----

typedef struct {
    sg_pipeline pip_planet;
    sg_pipeline pip_hex;
    sg_buffer index_buf;
    const LodUniformBlock* ub;
    int current_pip;  // 0=none, 1=planet, 2=hex
} WireState;

static void draw_node_wire(const LodNode* node, WireState* ws) {
    if (node->state != LOD_ACTIVE || node->gpu_vertex_count == 0 ||
        node->gpu_buffer.id == SG_INVALID_ID)
        return;

    int want = node->is_hex_mesh ? 2 : 1;
    if (want != ws->current_pip) {
        sg_apply_pipeline(want == 2 ? ws->pip_hex : ws->pip_planet);
        sg_apply_uniforms(ws->ub[0].slot,
            &(sg_range){ws->ub[0].data, ws->ub[0].size});
        sg_apply_uniforms(ws->ub[1].slot,
            &(sg_range){ws->ub[1].data, ws->ub[1].size});
        ws->current_pip = want;
    }

    sg_bindings bind = {0};
    bind.vertex_buffers[0] = node->gpu_buffer;
    bind.index_buffer = ws->index_buf;
    sg_apply_bindings(&bind);
    sg_draw(0, node->gpu_vertex_count * 2, 1);
}

static void render_node_wire(const LodTree* tree, int node_idx, WireState* ws) {
    const LodNode* node = &tree->nodes[node_idx];

    if (node->is_leaf) {
        draw_node_wire(node, ws);
        return;
    }

    bool all_ready = true;
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int ci = node->children[i];
        if (ci < 0) { all_ready = false; break; }
        const LodNode* child = &tree->nodes[ci];
        if (child->is_leaf &&
            (child->state != LOD_ACTIVE || child->gpu_vertex_count == 0 ||
             child->gpu_buffer.id == SG_INVALID_ID)) {
            all_ready = false;
            break;
        }
    }

    if (all_ready) {
        for (int i = 0; i < LOD_CHILDREN; i++)
            if (node->children[i] >= 0)
                render_node_wire(tree, node->children[i], ws);
    } else if (node->state == LOD_ACTIVE && node->gpu_vertex_count > 0 &&
               node->gpu_buffer.id != SG_INVALID_ID) {
        draw_node_wire(node, ws);
    } else {
        for (int i = 0; i < LOD_CHILDREN; i++)
            if (node->children[i] >= 0)
                render_node_wire(tree, node->children[i], ws);
    }
}

void lod_tree_render_wireframe(LodTree* tree,
                                sg_pipeline pip_planet, sg_pipeline pip_hex,
                                sg_buffer index_buf, const LodUniformBlock ub[2]) {
    WireState ws = {
        .pip_planet = pip_planet,
        .pip_hex = pip_hex,
        .index_buf = index_buf,
        .ub = ub,
        .current_pip = 0,
    };
    for (int i = 0; i < LOD_ROOT_COUNT; i++)
        render_node_wire(tree, tree->root_nodes[i], &ws);
}

// ---- Floating origin recentering ----

#define ORIGIN_RECENTER_THRESHOLD 5000000.0  // 5000km — well beyond single-planet range (r=796km, max surface distance ~2500km)

bool lod_tree_update_origin(LodTree* tree, const double cam_pos_d[3]) {
    double dx = cam_pos_d[0] - tree->world_origin[0];
    double dy = cam_pos_d[1] - tree->world_origin[1];
    double dz = cam_pos_d[2] - tree->world_origin[2];
    double dist_sq = dx*dx + dy*dy + dz*dz;

    if (dist_sq < ORIGIN_RECENTER_THRESHOLD * ORIGIN_RECENTER_THRESHOLD) {
        return false;  // Close enough, no recenter needed
    }

    // Recenter origin to camera position
    double new_origin[3] = {cam_pos_d[0], cam_pos_d[1], cam_pos_d[2]};

    LOG(GAME, "Floating origin recenter: (%.0f,%.0f,%.0f) -> (%.0f,%.0f,%.0f)\n",
        tree->world_origin[0], tree->world_origin[1], tree->world_origin[2],
        new_origin[0], new_origin[1], new_origin[2]);

    tree->world_origin[0] = new_origin[0];
    tree->world_origin[1] = new_origin[1];
    tree->world_origin[2] = new_origin[2];

    // Invalidate all node meshes — they need regeneration with new origin.
    // The LOD system will naturally rebuild everything within a few frames.
    for (int i = 0; i < tree->node_count; i++) {
        LodNode* node = &tree->nodes[i];

        // Cancel pending jobs
        if (node->pending_job) {
            MeshGenJob* job = (MeshGenJob*)node->pending_job;
            // Wait for completion (can't free while worker is using it)
            // The job will complete with old origin, but we'll re-submit
            // Just mark it — process_completed_jobs will pick it up
        }

        // Free CPU mesh data
        if (node->cpu_vertices) {
            free(node->cpu_vertices);
            node->cpu_vertices = NULL;
            node->cpu_vertex_count = 0;
        }

        // Destroy GPU buffer
        if (node->gpu_buffer.id != SG_INVALID_ID) {
            sg_destroy_buffer(node->gpu_buffer);
            node->gpu_buffer.id = SG_INVALID_ID;
            node->gpu_vertex_count = 0;
        }

        // Only reset state for non-generating nodes (let pending jobs finish)
        if (node->state != LOD_GENERATING) {
            node->state = LOD_UNLOADED;
        }
    }

    return true;
}

double lod_tree_terrain_height(const LodTree* tree, HMM_Vec3 world_pos) {
    HMM_Vec3 unit = vec3_normalize(world_pos);
    fnl_state continental = create_continental_noise(tree->seed);
    fnl_state mountain_n = create_mountain_noise(tree->seed);
    fnl_state warp = create_warp_noise(tree->seed);
    fnl_state detail = create_detail_noise(tree->seed);
    float h_m = sample_terrain_height_m(&continental, &mountain_n, &warp, &detail, unit);
    float effective_h_m = fmaxf(h_m, TERRAIN_SEA_LEVEL_M);

    // Add in double to avoid float quantization (6.25cm steps at 800km)
    return (double)tree->planet_radius + (double)effective_h_m;
}

int lod_tree_active_leaves(const LodTree* tree) {
    return tree->active_leaf_count;
}

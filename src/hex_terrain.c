#include "hex_terrain.h"
#include "hex_vertex.h"
#include "math_utils.h"
#include "planet.h"
#include "log_config.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#include "FastNoiseLite.h"
#include "terrain_noise.h"

// Debug logging to both console and file
#include <stdarg.h>
static void hex_log(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
    FILE* f = fopen("hex_debug.log", "a");
    if (f) {
        va_start(args, fmt);
        vfprintf(f, fmt, args);
        va_end(args);
        fclose(f);
    }
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---- Hex geometry constants ----
// Flat-topped hex: vertices at 0, 60, 120, 180, 240, 300 degrees
static const float HEX_COS[6] = {
    1.0f, 0.5f, -0.5f, -1.0f, -0.5f, 0.5f
};
static const float HEX_SIN[6] = {
    0.0f, 0.8660254f, 0.8660254f, 0.0f, -0.8660254f, -0.8660254f
};

// Axial neighbor offsets (flat-topped hex, even-q offset grid)
static const int AX_NEIGHBORS[6][2] = {
    { +1,  0 },   // direction 0: NE
    { +1, -1 },   // direction 1: SE
    {  0, -1 },   // direction 2: S
    { -1,  0 },   // direction 3: SW
    { -1, +1 },   // direction 4: NW
    {  0, +1 },   // direction 5: N
};

// Map neighbor direction to the hex edge index that faces that neighbor.
static const int DIR_TO_EDGE[6] = {
    0,  // dir 0 (NE neighbor) -> edge 0
    5,  // dir 1 (SE neighbor) -> edge 5
    4,  // dir 2 (S  neighbor) -> edge 4
    3,  // dir 3 (SW neighbor) -> edge 3
    2,  // dir 4 (NW neighbor) -> edge 2
    1,  // dir 5 (N  neighbor) -> edge 1
};

// Convert offset (col, row) to axial (q, r)
static void offset_to_axial(int col, int row, int* q, int* r) {
    *q = col;
    *r = row - (col - (col & 1)) / 2;
}

// Convert axial (q, r) to offset (col, row)
static void axial_to_offset(int q, int r, int* col, int* row) {
    *col = q;
    *row = r + (q - (q & 1)) / 2;
}

// Get neighbor offset coordinates for hex at (col, row) in direction dir (0-5)
static void hex_neighbor(int col, int row, int dir, int* ncol, int* nrow) {
    int q, r;
    offset_to_axial(col, row, &q, &r);
    int nq = q + AX_NEIGHBORS[dir][0];
    int nr = r + AX_NEIGHBORS[dir][1];
    axial_to_offset(nq, nr, ncol, nrow);
}

// ---- Terrain noise (shared with lod.c via terrain_noise.h) ----

// ---- Stone Arch (same as lod.c, for voxel filling) ----
#define ARCH_R_OUTER    100.0f   // outer semicircle radius (total width 200m, height 100m)
#define ARCH_R_INNER     80.0f   // inner semicircle radius (opening)
#define ARCH_HALF_D      25.0f   // half-depth = 25m (total 50m)

typedef struct {
    HMM_Vec3 center;     // unit sphere center
    HMM_Vec3 east, north; // tangent plane axes
} HtArchFrame;

static HtArchFrame ht_arch_make_frame(void) {
    HtArchFrame f;
    // Twilight zone spawn (lat 68.55, lon -106.52)
    f.center = vec3_normalize((HMM_Vec3){{-0.103983f, 0.930767f, -0.350514f}});
    HMM_Vec3 up = f.center;
    HMM_Vec3 arb = fabsf(up.Y) < 0.999f ? (HMM_Vec3){{0,1,0}} : (HMM_Vec3){{1,0,0}};
    f.east  = vec3_normalize(vec3_cross(arb, up));
    f.north = vec3_cross(up, f.east);
    return f;
}

static void ht_arch_project(const HtArchFrame* f, HMM_Vec3 unit_pos, float planet_radius,
                             float* lx, float* lz) {
    float d = vec3_dot(unit_pos, f->center);
    HMM_Vec3 off = vec3_sub(unit_pos, vec3_scale(f->center, d));
    *lx = vec3_dot(off, f->east)  * planet_radius;
    *lz = vec3_dot(off, f->north) * planet_radius;
}

// Query arch solid range at arch-local (lx, lz).
// Returns arch_bottom and arch_top in world layers. INT16_MIN = no arch.
static void ht_arch_query(float lx, float lz, float arch_base,
                           int* arch_bottom, int* arch_top) {
    if (fabsf(lz) > ARCH_HALF_D) {
        *arch_bottom = INT16_MIN; *arch_top = INT16_MIN; return;
    }
    float lx2 = lx * lx;
    if (lx2 > ARCH_R_OUTER * ARCH_R_OUTER) {
        *arch_bottom = INT16_MIN; *arch_top = INT16_MIN; return;
    }
    float outer_y = sqrtf(ARCH_R_OUTER * ARCH_R_OUTER - lx2);
    int top = (int)(arch_base + outer_y);
    int base = (int)arch_base;
    if (lx2 < ARCH_R_INNER * ARCH_R_INNER) {
        float inner_y = sqrtf(ARCH_R_INNER * ARCH_R_INNER - lx2);
        int bottom = (int)ceilf(arch_base + inner_y);
        if (bottom >= top) {
            *arch_bottom = INT16_MIN; *arch_top = INT16_MIN; return;
        }
        *arch_bottom = bottom;
        *arch_top = top;
    } else {
        *arch_bottom = base;
        *arch_top = top;
    }
}

static float ht_arch_compute_base(const HtArchFrame* f, int seed) {
    fnl_state continental = ht_create_continental_noise(seed);
    fnl_state mountain = ht_create_mountain_noise(seed);
    fnl_state warp = ht_create_warp_noise(seed);
    fnl_state detail = ht_create_detail_noise(seed);
    float h = ht_sample_height_m(&continental, &mountain, &warp, &detail, f->center);
    return fmaxf(h, TERRAIN_SEA_LEVEL_M) - 50.0f;
}

// ---- Public terrain height API (single source of truth for LOD + hex) ----

void hex_terrain_sample_height(int seed, float planet_radius, HMM_Vec3 unit_pos,
                                float* out_raw_h_m, float* out_effective_h_m) {
    fnl_state continental = ht_create_continental_noise(seed);
    fnl_state mountain = ht_create_mountain_noise(seed);
    fnl_state warp = ht_create_warp_noise(seed);
    fnl_state detail = ht_create_detail_noise(seed);
    float h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, unit_pos);
    float effective = fmaxf(h_m, TERRAIN_SEA_LEVEL_M);

    // Apply arch (and any future generated structures) as heightmap bump.
    // LOD can't represent the hollow underside, but the silhouette is correct at distance.
    HtArchFrame arch_f = ht_arch_make_frame();
    float arch_base = ht_arch_compute_base(&arch_f, seed);
    float alx, alz;
    ht_arch_project(&arch_f, unit_pos, planet_radius, &alx, &alz);
    if (fabsf(alz) <= ARCH_HALF_D) {
        float lx2 = alx * alx;
        if (lx2 <= ARCH_R_OUTER * ARCH_R_OUTER) {
            float outer_y = sqrtf(ARCH_R_OUTER * ARCH_R_OUTER - lx2);
            float arch_top = arch_base + outer_y;
            if (arch_top > effective) effective = arch_top;
        }
    }

    if (out_raw_h_m) *out_raw_h_m = h_m;
    if (out_effective_h_m) *out_effective_h_m = effective;
}

// Deterministic hash for dithering at biome transitions
static float dither_hash(int gcol, int grow, int seed) {
    uint32_t h = (uint32_t)(gcol * 374761393 + grow * 668265263 + seed * 1103515245);
    h = (h ^ (h >> 13)) * 1274126177u;
    h = h ^ (h >> 16);
    return (float)(h & 0xFFFFu) / 65535.0f;
}

// Voxel type by height with deterministic dithering at transitions.
// Each hex column gets a stable random threshold so biome boundaries
// show a natural mix of block types instead of a hard line.
static VoxelType ht_voxel_type(float height_m, int gcol, int grow, int seed) {
    float rel = height_m - TERRAIN_SEA_LEVEL_M;
    if (height_m < TERRAIN_SEA_LEVEL_M) return VOXEL_WATER;

    float d = dither_hash(gcol, grow, seed);
    float dz = 100.0f;  // dither zone half-width in meters

    // Sand -> Grass at 200m
    if (rel < 200.0f - dz) return VOXEL_SAND;
    if (rel < 200.0f + dz) {
        float t = (rel - (200.0f - dz)) / (2.0f * dz);
        return (d < t) ? VOXEL_GRASS : VOXEL_SAND;
    }
    // Grass -> Stone at 1500m
    if (rel < 1500.0f - dz) return VOXEL_GRASS;
    if (rel < 1500.0f + dz) {
        float t = (rel - (1500.0f - dz)) / (2.0f * dz);
        return (d < t) ? VOXEL_STONE : VOXEL_GRASS;
    }
    // Stone -> Ice at 4500m
    if (rel < 4500.0f - dz) return VOXEL_STONE;
    if (rel < 4500.0f + dz) {
        float t = (rel - (4500.0f - dz)) / (2.0f * dz);
        return (d < t) ? VOXEL_ICE : VOXEL_STONE;
    }
    return VOXEL_ICE;
}

// Voxel type at a specific layer depth below surface
static VoxelType ht_voxel_type_at_depth(VoxelType surface_type, int depth_layers) {
    if (depth_layers <= 0) return surface_type;
    if (depth_layers < 4) return VOXEL_DIRT;
    return VOXEL_STONE;
}

// ---- Terrain color (matches lod.c terrain_color_m for seamless LOD transition) ----
static HMM_Vec3 ht_terrain_color(float height_m) {
    float rel = height_m - TERRAIN_SEA_LEVEL_M;

    const HMM_Vec3 water_deep    = {{0.06f, 0.10f, 0.25f}};
    const HMM_Vec3 water_shallow = {{0.12f, 0.21f, 0.39f}};
    const HMM_Vec3 sand          = {{0.94f, 0.84f, 0.77f}};
    const HMM_Vec3 grass         = {{0.07f, 0.58f, 0.35f}};
    const HMM_Vec3 rock          = {{0.46f, 0.45f, 0.45f}};
    const HMM_Vec3 high_stone    = {{0.50f, 0.50f, 0.50f}};
    const HMM_Vec3 ice           = {{0.44f, 0.77f, 0.97f}};

    if (height_m < TERRAIN_SEA_LEVEL_M) {
        float depth_t = fminf(1.0f, -rel / 2000.0f);
        return vec3_lerp(water_shallow, water_deep, depth_t);
    }

    float blend = 100.0f;
    float t1 = ht_smoothstepf(200.0f  - blend, 200.0f  + blend, rel);
    float t2 = ht_smoothstepf(1500.0f - blend, 1500.0f + blend, rel);
    float t3 = ht_smoothstepf(3000.0f - blend, 3000.0f + blend, rel);
    float t4 = ht_smoothstepf(4500.0f - blend, 4500.0f + blend, rel);

    HMM_Vec3 c = sand;
    c = vec3_lerp(c, grass,      t1);
    c = vec3_lerp(c, rock,       t2);
    c = vec3_lerp(c, high_stone, t3);
    c = vec3_lerp(c, ice,        t4);
    return c;
}

// ---- Texture atlas mapping ----
#define HEX_ATLAS_TILES 10
#define ATLAS_WATER      0
#define ATLAS_SAND       1
#define ATLAS_DIRT       2
#define ATLAS_GRASS      3
#define ATLAS_STONE      4
#define ATLAS_ICE        5
#define ATLAS_SNOW       6
#define ATLAS_DIRT_GRASS  7
#define ATLAS_DIRT_SNOW   8
#define ATLAS_TORCH      9

// Transparent voxels: don't emit hex geometry, don't occlude neighbor faces
static inline bool is_transparent(uint8_t vtype) {
    return vtype == VOXEL_AIR || vtype == VOXEL_TORCH;
}

static int voxel_top_atlas(VoxelType type) {
    switch (type) {
        case VOXEL_WATER:   return ATLAS_WATER;
        case VOXEL_SAND:    return ATLAS_SAND;
        case VOXEL_DIRT:    return ATLAS_DIRT;
        case VOXEL_GRASS:   return ATLAS_GRASS;
        case VOXEL_STONE:   return ATLAS_STONE;
        case VOXEL_ICE:     return ATLAS_ICE;
        case VOXEL_BEDROCK: return ATLAS_STONE;
        case VOXEL_TORCH:   return ATLAS_TORCH;
        default:            return ATLAS_GRASS;
    }
}

static int voxel_side_atlas(VoxelType type) {
    switch (type) {
        case VOXEL_GRASS: return ATLAS_DIRT_GRASS;
        case VOXEL_ICE:   return ATLAS_DIRT_SNOW;
        default:          return voxel_top_atlas(type);
    }
}

static int voxel_bottom_atlas(VoxelType type) {
    switch (type) {
        case VOXEL_GRASS: return ATLAS_DIRT;
        default:          return voxel_top_atlas(type);
    }
}

// ---- Tangent frame computation ----

static void compute_tangent_frame(HexTerrain* ht, HMM_Vec3 camera_pos) {
    float cam_r = sqrtf(vec3_dot(camera_pos, camera_pos));
    if (cam_r < 1.0f) cam_r = 1.0f;
    HMM_Vec3 cam_dir = vec3_scale(camera_pos, 1.0f / cam_r);

    ht->tangent_up = cam_dir;
    float ground_r = ht->planet_radius + TERRAIN_SEA_LEVEL_M;
    ht->tangent_origin = vec3_scale(cam_dir, ground_r);

    HMM_Vec3 world_y = (HMM_Vec3){{0.0f, 1.0f, 0.0f}};
    if (fabsf(vec3_dot(cam_dir, world_y)) > 0.99f) {
        world_y = (HMM_Vec3){{1.0f, 0.0f, 0.0f}};
    }

    ht->tangent_east = vec3_normalize(vec3_cross(world_y, cam_dir));
    ht->tangent_north = vec3_normalize(vec3_cross(cam_dir, ht->tangent_east));
}

// ---- Hex position computation ----

static void hex_local_pos(int gcol, int grow, float* out_x, float* out_z) {
    *out_x = gcol * HEX_COL_SPACING;
    *out_z = ((gcol & 1) ? (grow + 0.5f) : (float)grow) * HEX_ROW_SPACING;
}

static HMM_Vec3 tangent_to_world(const HexTerrain* ht, float lx, float lz, float radius) {
    HMM_Vec3 world_dir = vec3_add(ht->tangent_origin,
        vec3_add(vec3_scale(ht->tangent_east, lx),
                 vec3_scale(ht->tangent_north, lz)));
    return vec3_scale(vec3_normalize(world_dir), radius);
}

// Inverse of tangent_to_world: world position → tangent-plane (lx, lz).
// Uses spherical projection to exactly invert the normalize() in the forward map.
// The flat projection (dot(P-T_o, T_e)) diverges at distance because it ignores
// the sphere curvature that normalize() applies.
static void world_to_tangent(const HexTerrain* ht, HMM_Vec3 world_pos,
                              float* out_lx, float* out_lz) {
    HMM_Vec3 dir = vec3_normalize(world_pos);
    float cos_theta = vec3_dot(dir, ht->tangent_up);
    if (cos_theta < 0.01f) {
        // Point is on the far side of the sphere — shouldn't happen for nearby hex
        *out_lx = 0.0f;
        *out_lz = 0.0f;
        return;
    }
    // ground_r = |tangent_origin| = planet_radius + sea_level
    float ground_r = ht->planet_radius + TERRAIN_SEA_LEVEL_M;
    float k = ground_r / cos_theta;
    *out_lx = k * vec3_dot(dir, ht->tangent_east);
    *out_lz = k * vec3_dot(dir, ht->tangent_north);
}

// Convert tangent-plane pixel position to hex grid coordinates
static void pixel_to_hex(float px, float pz, int* out_col, int* out_row) {
    float q_f = (2.0f / 3.0f) * px / HEX_RADIUS;
    float r_f = (-px / 3.0f + pz * 0.57735027f) / HEX_RADIUS;
    float s_f = -q_f - r_f;

    int qi = (int)roundf(q_f);
    int ri = (int)roundf(r_f);
    int si = (int)roundf(s_f);

    float dq = fabsf((float)qi - q_f);
    float dr = fabsf((float)ri - r_f);
    float ds = fabsf((float)si - s_f);

    if (dq > dr && dq > ds)      qi = -ri - si;
    else if (dr > ds)             ri = -qi - si;

    axial_to_offset(qi, ri, out_col, out_row);
}

// ---- Unit-sphere helpers (for voxel edit persistence) ----

// Convert hex grid (gcol, grow) to 64-bit unit-sphere direction.
// No polar singularity — doubles preserve full precision everywhere.
static void gcol_grow_to_unit(const HexTerrain* ht, int gcol, int grow,
                               double* out_ux, double* out_uy, double* out_uz) {
    float lx, lz;
    hex_local_pos(gcol, grow, &lx, &lz);
    HMM_Vec3 world_dir = vec3_normalize(
        vec3_add(ht->tangent_origin,
            vec3_add(vec3_scale(ht->tangent_east, lx),
                     vec3_scale(ht->tangent_north, lz))));
    // Promote to double for stable storage
    double dx = (double)world_dir.X;
    double dy = (double)world_dir.Y;
    double dz = (double)world_dir.Z;
    double len = sqrt(dx*dx + dy*dy + dz*dz);
    *out_ux = dx / len;
    *out_uy = dy / len;
    *out_uz = dz / len;
}

// ---- Mesh generation job ----

typedef struct HexMeshJob {
    int cx, cz;
    float planet_radius;
    int seed;

    HMM_Vec3 tangent_origin;
    HMM_Vec3 tangent_up;
    HMM_Vec3 tangent_east;
    HMM_Vec3 tangent_north;
    float origin[3];
    bool is_transition;

    // 3D voxel data (generated on worker thread)
    uint8_t* voxels;  // HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS (malloc'd)
    uint8_t* sky_map;   // Per-voxel sky light (0..SKY_MAX), same layout as voxels (malloc'd)
    uint8_t* torch_map; // Per-voxel torch light (0..TORCH_LIGHT_MAX), same layout (malloc'd)
    int base_layer;
    bool has_voxels;  // true = voxels already populated (remesh only, skip terrain gen)
    int16_t col_min_solid[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE];
    int16_t col_max_solid[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE];

    // Cached arch data (computed once per job)
    HtArchFrame arch_frame;
    float arch_base;

    // Output mesh
    HexVertex* vertices;
    int vertex_count;

    // Physics debug wireframe (hex prism outlines as float3 line pairs)
    float* wire_verts;
    int wire_count;    // number of floats (count/3 = vertex count, count/6 = line count)
    int wire_cap;

    int chunk_index;
    volatile int completed;
    int nb_count;  // number of neighbor voxel sets populated (0-8)
    bool is_initial_gen;  // true = first-time generation (no pre-existing voxels)

    // Edit cache (read-only pointer for edit queries on worker thread)
    const EditCache* edit_cache;

    // Neighbor chunk voxel data for cross-chunk BFS and geometry.
    // Indexed as [dx+1][dz+1] where dx,dz in {-1,0,1} (center [1][1] is self).
    // NULL if that neighbor doesn't have voxel data.
    uint8_t* nb_voxels[3][3];    // malloc'd copies of neighbor voxel arrays
    int nb_base_layer[3][3];     // neighbor base_layer values

    // Output: boundary torch light values (extracted after BFS for cross-chunk propagation)
    // [0]=col 0, [1]=col 31, [2]=row 0, [3]=row 31
    // Each is HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS bytes. NULL if all zeros.
    uint8_t* out_border_torch[4];

    // Input: neighbor boundary torch light values (facing our chunk edge)
    // [dx+1][dz+1] = pointer to the neighbor's facing border. NULL if unavailable.
    // Only populated for cardinal neighbors (dx=0 or dz=0, not both nonzero).
    uint8_t* nb_border_torch[3][3];
} HexMeshJob;

// Free all neighbor voxel data and border torch data from a mesh job
static void free_job_neighbors(HexMeshJob* job) {
    for (int dx = 0; dx < 3; dx++)
        for (int dz = 0; dz < 3; dz++) {
            if (dx == 1 && dz == 1) continue;
            if (job->nb_voxels[dx][dz]) {
                free(job->nb_voxels[dx][dz]);
                job->nb_voxels[dx][dz] = NULL;
            }
            if (job->nb_border_torch[dx][dz]) {
                free(job->nb_border_torch[dx][dz]);
                job->nb_border_torch[dx][dz] = NULL;
            }
        }
}

// Read a voxel from a neighbor chunk's data.
// gcol/grow are GLOBAL hex coords. Returns VOXEL_AIR if neighbor data not available.
static uint8_t job_read_neighbor_voxel(const HexMeshJob* job, int gcol, int grow, int layer) {
    int ncx = (gcol < 0) ? -(((-gcol - 1) / HEX_CHUNK_SIZE) + 1) : (gcol / HEX_CHUNK_SIZE);
    int ncz = (grow < 0) ? -(((-grow - 1) / HEX_CHUNK_SIZE) + 1) : (grow / HEX_CHUNK_SIZE);
    int dx = ncx - job->cx + 1;
    int dz = ncz - job->cz + 1;
    if (dx < 0 || dx > 2 || dz < 0 || dz > 2) return VOXEL_AIR;
    const uint8_t* nb = job->nb_voxels[dx][dz];
    if (!nb) return VOXEL_AIR;  // neighbor data not available

    int local_col = gcol - ncx * HEX_CHUNK_SIZE;
    int local_row = grow - ncz * HEX_CHUNK_SIZE;
    if (local_col < 0 || local_col >= HEX_CHUNK_SIZE ||
        local_row < 0 || local_row >= HEX_CHUNK_SIZE)
        return VOXEL_AIR;

    int nb_base = job->nb_base_layer[dx][dz];
    int local_layer = layer - nb_base;
    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) return VOXEL_AIR;
    return HEX_VOXEL(nb, local_col, local_row, local_layer);
}

// ---- Unit-sphere → hex grid under job's tangent frame ----

static void unit_to_gcol_grow(const HexMeshJob* job, double ux, double uy, double uz,
                               int* out_gcol, int* out_grow) {
    // Use float tangent-plane projection (sufficient for local hex lookup)
    HMM_Vec3 unit_pos = {{(float)ux, (float)uy, (float)uz}};

    HMM_Vec3 tangent_up = vec3_normalize(job->tangent_origin);
    float cos_theta = vec3_dot(unit_pos, tangent_up);
    if (cos_theta < 0.01f) { *out_gcol = 0; *out_grow = 0; return; }
    float ground_r = sqrtf(vec3_dot(job->tangent_origin, job->tangent_origin));
    float k = ground_r / cos_theta;
    float lx = k * vec3_dot(unit_pos, job->tangent_east);
    float lz = k * vec3_dot(unit_pos, job->tangent_north);

    pixel_to_hex(lx, lz, out_gcol, out_grow);
}

// Compute cube-map sector key range for a chunk (for edit queries).
static void compute_chunk_sector_bounds(const HexMeshJob* job,
                                         int8_t* out_face,
                                         int16_t* u_min, int16_t* u_max,
                                         int16_t* v_min, int16_t* v_max) {
    // Compute chunk center unit vector for cube-map face determination
    int center_gcol = job->cx * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE / 2;
    int center_grow = job->cz * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE / 2;
    float clx, clz;
    hex_local_pos(center_gcol, center_grow, &clx, &clz);
    HMM_Vec3 center_dir = vec3_normalize(
        vec3_add(job->tangent_origin,
            vec3_add(vec3_scale(job->tangent_east, clx),
                     vec3_scale(job->tangent_north, clz))));

    int8_t face;
    float center_u, center_v;
    edit_unit_to_cubemap(center_dir, &face, &center_u, &center_v);
    *out_face = face;

    // Project all 4 corners to cube-map UV on the same face
    float uv_min_u = center_u, uv_max_u = center_u;
    float uv_min_v = center_v, uv_max_v = center_v;

    int corners[4][2] = { {0, 0}, {HEX_CHUNK_SIZE - 1, 0},
                          {0, HEX_CHUNK_SIZE - 1}, {HEX_CHUNK_SIZE - 1, HEX_CHUNK_SIZE - 1} };
    for (int c = 0; c < 4; c++) {
        int gcol = job->cx * HEX_CHUNK_SIZE + corners[c][0];
        int grow = job->cz * HEX_CHUNK_SIZE + corners[c][1];
        float lx, lz;
        hex_local_pos(gcol, grow, &lx, &lz);
        HMM_Vec3 dir = vec3_normalize(
            vec3_add(job->tangent_origin,
                vec3_add(vec3_scale(job->tangent_east, lx),
                         vec3_scale(job->tangent_north, lz))));

        int8_t cf;
        float cu, cv;
        edit_unit_to_cubemap(dir, &cf, &cu, &cv);
        // If corner is on same face, use its UV; otherwise ignore (chunk is tiny)
        if (cf == face) {
            if (cu < uv_min_u) uv_min_u = cu;
            if (cu > uv_max_u) uv_max_u = cu;
            if (cv < uv_min_v) uv_min_v = cv;
            if (cv > uv_max_v) uv_max_v = cv;
        }
    }

    // Pad by 2 hex spacings to catch edits at chunk edges
    float pad = 2.0f * HEX_ROW_SPACING / job->planet_radius;
    uv_min_u -= pad; uv_max_u += pad;
    uv_min_v -= pad; uv_max_v += pad;

    // Convert UV to sector indices
    float sector_uv = EDIT_SECTOR_SIZE_M / job->planet_radius;
    *u_min = (int16_t)floorf(uv_min_u / sector_uv);
    *u_max = (int16_t)floorf(uv_max_u / sector_uv);
    *v_min = (int16_t)floorf(uv_min_v / sector_uv);
    *v_max = (int16_t)floorf(uv_max_v / sector_uv);
}

// ---- Emit triangle helper (textured) ----

typedef struct { float u, v; } HexUV;

static void hex_emit_tri(HexVertex** verts, int* count, int* cap,
                          HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                          HMM_Vec3 winding_normal,
                          HMM_Vec3 shading_normal,
                          HexUV uv0, HexUV uv1, HexUV uv2,
                          HMM_Vec3 vtx_color, float sky_light,
                          float torch_light) {
    HMM_Vec3 e1 = vec3_sub(p1, p0);
    HMM_Vec3 e2 = vec3_sub(p2, p0);
    HMM_Vec3 face_n = vec3_cross(e1, e2);
    if (vec3_dot(face_n, winding_normal) < 0.0f) {
        HMM_Vec3 tmp = p1; p1 = p2; p2 = tmp;
        HexUV utmp = uv1; uv1 = uv2; uv2 = utmp;
    }

    if (*count + 3 > *cap) {
        *cap = (*cap) * 2;
        *verts = (HexVertex*)realloc(*verts, *cap * sizeof(HexVertex));
    }
    HexVertex* v = &(*verts)[*count];

    for (int i = 0; i < 3; i++) {
        v[i].color[0] = vtx_color.X;
        v[i].color[1] = vtx_color.Y;
        v[i].color[2] = vtx_color.Z;
        v[i].sky_light = sky_light;
        v[i].torch_light = torch_light;
    }

    v[0].pos[0] = p0.X; v[0].pos[1] = p0.Y; v[0].pos[2] = p0.Z;
    v[0].normal[0] = shading_normal.X; v[0].normal[1] = shading_normal.Y; v[0].normal[2] = shading_normal.Z;
    v[0].uv[0] = uv0.u; v[0].uv[1] = uv0.v;

    v[1].pos[0] = p1.X; v[1].pos[1] = p1.Y; v[1].pos[2] = p1.Z;
    v[1].normal[0] = shading_normal.X; v[1].normal[1] = shading_normal.Y; v[1].normal[2] = shading_normal.Z;
    v[1].uv[0] = uv1.u; v[1].uv[1] = uv1.v;

    v[2].pos[0] = p2.X; v[2].pos[1] = p2.Y; v[2].pos[2] = p2.Z;
    v[2].normal[0] = shading_normal.X; v[2].normal[1] = shading_normal.Y; v[2].normal[2] = shading_normal.Z;
    v[2].uv[0] = uv2.u; v[2].uv[1] = uv2.v;

    *count += 3;
}

// ---- Wire push: emit a line segment (2 float3 vertices) ----
static void wire_push_line(float** verts, int* count, int* cap,
                            HMM_Vec3 a, HMM_Vec3 b) {
    if (*count + 6 > *cap) {
        *cap = (*cap == 0) ? 4096 : (*cap) * 2;
        *verts = (float*)realloc(*verts, *cap * sizeof(float));
    }
    float* v = &(*verts)[*count];
    v[0] = a.X; v[1] = a.Y; v[2] = a.Z;
    v[3] = b.X; v[4] = b.Y; v[5] = b.Z;
    *count += 6;
}

static HexUV hex_top_uv(float vx, float vz, float cx, float cz, int atlas_idx) {
    float local_u = (vx - cx) / (2.0f * HEX_RADIUS) + 0.5f;
    float local_v = (vz - cz) / (1.7320508f * HEX_RADIUS) + 0.5f;
    HexUV uv;
    uv.u = ((float)atlas_idx + local_u) / (float)HEX_ATLAS_TILES;
    uv.v = local_v;
    return uv;
}

// ---- Voxel access helpers for mesh job ----

// Get voxel from the job's 3D array (local coords)
static inline uint8_t job_get_voxel(const HexMeshJob* job, int col, int row, int layer) {
    if (layer < 0 || layer >= HEX_CHUNK_LAYERS) return VOXEL_AIR;
    return HEX_VOXEL(job->voxels, col, row, layer);
}

// Get neighbor voxel, handling cross-chunk boundaries.
// Uses actual neighbor chunk voxel data when available, falls back to noise.
static uint8_t job_get_neighbor_voxel(const HexMeshJob* job,
                                       fnl_state* continental, fnl_state* mountain,
                                       fnl_state* warp, fnl_state* detail,
                                       int col, int row, int layer, int dir) {
    int ncol_g, nrow_g;
    hex_neighbor(job->cx * HEX_CHUNK_SIZE + col,
                 job->cz * HEX_CHUNK_SIZE + row, dir, &ncol_g, &nrow_g);

    int local_col = ncol_g - job->cx * HEX_CHUNK_SIZE;
    int local_row = nrow_g - job->cz * HEX_CHUNK_SIZE;

    if (local_col >= 0 && local_col < HEX_CHUNK_SIZE &&
        local_row >= 0 && local_row < HEX_CHUNK_SIZE) {
        // Same chunk: direct lookup
        if (layer < 0 || layer >= HEX_CHUNK_LAYERS) return VOXEL_AIR;
        return HEX_VOXEL(job->voxels, local_col, local_row, layer);
    }

    // Cross-chunk boundary: try actual neighbor data first
    int world_layer = layer + job->base_layer;
    uint8_t nb = job_read_neighbor_voxel(job, ncol_g, nrow_g, world_layer);
    // job_read_neighbor_voxel returns VOXEL_AIR if neighbor data is NULL.
    // If we have actual neighbor data, trust it completely.
    int ncx = (ncol_g < 0) ? -(((-ncol_g - 1) / HEX_CHUNK_SIZE) + 1) : (ncol_g / HEX_CHUNK_SIZE);
    int ncz = (nrow_g < 0) ? -(((-nrow_g - 1) / HEX_CHUNK_SIZE) + 1) : (nrow_g / HEX_CHUNK_SIZE);
    int dx = ncx - job->cx + 1;
    int dz = ncz - job->cz + 1;
    if (dx >= 0 && dx <= 2 && dz >= 0 && dz <= 2 && job->nb_voxels[dx][dz]) {
        return nb;  // Actual neighbor data available — use it
    }

    // Fallback: re-sample terrain from noise (neighbor chunk not loaded)
    float nlx, nlz;
    hex_local_pos(ncol_g, nrow_g, &nlx, &nlz);
    HMM_Vec3 ndir = vec3_normalize(
        vec3_add(job->tangent_origin,
            vec3_add(vec3_scale(job->tangent_east, nlx),
                     vec3_scale(job->tangent_north, nlz))));
    float nh_m = ht_sample_height_m(continental, mountain, warp, detail, ndir);
    float nh_eff = fmaxf(nh_m, TERRAIN_SEA_LEVEL_M);
    int surface_layer = (int)ceilf(nh_eff / HEX_HEIGHT);

    // Check arch at this neighbor too
    float alx, alz;
    ht_arch_project(&job->arch_frame, vec3_normalize(ndir), job->planet_radius, &alx, &alz);
    int ab, at;
    ht_arch_query(alx, alz, job->arch_base, &ab, &at);

    int margin = 2;

    if (ab != INT16_MIN) {
        if (world_layer >= ab - margin && world_layer <= at + margin) {
            if (world_layer >= ab && world_layer <= at)
                return VOXEL_STONE;
            return VOXEL_AIR;
        }
    }

    if (world_layer >= surface_layer - margin && world_layer <= surface_layer + margin)
        return VOXEL_AIR;
    if (world_layer < surface_layer - margin) return VOXEL_STONE;
    return VOXEL_AIR;
}

// ---- BFS sky light propagation (Minecraft-style) ----
// Light level 0 = dark, SKY_MAX = full sunlight.
// Rules:
//   1. Sky column pre-pass: air blocks above the topmost solid in each
//      column get SKY_MAX (direct sunlight streaming down).
//   2. Downward from SKY_MAX through air: stays SKY_MAX (the "sky column
//      rule" — vertical shafts carry full sunlight to bedrock).
//   3. All other propagation: new = current - 1  (taxicab attenuation).
//   4. Only enqueue if the new value improves the neighbor's current value.
// These rules guarantee termination (monotonic decrease, bounded by SKY_MAX).

#define SKY_MAX  32
#define SKY_VOXEL(sky, col, row, layer) \
    ((sky)[(col) * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS + (row) * HEX_CHUNK_LAYERS + (layer)])

// Maximum BFS queue entries. With SKY_MAX=32, worst case per source is
// O(32^2) surface area ≈ ~4000 nodes.  For the whole chunk we cap at the
// full voxel count which is generous but safe.
#define SKY_BFS_CAP (HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS)

// Get the col_max_solid for a hex neighbor, handling cross-chunk boundary via noise.
static int sky_neighbor_max_solid(const HexMeshJob* job,
                                   fnl_state* continental, fnl_state* mountain,
                                   fnl_state* warp, fnl_state* detail,
                                   int gcol, int grow) {
    int local_col = gcol - job->cx * HEX_CHUNK_SIZE;
    int local_row = grow - job->cz * HEX_CHUNK_SIZE;
    if (local_col >= 0 && local_col < HEX_CHUNK_SIZE &&
        local_row >= 0 && local_row < HEX_CHUNK_SIZE) {
        return job->col_max_solid[local_col][local_row];
    }

    // Cross-chunk: sample noise to approximate surface height
    float lx, lz;
    hex_local_pos(gcol, grow, &lx, &lz);
    HMM_Vec3 ndir = vec3_normalize(
        vec3_add(job->tangent_origin,
            vec3_add(vec3_scale(job->tangent_east, lx),
                     vec3_scale(job->tangent_north, lz))));
    float nh_m = ht_sample_height_m(continental, mountain, warp, detail, ndir);
    float nh_eff = fmaxf(nh_m, TERRAIN_SEA_LEVEL_M);
    int surface_layer = (int)ceilf(nh_eff / HEX_HEIGHT);

    // Check arch
    float alx, alz;
    ht_arch_project(&job->arch_frame, ndir, job->planet_radius, &alx, &alz);
    int ab, at;
    ht_arch_query(alx, alz, job->arch_base, &ab, &at);
    if (ab != INT16_MIN && at > surface_layer)
        surface_layer = at;

    int local_max = surface_layer - job->base_layer;
    if (local_max < 0) return -1;
    if (local_max >= HEX_CHUNK_LAYERS) return HEX_CHUNK_LAYERS - 1;
    return local_max;
}

static void compute_sky_light(HexMeshJob* job,
                               fnl_state* continental, fnl_state* mountain,
                               fnl_state* warp, fnl_state* detail) {
    size_t voxel_count = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;
    uint8_t* sky = (uint8_t*)calloc(voxel_count, 1);
    job->sky_map = sky;

    // ---- Pass 1: Sky column fill ----
    // For each column, scan downward from the top of the chunk.
    // Fill all air blocks above (and including) the first unobstructed air
    // with SKY_MAX.  This handles both open-sky columns AND vertical shafts
    // that have been dug out (col_max_solid is recalculated from voxels).
    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            // Scan from the TOP of the chunk downward until we hit a solid (opaque) block
            for (int l = HEX_CHUNK_LAYERS - 1; l >= 0; l--) {
                if (!is_transparent(HEX_VOXEL(job->voxels, col, row, l))) break;
                SKY_VOXEL(sky, col, row, l) = SKY_MAX;
            }
        }
    }

    // Allocate BFS queue (needed for both Pass 1b seeding and Pass 3 BFS)
    int* queue = (int*)malloc(SKY_BFS_CAP * sizeof(int));
    int qhead = 0, qtail = 0;

    // ---- Pass 1b: Cross-chunk sky column seeding ----
    // For boundary cells in our chunk, check cross-chunk hex neighbors.
    // If a neighbor cell is sky-exposed (open sky column above), seed our
    // boundary cell with SKY_MAX-1 so sky light flows across chunk boundaries.
    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            // Only process boundary cells (edges of the chunk grid)
            if (col > 0 && col < HEX_CHUNK_SIZE - 1 &&
                row > 0 && row < HEX_CHUNK_SIZE - 1) continue;

            int gcol = job->cx * HEX_CHUNK_SIZE + col;
            int grow = job->cz * HEX_CHUNK_SIZE + row;

            for (int dir = 0; dir < 6; dir++) {
                int ncol_g, nrow_g;
                hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
                int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
                int nr = nrow_g - job->cz * HEX_CHUNK_SIZE;

                // Only process cross-chunk neighbors
                if (nc >= 0 && nc < HEX_CHUNK_SIZE &&
                    nr >= 0 && nr < HEX_CHUNK_SIZE) continue;

                // Find neighbor chunk data
                int ncx = (ncol_g < 0) ? -(((-ncol_g - 1) / HEX_CHUNK_SIZE) + 1) : (ncol_g / HEX_CHUNK_SIZE);
                int ncz = (nrow_g < 0) ? -(((-nrow_g - 1) / HEX_CHUNK_SIZE) + 1) : (nrow_g / HEX_CHUNK_SIZE);
                int dx = ncx - job->cx + 1;
                int dz = ncz - job->cz + 1;
                if (dx < 0 || dx > 2 || dz < 0 || dz > 2) continue;
                const uint8_t* nb = job->nb_voxels[dx][dz];
                if (!nb) continue;

                int nb_base = job->nb_base_layer[dx][dz];
                int lc = ncol_g - ncx * HEX_CHUNK_SIZE;
                int lr = nrow_g - ncz * HEX_CHUNK_SIZE;
                if (lc < 0 || lc >= HEX_CHUNK_SIZE || lr < 0 || lr >= HEX_CHUNK_SIZE) continue;

                // Find max solid in the neighbor column (scan from top)
                int nb_max_solid = -1;
                for (int sl = HEX_CHUNK_LAYERS - 1; sl >= 0; sl--) {
                    if (!is_transparent(HEX_VOXEL(nb, lc, lr, sl))) {
                        nb_max_solid = sl;
                        break;
                    }
                }

                // Seed our boundary cell from neighbor sky-lit layers (above max solid)
                int sky_start = nb_max_solid + 1;
                if (sky_start < 0) sky_start = 0;
                for (int nl = sky_start; nl < HEX_CHUNK_LAYERS; nl++) {
                    int world_layer = nl + nb_base;
                    int local_layer = world_layer - job->base_layer;
                    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) continue;
                    if (!is_transparent(HEX_VOXEL(job->voxels, col, row, local_layer))) continue;

                    uint8_t seed_val = SKY_MAX - 1;
                    if (seed_val > SKY_VOXEL(sky, col, row, local_layer)) {
                        SKY_VOXEL(sky, col, row, local_layer) = seed_val;
                        if (qtail < SKY_BFS_CAP) {
                            queue[qtail++] = col * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                           + row * HEX_CHUNK_LAYERS + local_layer;
                        }
                    }
                }
            }
        }
    }

    if (qtail > 0) {
        hex_log("[SKY_BFS] chunk(%d,%d) Pass 1b: %d cross-chunk sky seeds\n",
                job->cx, job->cz, qtail);
    }

    // ---- Pass 2: Seed BFS queue ----
    // Any sky-lit air block adjacent to a non-sky-lit air block is a seed.
    // Also seed from sky-lit blocks adjacent to columns with higher terrain
    // (cross-chunk boundary awareness).

    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            int max_s = job->col_max_solid[col][row];
            int min_s = job->col_min_solid[col][row];
            int my_top = (min_s > max_s) ? -1 : max_s;

            int gcol = job->cx * HEX_CHUNK_SIZE + col;
            int grow = job->cz * HEX_CHUNK_SIZE + row;

            // Find tallest neighbor column (for cross-chunk boundary seeding)
            int tallest_neighbor = my_top;
            for (int dir = 0; dir < 6; dir++) {
                int ncol_g, nrow_g;
                hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
                int n_max = sky_neighbor_max_solid(job, continental, mountain, warp, detail,
                                                    ncol_g, nrow_g);
                if (n_max > tallest_neighbor) tallest_neighbor = n_max;
            }

            // Seed sky-lit blocks that overlap a taller neighbor's shadow zone
            if (tallest_neighbor > my_top) {
                int seed_start = my_top + 1;
                if (seed_start < 0) seed_start = 0;
                int seed_end = tallest_neighbor;
                if (seed_end >= HEX_CHUNK_LAYERS) seed_end = HEX_CHUNK_LAYERS - 1;

                for (int l = seed_start; l <= seed_end; l++) {
                    if (SKY_VOXEL(sky, col, row, l) == SKY_MAX && qtail < SKY_BFS_CAP) {
                        queue[qtail++] = col * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                       + row * HEX_CHUNK_LAYERS + l;
                    }
                }
            }

            // Also seed any sky-lit block directly above a solid block
            // (the transition point where light needs to spread laterally)
            if (my_top >= 0 && my_top + 1 < HEX_CHUNK_LAYERS) {
                int l = my_top + 1;
                if (SKY_VOXEL(sky, col, row, l) == SKY_MAX && qtail < SKY_BFS_CAP) {
                    // Check if any lateral neighbor at this layer is non-sky-lit air
                    bool needs_seed = false;
                    for (int dir = 0; dir < 6 && !needs_seed; dir++) {
                        int ncol_g, nrow_g;
                        hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
                        int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
                        int nr = nrow_g - job->cz * HEX_CHUNK_SIZE;
                        if (nc >= 0 && nc < HEX_CHUNK_SIZE && nr >= 0 && nr < HEX_CHUNK_SIZE &&
                            is_transparent(HEX_VOXEL(job->voxels, nc, nr, l)) &&
                            SKY_VOXEL(sky, nc, nr, l) < SKY_MAX) {
                            needs_seed = true;
                        }
                    }
                    if (needs_seed) {
                        queue[qtail++] = col * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                       + row * HEX_CHUNK_LAYERS + l;
                    }
                }
            }
        }
    }

    // ---- Pass 3: BFS flood-fill ----
    while (qhead < qtail) {
        int idx = queue[qhead++];
        int c  = idx / (HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);
        int rem = idx % (HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);
        int r  = rem / HEX_CHUNK_LAYERS;
        int l  = rem % HEX_CHUNK_LAYERS;

        uint8_t current = SKY_VOXEL(sky, c, r, l);
        if (current <= 1) continue;  // Can't propagate further

        // Vertical neighbors (up and down)
        for (int dv = -1; dv <= 1; dv += 2) {
            int nl = l + dv;
            if (nl < 0 || nl >= HEX_CHUNK_LAYERS) continue;
            if (!is_transparent(HEX_VOXEL(job->voxels, c, r, nl))) continue;

            // Sky column rule: downward from SKY_MAX stays SKY_MAX
            uint8_t new_light;
            if (dv == -1 && current == SKY_MAX) {
                new_light = SKY_MAX;
            } else {
                new_light = current - 1;
            }

            if (new_light > SKY_VOXEL(sky, c, r, nl)) {
                SKY_VOXEL(sky, c, r, nl) = new_light;
                if (qtail < SKY_BFS_CAP) {
                    queue[qtail++] = c * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                   + r * HEX_CHUNK_LAYERS + nl;
                }
            }
        }

        // 6 lateral hex neighbors (cross-chunk aware)
        int gcol = job->cx * HEX_CHUNK_SIZE + c;
        int grow = job->cz * HEX_CHUNK_SIZE + r;
        for (int dir = 0; dir < 6; dir++) {
            int ncol_g, nrow_g;
            hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
            int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
            int nr_local = nrow_g - job->cz * HEX_CHUNK_SIZE;

            if (nc >= 0 && nc < HEX_CHUNK_SIZE && nr_local >= 0 && nr_local < HEX_CHUNK_SIZE) {
                // In-chunk neighbor: normal BFS
                if (!is_transparent(HEX_VOXEL(job->voxels, nc, nr_local, l))) continue;
                uint8_t new_light = current - 1;
                if (new_light > SKY_VOXEL(sky, nc, nr_local, l)) {
                    SKY_VOXEL(sky, nc, nr_local, l) = new_light;
                    if (qtail < SKY_BFS_CAP) {
                        queue[qtail++] = nc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                       + nr_local * HEX_CHUNK_LAYERS + l;
                    }
                }
            } else {
                // Cross-chunk neighbor: check if transparent via actual neighbor data,
                // then bounce light back to the originating in-chunk cell's other neighbors.
                // This allows light to "pass through" the boundary.
                int world_layer = l + job->base_layer;
                uint8_t nb_voxel = job_read_neighbor_voxel(job, ncol_g, nrow_g, world_layer);
                if (!is_transparent(nb_voxel)) continue;

                // The cross-chunk cell is transparent — propagate light to
                // in-chunk cells adjacent to the cross-chunk cell.
                // Use -1 cost (same as normal propagation) to avoid visible
                // seams at chunk boundaries.
                uint8_t new_light = current - 1;
                if (new_light < 1) continue;
                for (int rdir = 0; rdir < 6; rdir++) {
                    int rcol_g, rrow_g;
                    hex_neighbor(ncol_g, nrow_g, rdir, &rcol_g, &rrow_g);
                    int rc = rcol_g - job->cx * HEX_CHUNK_SIZE;
                    int rr = rrow_g - job->cz * HEX_CHUNK_SIZE;
                    if (rc < 0 || rc >= HEX_CHUNK_SIZE || rr < 0 || rr >= HEX_CHUNK_SIZE) continue;
                    if (rc == c && rr == r) continue;  // skip self
                    if (!is_transparent(HEX_VOXEL(job->voxels, rc, rr, l))) continue;
                    if (new_light > SKY_VOXEL(sky, rc, rr, l)) {
                        SKY_VOXEL(sky, rc, rr, l) = new_light;
                        if (qtail < SKY_BFS_CAP) {
                            queue[qtail++] = rc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                           + rr * HEX_CHUNK_LAYERS + l;
                        }
                    }
                }
            }
        }
    }

    free(queue);
}

// Read sky light directly from the sky_map for a face.
// Returns 0.0–1.0 (normalized from 0..SKY_MAX).
static float read_sky_light(const HexMeshJob* job, int col, int row, int layer) {
    if (!job->sky_map) return 1.0f;
    if (layer < 0 || layer >= HEX_CHUNK_LAYERS) return 1.0f;
    if (col < 0 || col >= HEX_CHUNK_SIZE || row < 0 || row >= HEX_CHUNK_SIZE) return 1.0f;
    return (float)SKY_VOXEL(job->sky_map, col, row, layer) / (float)SKY_MAX;
}

// ---- Torch BFS lighting ----

#define TORCH_LIGHT_MAX   8
#define TORCH_BFS_CAP     SKY_BFS_CAP
#define TORCH_VOXEL(map, col, row, layer) \
    ((map)[(col) * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS + (row) * HEX_CHUNK_LAYERS + (layer)])

static void compute_torch_light(HexMeshJob* job) {
    size_t voxel_count = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;
    uint8_t* torch = (uint8_t*)calloc(voxel_count, 1);
    job->torch_map = torch;

    int* queue = (int*)malloc(TORCH_BFS_CAP * sizeof(int));
    int qhead = 0, qtail = 0;

    // Seed from all VOXEL_TORCH positions in this chunk
    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            int min_s = job->col_min_solid[col][row];
            int max_s = job->col_max_solid[col][row];
            if (min_s > max_s) continue;
            for (int l = min_s; l <= max_s; l++) {
                if (HEX_VOXEL(job->voxels, col, row, l) == VOXEL_TORCH) {
                    TORCH_VOXEL(torch, col, row, l) = TORCH_LIGHT_MAX;
                    if (qtail < TORCH_BFS_CAP) {
                        queue[qtail++] = col * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                       + row * HEX_CHUNK_LAYERS + l;
                    }
                }
            }
        }
    }

    // Seed from cross-chunk torches: scan neighbor chunk border strips for VOXEL_TORCH.
    // Only scan the TORCH_LIGHT_MAX-wide border region facing our chunk.
    int cross_torch_seeds = 0;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dz = -1; dz <= 1; dz++) {
            if (dx == 0 && dz == 0) continue;
            const uint8_t* nb = job->nb_voxels[dx + 1][dz + 1];
            if (!nb) continue;
            int nb_base = job->nb_base_layer[dx + 1][dz + 1];
            int ncx = job->cx + dx;
            int ncz = job->cz + dz;
            // Determine scan range: only border region within torch range
            int c_lo = 0, c_hi = HEX_CHUNK_SIZE;
            int r_lo = 0, r_hi = HEX_CHUNK_SIZE;
            if (dx == -1) c_lo = HEX_CHUNK_SIZE - TORCH_LIGHT_MAX;  // right edge of left neighbor
            if (dx ==  1) c_hi = TORCH_LIGHT_MAX;                    // left edge of right neighbor
            if (dz == -1) r_lo = HEX_CHUNK_SIZE - TORCH_LIGHT_MAX;
            if (dz ==  1) r_hi = TORCH_LIGHT_MAX;
            if (c_lo < 0) c_lo = 0;
            if (r_lo < 0) r_lo = 0;

            for (int nc = c_lo; nc < c_hi; nc++) {
                for (int nr = r_lo; nr < r_hi; nr++) {
                    for (int nl = 0; nl < HEX_CHUNK_LAYERS; nl++) {
                        if (HEX_VOXEL(nb, nc, nr, nl) != VOXEL_TORCH) continue;
                        int torch_gcol = ncx * HEX_CHUNK_SIZE + nc;
                        int torch_grow = ncz * HEX_CHUNK_SIZE + nr;
                        int torch_wl = nl + nb_base;
                        // Seed in-chunk cells adjacent to this torch
                        for (int tdir = 0; tdir < 6; tdir++) {
                            int adj_g, adj_r;
                            hex_neighbor(torch_gcol, torch_grow, tdir, &adj_g, &adj_r);
                            int lc = adj_g - job->cx * HEX_CHUNK_SIZE;
                            int lr = adj_r - job->cz * HEX_CHUNK_SIZE;
                            if (lc < 0 || lc >= HEX_CHUNK_SIZE || lr < 0 || lr >= HEX_CHUNK_SIZE) continue;
                            int ll = torch_wl - job->base_layer;
                            if (ll < 0 || ll >= HEX_CHUNK_LAYERS) continue;
                            if (!is_transparent(HEX_VOXEL(job->voxels, lc, lr, ll))) continue;
                            uint8_t seed_val = TORCH_LIGHT_MAX - 1;
                            if (seed_val > TORCH_VOXEL(torch, lc, lr, ll)) {
                                TORCH_VOXEL(torch, lc, lr, ll) = seed_val;
                                cross_torch_seeds++;
                                if (qtail < TORCH_BFS_CAP) {
                                    queue[qtail++] = lc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                                   + lr * HEX_CHUNK_LAYERS + ll;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (cross_torch_seeds > 0 || qtail > 0) {
        hex_log("[TORCH_BFS] chunk(%d,%d) local_torches=%d cross_chunk_seeds=%d\n",
                job->cx, job->cz, qtail - cross_torch_seeds, cross_torch_seeds);
    }

    // BFS flood-fill (same structure as sky light Pass 3, but no sky-column rule)
    while (qhead < qtail) {
        int idx = queue[qhead++];
        int c   = idx / (HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);
        int rem = idx % (HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);
        int r   = rem / HEX_CHUNK_LAYERS;
        int l   = rem % HEX_CHUNK_LAYERS;

        uint8_t current = TORCH_VOXEL(torch, c, r, l);
        if (current <= 1) continue;

        // Vertical neighbors (up and down) — always decay by 1
        for (int dv = -1; dv <= 1; dv += 2) {
            int nl = l + dv;
            if (nl < 0 || nl >= HEX_CHUNK_LAYERS) continue;
            uint8_t nv = HEX_VOXEL(job->voxels, c, r, nl);
            if (!is_transparent(nv)) continue;

            uint8_t new_light = current - 1;
            if (new_light > TORCH_VOXEL(torch, c, r, nl)) {
                TORCH_VOXEL(torch, c, r, nl) = new_light;
                if (qtail < TORCH_BFS_CAP) {
                    queue[qtail++] = c * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                   + r * HEX_CHUNK_LAYERS + nl;
                }
            }
        }

        // 6 lateral hex neighbors (cross-chunk aware)
        int gcol = job->cx * HEX_CHUNK_SIZE + c;
        int grow = job->cz * HEX_CHUNK_SIZE + r;
        for (int dir = 0; dir < 6; dir++) {
            int ncol_g, nrow_g;
            hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
            int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
            int nr_local = nrow_g - job->cz * HEX_CHUNK_SIZE;

            if (nc >= 0 && nc < HEX_CHUNK_SIZE && nr_local >= 0 && nr_local < HEX_CHUNK_SIZE) {
                // In-chunk neighbor
                uint8_t nv = HEX_VOXEL(job->voxels, nc, nr_local, l);
                if (!is_transparent(nv)) continue;
                uint8_t new_light = current - 1;
                if (new_light > TORCH_VOXEL(torch, nc, nr_local, l)) {
                    TORCH_VOXEL(torch, nc, nr_local, l) = new_light;
                    if (qtail < TORCH_BFS_CAP) {
                        queue[qtail++] = nc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                       + nr_local * HEX_CHUNK_LAYERS + l;
                    }
                }
            } else {
                // Cross-chunk: propagate through transparent neighbor at -1 cost
                int world_layer = l + job->base_layer;
                uint8_t nb_voxel = job_read_neighbor_voxel(job, ncol_g, nrow_g, world_layer);
                if (!is_transparent(nb_voxel)) continue;
                uint8_t new_light = current - 1;
                if (new_light < 1) continue;
                for (int rdir = 0; rdir < 6; rdir++) {
                    int rcol_g, rrow_g;
                    hex_neighbor(ncol_g, nrow_g, rdir, &rcol_g, &rrow_g);
                    int rc = rcol_g - job->cx * HEX_CHUNK_SIZE;
                    int rr = rrow_g - job->cz * HEX_CHUNK_SIZE;
                    if (rc < 0 || rc >= HEX_CHUNK_SIZE || rr < 0 || rr >= HEX_CHUNK_SIZE) continue;
                    if (rc == c && rr == r) continue;
                    if (!is_transparent(HEX_VOXEL(job->voxels, rc, rr, l))) continue;
                    if (new_light > TORCH_VOXEL(torch, rc, rr, l)) {
                        TORCH_VOXEL(torch, rc, rr, l) = new_light;
                        if (qtail < TORCH_BFS_CAP) {
                            queue[qtail++] = rc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                           + rr * HEX_CHUNK_LAYERS + l;
                        }
                    }
                }
            }
        }
    }

    // Seed from neighbor boundary torch values (propagate light across chunk borders).
    // On the second pass, neighbor chunks have computed their BFS and stored border values.
    // We read the facing border and seed our edge cells with value - 1.
    // IMPORTANT: border values use the SOURCE chunk's local layer indexing.
    // We must convert via base_layer offset to our local layer space.
    {
        int border_seeds = 0;
        // Cardinal neighbors only: dx or dz is 0 (not both nonzero)
        static const int cardinal[4][2] = {{0,-1},{0,1},{-1,0},{1,0}}; // S,N,W,E
        for (int ci = 0; ci < 4; ci++) {
            int dx = cardinal[ci][0], dz = cardinal[ci][1];
            const uint8_t* nb_border = job->nb_border_torch[dx+1][dz+1];
            if (!nb_border) continue;

            // Layer offset: convert source local layer → our local layer
            int nb_base = job->nb_base_layer[dx+1][dz+1];
            int layer_off = nb_base - job->base_layer;

            for (int cell = 0; cell < HEX_CHUNK_SIZE; cell++) {
                for (int src_l = 0; src_l < HEX_CHUNK_LAYERS; src_l++) {
                    uint8_t nb_val = nb_border[cell * HEX_CHUNK_LAYERS + src_l];
                    if (nb_val <= 1) continue;
                    uint8_t seed_val = nb_val - 1;

                    // Convert source layer to our local layer
                    int l = src_l + layer_off;
                    if (l < 0 || l >= HEX_CHUNK_LAYERS) continue;

                    // Map to our edge cell
                    int c, r;
                    if (dx == -1) { c = 0; r = cell; }                // west neighbor → our col 0
                    else if (dx == 1) { c = HEX_CHUNK_SIZE - 1; r = cell; } // east → our col 31
                    else if (dz == -1) { c = cell; r = 0; }           // south → our row 0
                    else { c = cell; r = HEX_CHUNK_SIZE - 1; }        // north → our row 31

                    if (!is_transparent(HEX_VOXEL(job->voxels, c, r, l))) continue;
                    if (seed_val > TORCH_VOXEL(torch, c, r, l)) {
                        TORCH_VOXEL(torch, c, r, l) = seed_val;
                        border_seeds++;
                        if (qtail < TORCH_BFS_CAP) {
                            queue[qtail++] = c * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                           + r * HEX_CHUNK_LAYERS + l;
                        }
                    }
                }
            }
        }
        if (border_seeds > 0) {
            hex_log("[TORCH_BORDER] chunk(%d,%d) seeded %d cells from neighbor borders\n",
                    job->cx, job->cz, border_seeds);
        } else {
            // Debug: check if any neighbor had border data at all
            for (int ci = 0; ci < 4; ci++) {
                int dx = cardinal[ci][0], dz = cardinal[ci][1];
                const uint8_t* nb_border = job->nb_border_torch[dx+1][dz+1];
                if (nb_border) {
                    int nb_base = job->nb_base_layer[dx+1][dz+1];
                    int layer_off = nb_base - job->base_layer;
                    int nb_nonzero = 0;
                    for (int i = 0; i < HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS; i++)
                        if (nb_border[i] > 0) nb_nonzero++;
                    hex_log("[TORCH_BORDER_DBG] chunk(%d,%d) nb(%d,%d) has %d nonzero border cells, layer_off=%d (nb_base=%d, our_base=%d)\n",
                            job->cx, job->cz, dx, dz, nb_nonzero, layer_off, nb_base, job->base_layer);
                }
            }
        }
    }

    // Re-run BFS from border seeds (they were added to the queue above)
    while (qhead < qtail) {
        int idx = queue[qhead++];
        int c   = idx / (HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);
        int rem = idx % (HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);
        int r   = rem / HEX_CHUNK_LAYERS;
        int l   = rem % HEX_CHUNK_LAYERS;

        uint8_t current = TORCH_VOXEL(torch, c, r, l);
        if (current <= 1) continue;

        for (int dv = -1; dv <= 1; dv += 2) {
            int nl = l + dv;
            if (nl < 0 || nl >= HEX_CHUNK_LAYERS) continue;
            if (!is_transparent(HEX_VOXEL(job->voxels, c, r, nl))) continue;
            uint8_t new_light = current - 1;
            if (new_light > TORCH_VOXEL(torch, c, r, nl)) {
                TORCH_VOXEL(torch, c, r, nl) = new_light;
                if (qtail < TORCH_BFS_CAP)
                    queue[qtail++] = c * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                   + r * HEX_CHUNK_LAYERS + nl;
            }
        }
        int gcol = job->cx * HEX_CHUNK_SIZE + c;
        int grow = job->cz * HEX_CHUNK_SIZE + r;
        for (int dir = 0; dir < 6; dir++) {
            int ncol_g, nrow_g;
            hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
            int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
            int nr_local = nrow_g - job->cz * HEX_CHUNK_SIZE;
            if (nc < 0 || nc >= HEX_CHUNK_SIZE || nr_local < 0 || nr_local >= HEX_CHUNK_SIZE) continue;
            if (!is_transparent(HEX_VOXEL(job->voxels, nc, nr_local, l))) continue;
            uint8_t new_light = current - 1;
            if (new_light > TORCH_VOXEL(torch, nc, nr_local, l)) {
                TORCH_VOXEL(torch, nc, nr_local, l) = new_light;
                if (qtail < TORCH_BFS_CAP)
                    queue[qtail++] = nc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                   + nr_local * HEX_CHUNK_LAYERS + l;
            }
        }
    }

    // Extract boundary torch values for cross-chunk propagation
    {
        size_t border_size = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;
        for (int b = 0; b < 4; b++) {
            bool has_light = false;
            // Check if any border cell has torch light
            for (int cell = 0; cell < HEX_CHUNK_SIZE && !has_light; cell++) {
                for (int l = 0; l < HEX_CHUNK_LAYERS && !has_light; l++) {
                    int c, r;
                    if (b == 0) { c = 0; r = cell; }
                    else if (b == 1) { c = HEX_CHUNK_SIZE - 1; r = cell; }
                    else if (b == 2) { c = cell; r = 0; }
                    else { c = cell; r = HEX_CHUNK_SIZE - 1; }
                    if (TORCH_VOXEL(torch, c, r, l) > 0) has_light = true;
                }
            }
            if (has_light) {
                job->out_border_torch[b] = (uint8_t*)malloc(border_size);
                for (int cell = 0; cell < HEX_CHUNK_SIZE; cell++) {
                    int c, r;
                    if (b == 0) { c = 0; r = cell; }
                    else if (b == 1) { c = HEX_CHUNK_SIZE - 1; r = cell; }
                    else if (b == 2) { c = cell; r = 0; }
                    else { c = cell; r = HEX_CHUNK_SIZE - 1; }
                    memcpy(&job->out_border_torch[b][cell * HEX_CHUNK_LAYERS],
                           &torch[c * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS + r * HEX_CHUNK_LAYERS],
                           HEX_CHUNK_LAYERS);
                }
            }
        }
    }

    free(queue);
}

static float read_torch_light(const HexMeshJob* job, int col, int row, int layer) {
    if (!job->torch_map) return 0.0f;
    if (layer < 0 || layer >= HEX_CHUNK_LAYERS) return 0.0f;
    if (col < 0 || col >= HEX_CHUNK_SIZE || row < 0 || row >= HEX_CHUNK_SIZE) return 0.0f;
    return (float)TORCH_VOXEL(job->torch_map, col, row, layer) / (float)TORCH_LIGHT_MAX;
}

// Estimate sky light for a cross-chunk air block by scanning neighbor voxel column.
// gcol/grow are global hex coords, world_layer is the absolute layer.
static float estimate_nb_sky(const HexMeshJob* job, int gcol, int grow, int world_layer) {
    int ncx = (gcol < 0) ? -(((-gcol - 1) / HEX_CHUNK_SIZE) + 1) : (gcol / HEX_CHUNK_SIZE);
    int ncz = (grow < 0) ? -(((-grow - 1) / HEX_CHUNK_SIZE) + 1) : (grow / HEX_CHUNK_SIZE);
    int dx = ncx - job->cx + 1;
    int dz = ncz - job->cz + 1;
    if (dx < 0 || dx > 2 || dz < 0 || dz > 2) return 1.0f;
    const uint8_t* nb = job->nb_voxels[dx][dz];
    if (!nb) return 1.0f;

    int nb_base = job->nb_base_layer[dx][dz];
    int lc = gcol - ncx * HEX_CHUNK_SIZE;
    int lr = grow - ncz * HEX_CHUNK_SIZE;
    int ll = world_layer - nb_base;
    if (lc < 0 || lc >= HEX_CHUNK_SIZE || lr < 0 || lr >= HEX_CHUNK_SIZE) return 1.0f;
    if (ll < 0 || ll >= HEX_CHUNK_LAYERS) return 1.0f;

    // Scan upward from this cell: if all air above → sky column → full light
    for (int sl = ll; sl < HEX_CHUNK_LAYERS; sl++) {
        if (!is_transparent(HEX_VOXEL(nb, lc, lr, sl))) {
            // Covered by solid: light can only arrive laterally.
            // Use conservative dim value — accurate BFS is only in our own chunk.
            return 0.15f;
        }
    }
    return 1.0f;  // open sky column
}

// Estimate torch light for a cross-chunk air block by scanning nearby cells for torches.
static float estimate_nb_torch(const HexMeshJob* job, int gcol, int grow, int world_layer) {
    int ncx = (gcol < 0) ? -(((-gcol - 1) / HEX_CHUNK_SIZE) + 1) : (gcol / HEX_CHUNK_SIZE);
    int ncz = (grow < 0) ? -(((-grow - 1) / HEX_CHUNK_SIZE) + 1) : (grow / HEX_CHUNK_SIZE);
    int dx = ncx - job->cx + 1;
    int dz = ncz - job->cz + 1;
    if (dx < 0 || dx > 2 || dz < 0 || dz > 2) return 0.0f;
    const uint8_t* nb = job->nb_voxels[dx][dz];
    if (!nb) return 0.0f;

    int nb_base = job->nb_base_layer[dx][dz];
    int lc = gcol - ncx * HEX_CHUNK_SIZE;
    int lr = grow - ncz * HEX_CHUNK_SIZE;
    int ll = world_layer - nb_base;
    if (lc < 0 || lc >= HEX_CHUNK_SIZE || lr < 0 || lr >= HEX_CHUNK_SIZE) return 0.0f;
    if (ll < 0 || ll >= HEX_CHUNK_LAYERS) return 0.0f;

    // Check if the cell itself is a torch
    if (HEX_VOXEL(nb, lc, lr, ll) == VOXEL_TORCH)
        return 1.0f;

    // Check immediate neighbors (6 lateral + 2 vertical) for torches
    float best = 0.0f;
    for (int dv = -1; dv <= 1; dv++) {
        int sl = ll + dv;
        if (sl < 0 || sl >= HEX_CHUNK_LAYERS) continue;
        if (HEX_VOXEL(nb, lc, lr, sl) == VOXEL_TORCH) {
            float v = (float)(TORCH_LIGHT_MAX - 1) / (float)TORCH_LIGHT_MAX;
            if (v > best) best = v;
        }
    }
    // Check lateral neighbors within the same neighbor chunk
    int ng_gcol = gcol, ng_grow = grow;
    for (int dir = 0; dir < 6; dir++) {
        int hcol, hrow;
        hex_neighbor(ng_gcol, ng_grow, dir, &hcol, &hrow);
        int hlc = hcol - ncx * HEX_CHUNK_SIZE;
        int hlr = hrow - ncz * HEX_CHUNK_SIZE;
        if (hlc < 0 || hlc >= HEX_CHUNK_SIZE || hlr < 0 || hlr >= HEX_CHUNK_SIZE) continue;
        if (ll >= 0 && ll < HEX_CHUNK_LAYERS && HEX_VOXEL(nb, hlc, hlr, ll) == VOXEL_TORCH) {
            float v = (float)(TORCH_LIGHT_MAX - 1) / (float)TORCH_LIGHT_MAX;
            if (v > best) best = v;
        }
    }
    return best;
}

// ---- Timing helper (thread-safe) ----
static double timer_ms(void) {
#ifdef _WIN32
    LARGE_INTEGER freq, now;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&now);
    return (double)now.QuadPart * 1000.0 / (double)freq.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1000000.0;
#endif
}

// ---- Core mesh generation for a chunk ----

static void generate_chunk_mesh(HexMeshJob* job) {
    double t_start = timer_ms();

    fnl_state continental = ht_create_continental_noise(job->seed);
    fnl_state mountain = ht_create_mountain_noise(job->seed);
    fnl_state warp = ht_create_warp_noise(job->seed);
    fnl_state detail = ht_create_detail_noise(job->seed);

    // Precompute arch data once per job (used in voxel fill + cross-chunk queries)
    job->arch_frame = ht_arch_make_frame();
    job->arch_base = ht_arch_compute_base(&job->arch_frame, job->seed);

    // Skip voxel generation if job already has populated voxels (remesh after break/place)
    if (!job->has_voxels) {
        // Phase 1: Compute base_layer by sampling terrain at chunk corners + center
        {
            int min_h = INT32_MAX;
            int samples[][2] = {{0,0}, {31,0}, {0,31}, {31,31}, {16,16}};
            for (int s = 0; s < 5; s++) {
                int gcol = job->cx * HEX_CHUNK_SIZE + samples[s][0];
                int grow = job->cz * HEX_CHUNK_SIZE + samples[s][1];
                float lx, lz;
                hex_local_pos(gcol, grow, &lx, &lz);
                HMM_Vec3 dir = vec3_normalize(
                    vec3_add(job->tangent_origin,
                        vec3_add(vec3_scale(job->tangent_east, lx),
                                 vec3_scale(job->tangent_north, lz))));
                float h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, dir);
                int h_layers = (int)floorf(fmaxf(h_m, TERRAIN_SEA_LEVEL_M) / HEX_HEIGHT);
                if (h_layers < min_h) min_h = h_layers;
            }
            job->base_layer = min_h - HEX_BASE_PADDING;
        }

        // Phase 2: Fill 3D voxel array
        memset(job->voxels, VOXEL_AIR,
               (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS);

        for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
            for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
                int gcol = job->cx * HEX_CHUNK_SIZE + col;
                int grow = job->cz * HEX_CHUNK_SIZE + row;

                float lx, lz;
                hex_local_pos(gcol, grow, &lx, &lz);

                HMM_Vec3 world_dir = vec3_add(job->tangent_origin,
                    vec3_add(vec3_scale(job->tangent_east, lx),
                             vec3_scale(job->tangent_north, lz)));
                HMM_Vec3 unit = vec3_normalize(world_dir);

                float h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, unit);
                float effective_h = fmaxf(h_m, TERRAIN_SEA_LEVEL_M);
                int surface_layer = (int)ceilf(effective_h / HEX_HEIGHT);
                VoxelType surface_type = ht_voxel_type(h_m, gcol, grow, job->seed); // biome from raw h_m

                int min_solid = HEX_CHUNK_LAYERS;
                int max_solid = -1;

                for (int l = 0; l < HEX_CHUNK_LAYERS; l++) {
                    int world_layer = l + job->base_layer;
                    if (world_layer <= surface_layer) {
                        int depth = surface_layer - world_layer;
                        VoxelType vtype = ht_voxel_type_at_depth(surface_type, depth);
                        HEX_VOXEL(job->voxels, col, row, l) = (uint8_t)vtype;
                        if (l < min_solid) min_solid = l;
                        if (l > max_solid) max_solid = l;
                    }
                }

                job->col_min_solid[col][row] = (int16_t)min_solid;
                job->col_max_solid[col][row] = (int16_t)max_solid;
            }
        }

        // Phase 2b: Fill arch voxels (stone arch at spawn)
        {
            for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
                for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
                    int gcol = job->cx * HEX_CHUNK_SIZE + col;
                    int grow = job->cz * HEX_CHUNK_SIZE + row;

                    float lx, lz;
                    hex_local_pos(gcol, grow, &lx, &lz);

                    HMM_Vec3 world_dir = vec3_add(job->tangent_origin,
                        vec3_add(vec3_scale(job->tangent_east, lx),
                                 vec3_scale(job->tangent_north, lz)));
                    HMM_Vec3 unit = vec3_normalize(world_dir);

                    float alx, alz;
                    ht_arch_project(&job->arch_frame, unit, job->planet_radius, &alx, &alz);

                    int ab, at;
                    ht_arch_query(alx, alz, job->arch_base, &ab, &at);
                    if (ab == INT16_MIN) continue;

                    int min_solid = job->col_min_solid[col][row];
                    int max_solid = job->col_max_solid[col][row];

                    for (int wl = ab; wl <= at; wl++) {
                        int l = wl - job->base_layer;
                        if (l < 0 || l >= HEX_CHUNK_LAYERS) continue;
                        HEX_VOXEL(job->voxels, col, row, l) = (uint8_t)VOXEL_STONE;
                        if (l < min_solid) min_solid = l;
                        if (l > max_solid) max_solid = l;
                    }

                    job->col_min_solid[col][row] = (int16_t)min_solid;
                    job->col_max_solid[col][row] = (int16_t)max_solid;
                }
            }
        }

        // Phase 2b2: Apply persisted voxel edits from disk cache
        if (job->edit_cache) {
            int8_t qface;
            int16_t qu_min, qu_max, qv_min, qv_max;
            compute_chunk_sector_bounds(job, &qface, &qu_min, &qu_max, &qv_min, &qv_max);

            VoxelEdit chunk_edits[512];
            int num_edits = edit_cache_query_area(job->edit_cache,
                qface, qu_min, qu_max, qv_min, qv_max,
                chunk_edits, 512);

            for (int e = 0; e < num_edits; e++) {
                int egcol, egrow;
                unit_to_gcol_grow(job, chunk_edits[e].ux, chunk_edits[e].uy,
                                   chunk_edits[e].uz, &egcol, &egrow);

                int local_col = egcol - job->cx * HEX_CHUNK_SIZE;
                int local_row = egrow - job->cz * HEX_CHUNK_SIZE;
                if (local_col < 0 || local_col >= HEX_CHUNK_SIZE ||
                    local_row < 0 || local_row >= HEX_CHUNK_SIZE) continue;

                int local_layer = chunk_edits[e].layer - job->base_layer;
                if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) continue;

                HEX_VOXEL(job->voxels, local_col, local_row, local_layer) =
                    chunk_edits[e].voxel_type;

                // Update column cache
                if (chunk_edits[e].voxel_type != VOXEL_AIR) {
                    if (local_layer < job->col_min_solid[local_col][local_row])
                        job->col_min_solid[local_col][local_row] = (int16_t)local_layer;
                    if (local_layer > job->col_max_solid[local_col][local_row])
                        job->col_max_solid[local_col][local_row] = (int16_t)local_layer;
                }
            }
        }
    }

    double t_voxels = timer_ms();

    // Phase 2c: Fill bedrock floor at bottom of chunk (unbreakable)
    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            for (int l = 0; l < HEX_BEDROCK_LAYERS; l++) {
                HEX_VOXEL(job->voxels, col, row, l) = (uint8_t)VOXEL_BEDROCK;
            }
            // Update column cache to include bedrock
            if (job->col_min_solid[col][row] > 0)
                job->col_min_solid[col][row] = 0;
            if (job->col_max_solid[col][row] < HEX_BEDROCK_LAYERS - 1)
                job->col_max_solid[col][row] = HEX_BEDROCK_LAYERS - 1;
        }
    }

    // Phase 2d: BFS sky light propagation
    double t_bedrock = timer_ms();
    compute_sky_light(job, &continental, &mountain, &warp, &detail);
    // Phase 2e: BFS torch light propagation
    compute_torch_light(job);

    // Log torch light grid around any torch in this chunk (single atomic write)
    if (job->torch_map) {
        for (int tc = 0; tc < HEX_CHUNK_SIZE; tc++) {
            for (int tr = 0; tr < HEX_CHUNK_SIZE; tr++) {
                int min_s = job->col_min_solid[tc][tr];
                int max_s = job->col_max_solid[tc][tr];
                for (int tl = min_s; tl <= max_s; tl++) {
                    if (HEX_VOXEL(job->voxels, tc, tr, tl) == VOXEL_TORCH) {
                        int gcol_t = job->cx * HEX_CHUNK_SIZE + tc;
                        int grow_t = job->cz * HEX_CHUNK_SIZE + tr;
                        // Build entire grid into buffer, then write atomically
                        char buf[4096];
                        int pos = 0;
                        pos += snprintf(buf + pos, sizeof(buf) - pos,
                            "\n[TORCH_GRID] chunk(%d,%d) torch@local(%d,%d) global(%d,%d) layer=%d nb=%d\n",
                            job->cx, job->cz, tc, tr, gcol_t, grow_t, tl + job->base_layer, job->nb_count);
                        pos += snprintf(buf + pos, sizeof(buf) - pos, "     ");
                        for (int c = tc - 8; c <= tc + 8; c++)
                            pos += snprintf(buf + pos, sizeof(buf) - pos, " %3d", c);
                        pos += snprintf(buf + pos, sizeof(buf) - pos, "\n");
                        for (int r = tr - 8; r <= tr + 8; r++) {
                            pos += snprintf(buf + pos, sizeof(buf) - pos, "r%3d:", r);
                            for (int c = tc - 8; c <= tc + 8; c++) {
                                if (c >= 0 && c < HEX_CHUNK_SIZE && r >= 0 && r < HEX_CHUNK_SIZE) {
                                    int val = TORCH_VOXEL(job->torch_map, c, r, tl);
                                    pos += snprintf(buf + pos, sizeof(buf) - pos, " %3d", val);
                                } else {
                                    pos += snprintf(buf + pos, sizeof(buf) - pos, "   .");
                                }
                                if (pos >= (int)sizeof(buf) - 20) break;
                            }
                            pos += snprintf(buf + pos, sizeof(buf) - pos, "\n");
                            if (pos >= (int)sizeof(buf) - 100) break;
                        }
                        hex_log("%s", buf);
                        goto done_torch_log;
                    }
                }
            }
        }
        done_torch_log:;

        // Also log grid for chunks that received border seeds but have NO local torch
        // (shows cross-chunk propagation into neighbor chunks)
        bool has_local_torch = false;
        for (int tc2 = 0; tc2 < HEX_CHUNK_SIZE && !has_local_torch; tc2++)
            for (int tr2 = 0; tr2 < HEX_CHUNK_SIZE && !has_local_torch; tr2++) {
                int min_s = job->col_min_solid[tc2][tr2];
                int max_s = job->col_max_solid[tc2][tr2];
                for (int tl2 = min_s; tl2 <= max_s; tl2++)
                    if (HEX_VOXEL(job->voxels, tc2, tr2, tl2) == VOXEL_TORCH) { has_local_torch = true; break; }
            }
        if (!has_local_torch) {
            // Find the brightest torch-lit cell to center the grid on
            int best_c = -1, best_r = -1, best_l = -1, best_v = 0;
            for (int c = 0; c < HEX_CHUNK_SIZE; c++)
                for (int r = 0; r < HEX_CHUNK_SIZE; r++) {
                    int min_s = job->col_min_solid[c][r];
                    int max_s = job->col_max_solid[c][r];
                    for (int l = min_s; l <= max_s; l++) {
                        int v = TORCH_VOXEL(job->torch_map, c, r, l);
                        if (v > best_v) { best_v = v; best_c = c; best_r = r; best_l = l; }
                    }
                }
            if (best_v > 0) {
                char buf[4096];
                int pos = 0;
                pos += snprintf(buf + pos, sizeof(buf) - pos,
                    "\n[TORCH_RECV] chunk(%d,%d) brightest=%d @local(%d,%d) layer=%d nb=%d\n",
                    job->cx, job->cz, best_v, best_c, best_r, best_l + job->base_layer, job->nb_count);
                pos += snprintf(buf + pos, sizeof(buf) - pos, "     ");
                for (int c = best_c - 8; c <= best_c + 8; c++)
                    pos += snprintf(buf + pos, sizeof(buf) - pos, " %3d", c);
                pos += snprintf(buf + pos, sizeof(buf) - pos, "\n");
                for (int r = best_r - 8; r <= best_r + 8; r++) {
                    pos += snprintf(buf + pos, sizeof(buf) - pos, "r%3d:", r);
                    for (int c = best_c - 8; c <= best_c + 8; c++) {
                        if (c >= 0 && c < HEX_CHUNK_SIZE && r >= 0 && r < HEX_CHUNK_SIZE) {
                            pos += snprintf(buf + pos, sizeof(buf) - pos, " %3d",
                                TORCH_VOXEL(job->torch_map, c, r, best_l));
                        } else {
                            pos += snprintf(buf + pos, sizeof(buf) - pos, "   .");
                        }
                        if (pos >= (int)sizeof(buf) - 20) break;
                    }
                    pos += snprintf(buf + pos, sizeof(buf) - pos, "\n");
                    if (pos >= (int)sizeof(buf) - 100) break;
                }
                hex_log("%s", buf);
            }
        }
    }

    double t_skylight = timer_ms();

    // Phase 3: Generate mesh geometry
    int cap = HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * 48; // generous initial
    job->vertices = (HexVertex*)malloc(cap * sizeof(HexVertex));
    job->vertex_count = 0;

    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            int gcol = job->cx * HEX_CHUNK_SIZE + col;
            int grow = job->cz * HEX_CHUNK_SIZE + row;

            float lx, lz;
            hex_local_pos(gcol, grow, &lx, &lz);

            HMM_Vec3 center_dir = vec3_normalize(
                vec3_add(job->tangent_origin,
                    vec3_add(vec3_scale(job->tangent_east, lx),
                             vec3_scale(job->tangent_north, lz))));
            HMM_Vec3 local_up = center_dir;
            HMM_Vec3 local_down = vec3_scale(center_dir, -1.0f);

            // Compute per-column terrain color (no brightness bake — shader handles lighting)
            float col_h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, center_dir);
            HMM_Vec3 col_color = ht_terrain_color(col_h_m);

            // Compute terrain slope normal from noise at LOD-triangle-scale offsets (~10m).
            // Used as shading normal so the fragment shader's NdotL + AO match the LOD mesh.
            HMM_Vec3 slope_normal = local_up;
            {
                const float d = 10.0f;  // ~LOD depth-12 vertex spacing
                HMM_Vec3 dir_e = vec3_normalize(vec3_add(job->tangent_origin,
                    vec3_add(vec3_scale(job->tangent_east, lx + d),
                             vec3_scale(job->tangent_north, lz))));
                HMM_Vec3 dir_n = vec3_normalize(vec3_add(job->tangent_origin,
                    vec3_add(vec3_scale(job->tangent_east, lx),
                             vec3_scale(job->tangent_north, lz + d))));
                float h_e = ht_sample_height_m(&continental, &mountain, &warp, &detail, dir_e);
                float h_n = ht_sample_height_m(&continental, &mountain, &warp, &detail, dir_n);

                float grad_e = (h_e - col_h_m) / d;
                float grad_n = (h_n - col_h_m) / d;

                slope_normal = vec3_normalize(vec3_sub(local_up,
                    vec3_add(vec3_scale(job->tangent_east, grad_e),
                             vec3_scale(job->tangent_north, grad_n))));
                if (vec3_dot(slope_normal, local_up) < 0.1f) slope_normal = local_up;
            }

            if (job->is_transition) {
                // Transition zone: smooth per-vertex surface, top face only (no 3D)
                int min_s = job->col_min_solid[col][row];
                int max_s = job->col_max_solid[col][row];
                if (min_s > max_s) continue;

                VoxelType surface_type = (VoxelType)HEX_VOXEL(job->voxels, col, row, max_s);
                int top_atlas = voxel_top_atlas(surface_type);

                HMM_Vec3 hex_verts_top[6];
                HexUV top_uvs[6];
                for (int i = 0; i < 6; i++) {
                    float vx = lx + HEX_RADIUS * HEX_COS[i];
                    float vz = lz + HEX_RADIUS * HEX_SIN[i];
                    HMM_Vec3 vdir = vec3_normalize(
                        vec3_add(job->tangent_origin,
                            vec3_add(vec3_scale(job->tangent_east, vx),
                                     vec3_scale(job->tangent_north, vz))));
                    float v_h = ht_sample_height_m(&continental, &mountain, &warp, &detail, vdir);
                    float v_eff = fmaxf(v_h, TERRAIN_SEA_LEVEL_M);
                    float v_r = job->planet_radius + v_eff + HEX_SURFACE_BIAS;
                    hex_verts_top[i] = vec3_scale(vdir, v_r);
                    top_uvs[i] = hex_top_uv(vx, vz, lx, lz, top_atlas);
                }

                for (int i = 0; i < 4; i++) {
                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        hex_verts_top[0], hex_verts_top[i + 1], hex_verts_top[i + 2],
                        local_up, slope_normal, top_uvs[0], top_uvs[i + 1], top_uvs[i + 2],
                        col_color, 1.0f, 0.0f);
                }
                continue;
            }

            // ---- Inner (voxelized) chunk: iterate all solid layers ----
            int min_s = job->col_min_solid[col][row];
            int max_s = job->col_max_solid[col][row];
            if (min_s > max_s) continue;

            // Precompute 6 hex corner directions (reused across layers)
            HMM_Vec3 hex_dirs[6];
            float hex_vx[6], hex_vz[6];
            for (int i = 0; i < 6; i++) {
                hex_vx[i] = lx + HEX_RADIUS * HEX_COS[i];
                hex_vz[i] = lz + HEX_RADIUS * HEX_SIN[i];
                hex_dirs[i] = vec3_normalize(
                    vec3_add(job->tangent_origin,
                        vec3_add(vec3_scale(job->tangent_east, hex_vx[i]),
                                 vec3_scale(job->tangent_north, hex_vz[i]))));
            }

            // Precompute wall outward directions for 6 edges
            HMM_Vec3 wall_outs[6];
            for (int dir = 0; dir < 6; dir++) {
                int edge = DIR_TO_EDGE[dir];
                int vi0 = edge;
                int vi1 = (edge + 1) % 6;
                float mx = (hex_vx[vi0] + hex_vx[vi1]) * 0.5f;
                float mz = (hex_vz[vi0] + hex_vz[vi1]) * 0.5f;
                wall_outs[dir] = vec3_normalize(
                    vec3_sub(
                        vec3_add(job->tangent_origin,
                            vec3_add(vec3_scale(job->tangent_east, mx),
                                     vec3_scale(job->tangent_north, mz))),
                        vec3_add(job->tangent_origin,
                            vec3_add(vec3_scale(job->tangent_east, lx),
                                     vec3_scale(job->tangent_north, lz)))
                    ));
            }

            for (int l = min_s; l <= max_s; l++) {
                uint8_t vtype = HEX_VOXEL(job->voxels, col, row, l);
                if (is_transparent(vtype)) continue;

                int world_layer = l + job->base_layer;
                float bot_r = job->planet_radius + world_layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
                float top_r = bot_r + HEX_HEIGHT;

                // --- TOP FACE: emit if layer above is transparent ---
                uint8_t above = job_get_voxel(job, col, row, l + 1);
                if (is_transparent(above)) {
                    // Column sky light from the air block above this face
                    float face_sky = read_sky_light(job, col, row, l + 1);
                    if (face_sky < 0.05f) face_sky = 0.05f;

                    int atlas = voxel_top_atlas((VoxelType)vtype);
                    HMM_Vec3 verts[6];
                    HexUV uvs[6];
                    for (int i = 0; i < 6; i++) {
                        verts[i] = vec3_scale(hex_dirs[i], top_r);
                        uvs[i] = hex_top_uv(hex_vx[i], hex_vz[i], lx, lz, atlas);
                    }
                    float face_torch = read_torch_light(job, col, row, l + 1);
                    for (int i = 0; i < 4; i++) {
                        hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                            verts[0], verts[i + 1], verts[i + 2],
                            local_up, slope_normal, uvs[0], uvs[i + 1], uvs[i + 2],
                            col_color, face_sky, face_torch);
                    }
                }

                // --- BOTTOM FACE: emit if layer below is transparent ---
                uint8_t below = job_get_voxel(job, col, row, l - 1);
                if (is_transparent(below)) {
                    // Column sky light from the air block below, with mild ceiling penalty
                    float face_sky = read_sky_light(job, col, row, l - 1) * 0.85f;
                    if (face_sky < 0.05f) face_sky = 0.05f;

                    int atlas = voxel_bottom_atlas((VoxelType)vtype);
                    HMM_Vec3 verts[6];
                    HexUV uvs[6];
                    for (int i = 0; i < 6; i++) {
                        verts[i] = vec3_scale(hex_dirs[i], bot_r);
                        uvs[i] = hex_top_uv(hex_vx[i], hex_vz[i], lx, lz, atlas);
                    }
                    // Reversed winding for bottom face
                    float face_torch = read_torch_light(job, col, row, l - 1);
                    for (int i = 0; i < 4; i++) {
                        hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                            verts[0], verts[i + 2], verts[i + 1],
                            local_down, local_down, uvs[0], uvs[i + 2], uvs[i + 1],
                            col_color, face_sky, face_torch);
                    }
                }

                // --- SIDE FACES: 6 hex neighbors ---
                for (int dir = 0; dir < 6; dir++) {
                    uint8_t nv = job_get_neighbor_voxel(job, &continental, &mountain, &warp, &detail,
                                                         col, row, l, dir);
                    if (!is_transparent(nv)) continue;

                    // Column sky light from the air block the side faces into
                    float face_sky;
                    float face_torch;
                    {
                        int ncol_g, nrow_g;
                        hex_neighbor(job->cx * HEX_CHUNK_SIZE + col,
                                     job->cz * HEX_CHUNK_SIZE + row, dir, &ncol_g, &nrow_g);
                        int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
                        int nr_local = nrow_g - job->cz * HEX_CHUNK_SIZE;
                        if (nc >= 0 && nc < HEX_CHUNK_SIZE &&
                            nr_local >= 0 && nr_local < HEX_CHUNK_SIZE) {
                            face_sky = read_sky_light(job, nc, nr_local, l);
                            face_torch = read_torch_light(job, nc, nr_local, l);
                        } else {
                            // Cross-chunk: estimate from neighbor voxel data
                            int world_layer = l + job->base_layer;
                            face_sky = estimate_nb_sky(job, ncol_g, nrow_g, world_layer);
                            face_torch = estimate_nb_torch(job, ncol_g, nrow_g, world_layer);
                        }
                    }
                    if (face_sky < 0.05f) face_sky = 0.05f;

                    int edge = DIR_TO_EDGE[dir];
                    int vi0 = edge;
                    int vi1 = (edge + 1) % 6;

                    HMM_Vec3 p0b = vec3_scale(hex_dirs[vi0], bot_r);
                    HMM_Vec3 p1b = vec3_scale(hex_dirs[vi1], bot_r);
                    HMM_Vec3 p0t = vec3_scale(hex_dirs[vi0], top_r);
                    HMM_Vec3 p1t = vec3_scale(hex_dirs[vi1], top_r);

                    int wall_atlas = voxel_side_atlas((VoxelType)vtype);
                    HexUV wuv00 = { ((float)wall_atlas + 0.0f) / (float)HEX_ATLAS_TILES, 1.0f };
                    HexUV wuv10 = { ((float)wall_atlas + 1.0f) / (float)HEX_ATLAS_TILES, 1.0f };
                    HexUV wuv01 = { ((float)wall_atlas + 0.0f) / (float)HEX_ATLAS_TILES, 0.0f };
                    HexUV wuv11 = { ((float)wall_atlas + 1.0f) / (float)HEX_ATLAS_TILES, 0.0f };

                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        p0b, p1b, p1t, wall_outs[dir], local_up, wuv00, wuv10, wuv11,
                        col_color, face_sky, face_torch);
                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        p0b, p1t, p0t, wall_outs[dir], local_up, wuv00, wuv11, wuv01,
                        col_color, face_sky, face_torch);
                }
            }
        }
    }

    // Phase 3b: Generate hex prism outline wireframe (physics debug)
    job->wire_verts = NULL;
    job->wire_count = 0;
    job->wire_cap = 0;

    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            int min_s = job->col_min_solid[col][row];
            int max_s = job->col_max_solid[col][row];
            if (min_s > max_s) continue;  // all-air column

            int gcol = job->cx * HEX_CHUNK_SIZE + col;
            int grow = job->cz * HEX_CHUNK_SIZE + row;

            float lx, lz;
            hex_local_pos(gcol, grow, &lx, &lz);

            HMM_Vec3 hex_wire_dirs[6];
            for (int i = 0; i < 6; i++) {
                float vx = lx + HEX_RADIUS * HEX_COS[i];
                float vz = lz + HEX_RADIUS * HEX_SIN[i];
                hex_wire_dirs[i] = vec3_normalize(
                    vec3_add(job->tangent_origin,
                        vec3_add(vec3_scale(job->tangent_east, vx),
                                 vec3_scale(job->tangent_north, vz))));
            }

            for (int l = min_s; l <= max_s; l++) {
                uint8_t vtype = HEX_VOXEL(job->voxels, col, row, l);
                if (is_transparent(vtype)) continue;

                int world_layer = l + job->base_layer;
                float bot_r = job->planet_radius + world_layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
                float top_r = bot_r + HEX_HEIGHT;

                // Top face outline: hexagon at top_r
                uint8_t above = (l + 1 < HEX_CHUNK_LAYERS)
                    ? HEX_VOXEL(job->voxels, col, row, l + 1) : VOXEL_AIR;
                if (is_transparent(above)) {
                    for (int i = 0; i < 6; i++) {
                        HMM_Vec3 a = vec3_scale(hex_wire_dirs[i], top_r);
                        HMM_Vec3 b = vec3_scale(hex_wire_dirs[(i + 1) % 6], top_r);
                        wire_push_line(&job->wire_verts, &job->wire_count,
                                       &job->wire_cap, a, b);
                    }
                }

                // Bottom face outline: hexagon at bot_r
                uint8_t below = (l - 1 >= 0)
                    ? HEX_VOXEL(job->voxels, col, row, l - 1) : VOXEL_AIR;
                if (is_transparent(below)) {
                    for (int i = 0; i < 6; i++) {
                        HMM_Vec3 a = vec3_scale(hex_wire_dirs[i], bot_r);
                        HMM_Vec3 b = vec3_scale(hex_wire_dirs[(i + 1) % 6], bot_r);
                        wire_push_line(&job->wire_verts, &job->wire_count,
                                       &job->wire_cap, a, b);
                    }
                }

                // Vertical edges for exposed sides
                for (int dir = 0; dir < 6; dir++) {
                    int ncol_g, nrow_g;
                    hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
                    int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
                    int nr = nrow_g - job->cz * HEX_CHUNK_SIZE;
                    bool exposed = true;
                    if (nc >= 0 && nc < HEX_CHUNK_SIZE && nr >= 0 && nr < HEX_CHUNK_SIZE) {
                        exposed = is_transparent(HEX_VOXEL(job->voxels, nc, nr, l));
                    }
                    if (exposed) {
                        int edge = DIR_TO_EDGE[dir];
                        int vi0 = edge;
                        int vi1 = (edge + 1) % 6;
                        // Two vertical edges for this side
                        wire_push_line(&job->wire_verts, &job->wire_count,
                                       &job->wire_cap,
                                       vec3_scale(hex_wire_dirs[vi0], bot_r),
                                       vec3_scale(hex_wire_dirs[vi0], top_r));
                        wire_push_line(&job->wire_verts, &job->wire_count,
                                       &job->wire_cap,
                                       vec3_scale(hex_wire_dirs[vi1], bot_r),
                                       vec3_scale(hex_wire_dirs[vi1], top_r));
                    }
                }
            }
        }
    }

    // Phase 4: Apply floating origin offset
    for (int i = 0; i < job->vertex_count; i++) {
        job->vertices[i].pos[0] -= job->origin[0];
        job->vertices[i].pos[1] -= job->origin[1];
        job->vertices[i].pos[2] -= job->origin[2];
    }
    // Also offset wireframe vertices
    for (int i = 0; i < job->wire_count; i += 3) {
        job->wire_verts[i + 0] -= job->origin[0];
        job->wire_verts[i + 1] -= job->origin[1];
        job->wire_verts[i + 2] -= job->origin[2];
    }

    // Free light maps and neighbor data (only needed during mesh gen, not after)
    if (job->sky_map) {
        free(job->sky_map);
        job->sky_map = NULL;
    }
    if (job->torch_map) {
        free(job->torch_map);
        job->torch_map = NULL;
    }
    free_job_neighbors(job);

    double t_end = timer_ms();

    // Benchmark: print timing for edit-dirty remeshes (has_voxels=true)
    if (job->has_voxels) {
        printf("[HEX BENCH] Remesh chunk(%d,%d): voxels=%.1fms bedrock=%.1fms skylight=%.1fms mesh=%.1fms TOTAL=%.1fms (%d verts)\n",
               job->cx, job->cz,
               t_voxels - t_start,
               t_bedrock - t_voxels,
               t_skylight - t_bedrock,
               t_end - t_skylight,
               t_end - t_start,
               job->vertex_count);
        fflush(stdout);
    }

    job->completed = 1;
}

// Worker thread entry point
static void hex_mesh_worker(void* data) {
    HexMeshJob* job = (HexMeshJob*)data;
    generate_chunk_mesh(job);
}

// ---- Orphaned job cleanup ----
#define HEX_MAX_ORPHANS 128
static void* s_orphaned_jobs[HEX_MAX_ORPHANS];
static int s_orphan_count = 0;

static void orphan_job(void* job) {
    if (s_orphan_count < HEX_MAX_ORPHANS) {
        s_orphaned_jobs[s_orphan_count++] = job;
    }
}

static void sweep_orphans(void) {
    int write = 0;
    for (int i = 0; i < s_orphan_count; i++) {
        HexMeshJob* job = (HexMeshJob*)s_orphaned_jobs[i];
        if (job->completed) {
            if (job->vertices) free(job->vertices);
            if (job->voxels) free(job->voxels);
            if (job->sky_map) free(job->sky_map);
            if (job->torch_map) free(job->torch_map);
            if (job->wire_verts) free(job->wire_verts);
            for (int b = 0; b < 4; b++) { if (job->out_border_torch[b]) free(job->out_border_torch[b]); }
            free_job_neighbors(job);
            free(job);
        } else {
            s_orphaned_jobs[write++] = s_orphaned_jobs[i];
        }
    }
    s_orphan_count = write;
}

// ---- Chunk management ----

static int find_chunk(HexTerrain* ht, int cx, int cz) {
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (ht->chunks[i].active && ht->chunks[i].cx == cx && ht->chunks[i].cz == cz) {
            return i;
        }
    }
    return -1;
}

static int find_chunk_const(const HexTerrain* ht, int cx, int cz) {
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (ht->chunks[i].active && ht->chunks[i].cx == cx && ht->chunks[i].cz == cz) {
            return i;
        }
    }
    return -1;
}

static int find_free_chunk(HexTerrain* ht) {
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (!ht->chunks[i].active) return i;
    }
    return -1;
}

// Copy neighbor chunk voxel data into a mesh job for cross-chunk BFS and geometry.
// Called on the main thread before job submission.
static void populate_job_neighbors(const HexTerrain* ht, HexMeshJob* job) {
    size_t voxel_size = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;
    size_t border_size = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;
    int nb_count = 0;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dz = -1; dz <= 1; dz++) {
            if (dx == 0 && dz == 0) continue;
            int ncx = job->cx + dx;
            int ncz = job->cz + dz;
            int ni = find_chunk_const(ht, ncx, ncz);
            if (ni >= 0 && ht->chunks[ni].voxels) {
                job->nb_voxels[dx + 1][dz + 1] = (uint8_t*)malloc(voxel_size);
                memcpy(job->nb_voxels[dx + 1][dz + 1], ht->chunks[ni].voxels, voxel_size);
                job->nb_base_layer[dx + 1][dz + 1] = ht->chunks[ni].base_layer;
                nb_count++;

                // Copy neighbor's facing border torch values (cardinal neighbors only)
                if (dx == 0 || dz == 0) {
                    // Which border of the neighbor faces us?
                    // dx=-1: neighbor is west, their east border (border[1]) faces our west
                    // dx=+1: neighbor is east, their west border (border[0]) faces our east
                    // dz=-1: neighbor is south, their north border (border[3]) faces our south
                    // dz=+1: neighbor is north, their south border (border[2]) faces our north
                    int bi = -1;
                    if (dx == -1) bi = 1;
                    else if (dx == 1) bi = 0;
                    else if (dz == -1) bi = 3;
                    else if (dz == 1) bi = 2;

                    if (bi >= 0 && ht->chunks[ni].border_torch[bi]) {
                        job->nb_border_torch[dx+1][dz+1] = (uint8_t*)malloc(border_size);
                        memcpy(job->nb_border_torch[dx+1][dz+1],
                               ht->chunks[ni].border_torch[bi], border_size);
                    }
                }
            }
        }
    }
    job->nb_count = nb_count;
    hex_log("[NEIGHBORS] chunk(%d,%d) populated %d/8 neighbor voxel sets\n",
            job->cx, job->cz, nb_count);
}

static void free_chunk(HexChunk* chunk) {
    if (chunk->cpu_vertices) {
        free(chunk->cpu_vertices);
        chunk->cpu_vertices = NULL;
    }
    chunk->cpu_vertex_count = 0;

    if (chunk->gpu_buffer.id != SG_INVALID_ID) {
        sg_destroy_buffer(chunk->gpu_buffer);
        chunk->gpu_buffer.id = SG_INVALID_ID;
    }
    chunk->gpu_vertex_count = 0;

    // Free wireframe data
    if (chunk->cpu_wire_verts) {
        free(chunk->cpu_wire_verts);
        chunk->cpu_wire_verts = NULL;
    }
    chunk->cpu_wire_count = 0;
    if (chunk->wire_buf.id != SG_INVALID_ID) {
        sg_destroy_buffer(chunk->wire_buf);
        chunk->wire_buf = (sg_buffer){0};
    }
    chunk->wire_vertex_count = 0;

    if (chunk->voxels) {
        free(chunk->voxels);
        chunk->voxels = NULL;
    }

    // Free boundary torch light data
    for (int b = 0; b < 4; b++) {
        if (chunk->border_torch[b]) {
            free(chunk->border_torch[b]);
            chunk->border_torch[b] = NULL;
        }
    }

    if (chunk->pending_job) {
        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        if (job->completed) {
            if (job->vertices) free(job->vertices);
            if (job->voxels) free(job->voxels);
            if (job->sky_map) free(job->sky_map);
            if (job->torch_map) free(job->torch_map);
            if (job->wire_verts) free(job->wire_verts);
            for (int b = 0; b < 4; b++) { if (job->out_border_torch[b]) free(job->out_border_torch[b]); }
            free_job_neighbors(job);
            free(job);
        } else {
            orphan_job(job);
        }
        chunk->pending_job = NULL;
    }

    chunk->active = false;
    chunk->dirty = false;
    chunk->generating = false;
    chunk->torch_scanned = false;
}

// Update column min/max solid caches after voxel modification
static void update_column_cache(HexChunk* chunk, int col, int row) {
    if (!chunk->voxels) return;
    int min_s = HEX_CHUNK_LAYERS;
    int max_s = -1;
    for (int l = 0; l < HEX_CHUNK_LAYERS; l++) {
        if (HEX_VOXEL(chunk->voxels, col, row, l) != VOXEL_AIR) {
            if (l < min_s) min_s = l;
            if (l > max_s) max_s = l;
        }
    }
    chunk->col_min_solid[col][row] = (int16_t)min_s;
    chunk->col_max_solid[col][row] = (int16_t)max_s;
}

// ---- Public API ----

void hex_terrain_init(HexTerrain* ht, float planet_radius, float layer_thickness,
                      int sea_level, int seed, JobSystem* jobs) {
    memset(ht, 0, sizeof(HexTerrain));
    ht->planet_radius = planet_radius;
    ht->layer_thickness = layer_thickness;
    ht->sea_level = sea_level;
    ht->seed = seed;
    ht->jobs = jobs;

    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        ht->chunks[i].active = false;
        ht->chunks[i].gpu_buffer.id = SG_INVALID_ID;
        ht->chunks[i].pending_job = NULL;
        ht->chunks[i].voxels = NULL;
    }

    // Initialize voxel edit persistence
    edit_cache_init(&ht->edits, seed, planet_radius, "cache/edits");

    printf("[HEX] 3D voxel terrain initialized: radius=%.0f, range=%.0fm, chunk=%dx%dx%d\n",
        planet_radius, HEX_RANGE, HEX_CHUNK_SIZE, HEX_CHUNK_SIZE, HEX_CHUNK_LAYERS);
    fflush(stdout);
}

void hex_terrain_destroy(HexTerrain* ht) {
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (ht->chunks[i].pending_job) {
            HexMeshJob* job = (HexMeshJob*)ht->chunks[i].pending_job;
            while (!job->completed) { /* spin-wait */ }
        }
        if (ht->chunks[i].active) {
            free_chunk(&ht->chunks[i]);
        }
    }
    for (int i = 0; i < s_orphan_count; i++) {
        HexMeshJob* job = (HexMeshJob*)s_orphaned_jobs[i];
        while (!job->completed) { /* spin-wait */ }
        if (job->vertices) free(job->vertices);
        if (job->voxels) free(job->voxels);
        if (job->sky_map) free(job->sky_map);
        if (job->torch_map) free(job->torch_map);
        if (job->wire_verts) free(job->wire_verts);
        for (int b = 0; b < 4; b++) { if (job->out_border_torch[b]) free(job->out_border_torch[b]); }
        free_job_neighbors(job);
        free(job);
    }
    s_orphan_count = 0;
    edit_cache_flush_all(&ht->edits);
    edit_cache_destroy(&ht->edits);
}

#define FRAME_REANCHOR_THRESHOLD 10000.0f

void hex_terrain_update(HexTerrain* ht, HMM_Vec3 camera_pos,
                        const double world_origin[3]) {
    sweep_orphans();

    ht->camera_pos = camera_pos;
    ht->world_origin[0] = world_origin[0];
    ht->world_origin[1] = world_origin[1];
    ht->world_origin[2] = world_origin[2];

    // Disable hex terrain when camera is high above ground
    float cam_r = sqrtf(vec3_dot(camera_pos, camera_pos));
    HMM_Vec3 cam_dir = vec3_scale(camera_pos, 1.0f / fmaxf(cam_r, 1.0f));
    fnl_state _c = ht_create_continental_noise(ht->seed);
    fnl_state _m = ht_create_mountain_noise(ht->seed);
    fnl_state _w = ht_create_warp_noise(ht->seed);
    fnl_state _d = ht_create_detail_noise(ht->seed);
    float ground_h = ht_sample_height_m(&_c, &_m, &_w, &_d, cam_dir);
    float ground_r = ht->planet_radius + fmaxf(ground_h, TERRAIN_SEA_LEVEL_M);
    float altitude = cam_r - ground_r;
    if (altitude > HEX_MAX_DRAW_ALT) {
        for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
            if (ht->chunks[i].active) free_chunk(&ht->chunks[i]);
        }
        ht->frame_valid = false;
        ht->active_count = 0;
        return;
    }

    // Update voxel edit cache (preload nearby sectors, periodic flush)
    edit_cache_update(&ht->edits, cam_dir);

    // Check if tangent frame needs re-anchoring
    bool reanchor = false;
    if (!ht->frame_valid) {
        reanchor = true;
    } else {
        HMM_Vec3 drift = vec3_sub(camera_pos, ht->frame_anchor);
        float drift_sq = vec3_dot(drift, drift);
        if (drift_sq > FRAME_REANCHOR_THRESHOLD * FRAME_REANCHOR_THRESHOLD) {
            reanchor = true;
        }
    }

    if (reanchor) {
        compute_tangent_frame(ht, camera_pos);
        ht->frame_anchor = camera_pos;
        ht->frame_valid = true;
        for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
            if (ht->chunks[i].active) free_chunk(&ht->chunks[i]);
        }
        ht->initial_gen_pending = 0;
        ht->initial_load_complete = false;
    }

    // Process completed mesh generation jobs
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (!chunk->active || !chunk->generating || !chunk->pending_job) continue;

        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        if (!job->completed) continue;

        // If chunk was re-dirtied during generation (break/place happened),
        // keep the chunk's current voxels (which have edits) and discard
        // the job's stale voxel data. Still take the mesh as a visual fallback
        // until the re-generation completes.
        bool re_dirtied = chunk->dirty;

        // Transfer results: mesh (always — better than nothing)
        chunk->cpu_vertices = job->vertices;
        chunk->cpu_vertex_count = job->vertex_count;
        job->vertices = NULL;

        // Transfer results: wireframe
        if (chunk->cpu_wire_verts) free(chunk->cpu_wire_verts);
        chunk->cpu_wire_verts = job->wire_verts;
        chunk->cpu_wire_count = job->wire_count;
        job->wire_verts = NULL;

        // Transfer results: voxel data — only if chunk wasn't edited during generation
        if (re_dirtied) {
            // Chunk voxels were modified by break/place; discard job's stale copy
            free(job->voxels);
            job->voxels = NULL;
        } else {
            if (chunk->voxels) free(chunk->voxels);
            chunk->voxels = job->voxels;
            job->voxels = NULL;
            chunk->base_layer = job->base_layer;
            memcpy(chunk->col_min_solid, job->col_min_solid, sizeof(chunk->col_min_solid));
            memcpy(chunk->col_max_solid, job->col_max_solid, sizeof(chunk->col_max_solid));
            chunk->torch_scanned = false;  // Re-scan for torch instances
        }

        // Transfer boundary torch light values (for cross-chunk propagation)
        for (int b = 0; b < 4; b++) {
            if (chunk->border_torch[b]) { free(chunk->border_torch[b]); chunk->border_torch[b] = NULL; }
            chunk->border_torch[b] = job->out_border_torch[b];
            job->out_border_torch[b] = NULL;
        }

        // Track initial-gen completion for coordinated second pass
        if (job->is_initial_gen) {
            ht->initial_gen_pending--;
        }

        free_job_neighbors(job);
        free(job);
        chunk->pending_job = NULL;
        chunk->generating = false;
        // Don't clear dirty — it was already cleared at submission time.
        // If it's true now, something re-dirtied it and we need another pass.
    }

    // Coordinated second pass: once ALL initial-gen chunks have completed
    // and no more chunks need activation, mark all active chunks dirty.
    // This ensures the second pass has full neighbor data for cross-chunk BFS.
    if (!ht->initial_load_complete && ht->initial_gen_pending == 0) {
        // Check that at least one active chunk exists with voxel data,
        // and no active chunk still lacks voxel data (waiting for generation).
        int active_with_voxels = 0;
        bool any_missing = false;
        for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
            if (!ht->chunks[i].active) continue;
            if (ht->chunks[i].voxels) {
                active_with_voxels++;
            } else {
                any_missing = true;
                break;
            }
        }
        if (active_with_voxels > 0 && !any_missing) {
            for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
                if (ht->chunks[i].active && !ht->chunks[i].generating) {
                    ht->chunks[i].dirty = true;
                }
            }
            ht->initial_load_complete = true;
            hex_log("[SECOND_PASS] All %d initial chunks complete — marking all dirty for coordinated rebuild\n", active_with_voxels);
        }
    }

    // Determine chunk loading range
    HMM_Vec3 cam_offset = vec3_sub(camera_pos, ht->tangent_origin);
    float cam_lx = vec3_dot(cam_offset, ht->tangent_east);
    float cam_lz = vec3_dot(cam_offset, ht->tangent_north);

    int cam_cx = (int)floorf(cam_lx / HEX_CHUNK_WIDTH);
    int cam_cz = (int)floorf(cam_lz / HEX_CHUNK_DEPTH);
    int range_cx = (int)ceilf(HEX_RANGE / HEX_CHUNK_WIDTH) + 1;
    int range_cz = (int)ceilf(HEX_RANGE / HEX_CHUNK_DEPTH) + 1;

    // Deactivate out-of-range chunks
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (!chunk->active) continue;
        if (abs(chunk->cx - cam_cx) > range_cx || abs(chunk->cz - cam_cz) > range_cz) {
            free_chunk(chunk);
        }
    }

    // Activate/generate chunks within range
    int jobs_submitted = 0;
    int new_activations = 0;
    for (int cx = cam_cx - range_cx; cx <= cam_cx + range_cx; cx++) {
        for (int cz = cam_cz - range_cz; cz <= cam_cz + range_cz; cz++) {
            if (new_activations >= HEX_MAX_ACTIVATIONS || jobs_submitted >= 8) break;

            float chunk_center_x = (cx * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE / 2) * HEX_COL_SPACING;
            float chunk_center_z = (cz * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE / 2) * HEX_ROW_SPACING;
            float dx = chunk_center_x - cam_lx;
            float dz = chunk_center_z - cam_lz;
            float dist = sqrtf(dx * dx + dz * dz);
            if (dist > HEX_RANGE + HEX_CHUNK_WIDTH) continue;

            int idx = find_chunk(ht, cx, cz);

            if (idx >= 0) {
                HexChunk* chunk = &ht->chunks[idx];
                if (chunk->gpu_vertex_count == 0 && !chunk->generating &&
                    chunk->cpu_vertex_count == 0) {
                    chunk->dirty = true;
                }
                // Force voxel prism mode if stuck in transition
                if (chunk->is_transition) {
                    chunk->is_transition = false;
                    if (!chunk->generating) chunk->dirty = true;
                }
            } else {
                if (new_activations >= HEX_MAX_ACTIVATIONS) continue;
                idx = find_free_chunk(ht);
                if (idx < 0) continue;

                HexChunk* chunk = &ht->chunks[idx];
                memset(chunk, 0, sizeof(HexChunk));
                chunk->cx = cx;
                chunk->cz = cz;
                chunk->active = true;
                chunk->dirty = true;
                chunk->is_transition = false;  // Always voxel prisms
                chunk->gpu_buffer.id = SG_INVALID_ID;
                chunk->wire_buf.id = SG_INVALID_ID;
                chunk->voxels = NULL;
                new_activations++;
            }

            // Submit mesh gen job for dirty chunks
            if (idx >= 0 && idx < HEX_MAX_CHUNKS && jobs_submitted < 8) {
                HexChunk* chunk = &ht->chunks[idx];
                if (chunk->dirty && !chunk->generating) {
                    size_t voxel_size = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;

                    HexMeshJob* job = (HexMeshJob*)calloc(1, sizeof(HexMeshJob));
                    job->cx = chunk->cx;
                    job->cz = chunk->cz;
                    job->planet_radius = ht->planet_radius;
                    job->seed = ht->seed;
                    job->tangent_origin = ht->tangent_origin;
                    job->tangent_up = ht->tangent_up;
                    job->tangent_east = ht->tangent_east;
                    job->tangent_north = ht->tangent_north;
                    job->origin[0] = (float)ht->world_origin[0];
                    job->origin[1] = (float)ht->world_origin[1];
                    job->origin[2] = (float)ht->world_origin[2];
                    job->chunk_index = idx;
                    job->is_transition = chunk->is_transition;
                    job->edit_cache = &ht->edits;

                    // If chunk already has voxel data (dirty from break/place),
                    // copy it so the worker rebuilds mesh from existing voxels
                    if (chunk->voxels) {
                        job->voxels = (uint8_t*)malloc(voxel_size);
                        memcpy(job->voxels, chunk->voxels, voxel_size);
                        job->base_layer = chunk->base_layer;
                        job->has_voxels = true;
                        memcpy(job->col_min_solid, chunk->col_min_solid, sizeof(job->col_min_solid));
                        memcpy(job->col_max_solid, chunk->col_max_solid, sizeof(job->col_max_solid));
                    } else {
                        job->voxels = (uint8_t*)malloc(voxel_size);
                        job->has_voxels = false;
                        job->is_initial_gen = true;
                        // Will be filled by generate_chunk_mesh
                    }

                    job->vertices = NULL;
                    job->vertex_count = 0;
                    job->completed = 0;
                    populate_job_neighbors(ht, job);

                    if (job_system_try_submit(ht->jobs, hex_mesh_worker, job)) {
                        chunk->pending_job = job;
                        chunk->generating = true;
                        chunk->dirty = false;
                        chunk->edit_dirty = false;
                        jobs_submitted++;
                        if (job->is_initial_gen) ht->initial_gen_pending++;
                    } else {
                        free(job->voxels);
                        free_job_neighbors(job);
                        free(job);
                    }
                }
            }
        }
        if (new_activations >= HEX_MAX_ACTIVATIONS || jobs_submitted >= 8) break;
    }

    // Second pass: submit jobs for remaining dirty chunks.
    // Priority 1: edit_dirty chunks (player break/place — must be immediate)
    // Priority 2: regular dirty chunks
    for (int priority = 0; priority < 2 && jobs_submitted < 16; priority++) {
        for (int i = 0; i < HEX_MAX_CHUNKS && jobs_submitted < 16; i++) {
            HexChunk* chunk = &ht->chunks[i];
            if (!chunk->active || chunk->generating) continue;

            bool should_submit = false;
            if (priority == 0) {
                should_submit = chunk->edit_dirty;
            } else {
                should_submit = chunk->dirty && !chunk->edit_dirty;
            }
            if (!should_submit) continue;

            size_t voxel_size = (size_t)HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS;
            HexMeshJob* job = (HexMeshJob*)calloc(1, sizeof(HexMeshJob));
            job->cx = chunk->cx;
            job->cz = chunk->cz;
            job->planet_radius = ht->planet_radius;
            job->seed = ht->seed;
            job->tangent_origin = ht->tangent_origin;
            job->tangent_up = ht->tangent_up;
            job->tangent_east = ht->tangent_east;
            job->tangent_north = ht->tangent_north;
            job->origin[0] = (float)ht->world_origin[0];
            job->origin[1] = (float)ht->world_origin[1];
            job->origin[2] = (float)ht->world_origin[2];
            job->chunk_index = i;
            job->is_transition = chunk->is_transition;
            job->edit_cache = &ht->edits;

            if (chunk->voxels) {
                job->voxels = (uint8_t*)malloc(voxel_size);
                memcpy(job->voxels, chunk->voxels, voxel_size);
                job->base_layer = chunk->base_layer;
                job->has_voxels = true;
                memcpy(job->col_min_solid, chunk->col_min_solid, sizeof(job->col_min_solid));
                memcpy(job->col_max_solid, chunk->col_max_solid, sizeof(job->col_max_solid));
            } else {
                job->voxels = (uint8_t*)malloc(voxel_size);
                job->has_voxels = false;
                job->is_initial_gen = true;
            }

            job->vertices = NULL;
            job->vertex_count = 0;
            job->completed = 0;
            populate_job_neighbors(ht, job);

            if (job_system_try_submit(ht->jobs, hex_mesh_worker, job)) {
                hex_log("[JOB_SUBMIT] chunk(%d,%d) priority=%d edit_dirty=%d has_voxels=%d\n",
                        chunk->cx, chunk->cz, priority, chunk->edit_dirty, job->has_voxels);
                chunk->pending_job = job;
                chunk->generating = true;
                chunk->dirty = false;
                chunk->edit_dirty = false;
                jobs_submitted++;
                if (job->is_initial_gen) ht->initial_gen_pending++;
            } else {
                free(job->voxels);
                free_job_neighbors(job);
                free(job);
            }
        }
    }

    // Watchdog: force-dirty any active chunks stuck with no mesh
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (chunk->active && !chunk->generating && !chunk->dirty &&
            chunk->gpu_vertex_count == 0 && chunk->cpu_vertex_count == 0) {
            chunk->dirty = true;
        }
    }

    int active = 0;
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (ht->chunks[i].active) active++;
    }
    ht->active_count = active;
}

void hex_terrain_upload_meshes(HexTerrain* ht) {
    int uploads = 0;
    for (int i = 0; i < HEX_MAX_CHUNKS && uploads < HEX_MAX_UPLOADS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (!chunk->active) continue;
        if (!chunk->cpu_vertices || chunk->cpu_vertex_count == 0) continue;

        if (chunk->gpu_buffer.id != SG_INVALID_ID) {
            sg_destroy_buffer(chunk->gpu_buffer);
        }

        chunk->gpu_buffer = sg_make_buffer(&(sg_buffer_desc){
            .data = (sg_range){
                chunk->cpu_vertices,
                (size_t)chunk->cpu_vertex_count * sizeof(HexVertex)
            },
            .label = "hex-chunk-vertices",
        });

        if (chunk->gpu_buffer.id == SG_INVALID_ID) continue;

        chunk->gpu_vertex_count = chunk->cpu_vertex_count;
        uploads++;

        free(chunk->cpu_vertices);
        chunk->cpu_vertices = NULL;
        chunk->cpu_vertex_count = 0;

        // Upload wireframe alongside mesh
        if (chunk->wire_buf.id != SG_INVALID_ID) {
            sg_destroy_buffer(chunk->wire_buf);
            chunk->wire_buf = (sg_buffer){0};
        }
        if (chunk->cpu_wire_verts && chunk->cpu_wire_count > 0) {
            chunk->wire_buf = sg_make_buffer(&(sg_buffer_desc){
                .data = (sg_range){
                    chunk->cpu_wire_verts,
                    (size_t)chunk->cpu_wire_count * sizeof(float)
                },
                .label = "hex-wire-verts",
            });
            chunk->wire_vertex_count = chunk->cpu_wire_count / 3;
        }
        if (chunk->cpu_wire_verts) {
            free(chunk->cpu_wire_verts);
            chunk->cpu_wire_verts = NULL;
            chunk->cpu_wire_count = 0;
        }
    }
}

void hex_terrain_render(HexTerrain* ht, sg_pipeline pip,
                        sg_view atlas_view, sg_sampler atlas_smp) {
    (void)pip;
    sg_bindings bind = {0};
    bind.views[0] = atlas_view;
    bind.samplers[0] = atlas_smp;
    int total_verts = 0;
    int rendered = 0;

    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (!chunk->active) continue;
        if (chunk->gpu_vertex_count == 0) continue;
        if (chunk->gpu_buffer.id == SG_INVALID_ID) continue;

        bind.vertex_buffers[0] = chunk->gpu_buffer;
        sg_apply_bindings(&bind);
        sg_draw(0, chunk->gpu_vertex_count, 1);
        total_verts += chunk->gpu_vertex_count;
        rendered++;
    }

    ht->total_vertex_count = total_verts;
    ht->chunks_rendered = rendered;
}

bool hex_terrain_covers_position(const HexTerrain* ht, HMM_Vec3 world_pos) {
    HMM_Vec3 diff = vec3_sub(world_pos, ht->camera_pos);
    float dist_sq = vec3_dot(diff, diff);
    return dist_sq < HEX_RANGE * HEX_RANGE;
}

float hex_terrain_get_range(void) {
    return HEX_RANGE;
}

float hex_terrain_effective_range(const HexTerrain* ht) {
    if (!ht->frame_valid) return 0.0f;

    // Find the nearest unmeshed chunk — any hole in coverage means we can't
    // suppress LOD patches beyond that distance.
    float min_unmeshed_dist = HEX_RANGE * 2.0f; // sentinel = very far
    int total = 0, meshed = 0;

    // Camera position in tangent-plane coords (computed once)
    HMM_Vec3 cam_off = vec3_sub(ht->camera_pos, ht->tangent_origin);
    float cam_tx = vec3_dot(cam_off, ht->tangent_east);
    float cam_tz = vec3_dot(cam_off, ht->tangent_north);

    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (!ht->chunks[i].active) continue;
        total++;

        // Chunk center in tangent-plane meters
        float cx_m = (ht->chunks[i].cx * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE * 0.5f) * HEX_COL_SPACING;
        float cz_m = (ht->chunks[i].cz * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE * 0.5f) * HEX_ROW_SPACING;
        float dx = cx_m - cam_tx;
        float dz = cz_m - cam_tz;
        float dist = sqrtf(dx * dx + dz * dz);

        if (ht->chunks[i].gpu_vertex_count > 0) {
            meshed++;
        } else if (dist < min_unmeshed_dist) {
            min_unmeshed_dist = dist;
        }
    }

    if (total == 0) return 0.0f;

    // All chunks meshed: suppress the full inner range
    if (meshed == total) return HEX_INNER_RANGE;

    // Suppress up to the nearest unmeshed chunk minus a safety margin
    float safe_range = min_unmeshed_dist - HEX_CHUNK_WIDTH * 2.0f;
    if (safe_range < 100.0f) return 0.0f;  // Too close, don't suppress at all
    return fminf(safe_range, HEX_INNER_RANGE);
}

// ---- Raycast (3D voxel-aware) ----

HexHitResult hex_terrain_raycast(const HexTerrain* ht, HMM_Vec3 ray_origin,
                                  HMM_Vec3 ray_dir, float max_dist) {
    HexHitResult result = { .valid = false };
    if (!ht->frame_valid) return result;

    float step = 0.1f;
    int prev_gcol = -99999, prev_grow = -99999, prev_layer = -99999;

    for (float t = 0.0f; t <= max_dist; t += step) {
        HMM_Vec3 world_pt = vec3_add(ray_origin, vec3_scale(ray_dir, t));
        float point_r = sqrtf(vec3_dot(world_pt, world_pt));

        // Spherical inverse: world position → tangent-plane hex coords
        float lx, lz;
        world_to_tangent(ht, world_pt, &lx, &lz);

        int gcol, grow;
        pixel_to_hex(lx, lz, &gcol, &grow);

        // Compute layer from radius
        float altitude = point_r - ht->planet_radius;
        int world_layer = (int)floorf(altitude / HEX_HEIGHT);

        // Skip if same voxel
        if (gcol == prev_gcol && grow == prev_grow && world_layer == prev_layer) continue;

        // Look up voxel
        uint8_t vtype = hex_terrain_get_voxel(ht, gcol, grow, world_layer);

        if (vtype != VOXEL_AIR) {
            // Determine hit face
            int face;
            int p_gcol, p_grow, p_layer;
            if (world_layer != prev_layer && gcol == prev_gcol && grow == prev_grow) {
                // Vertical transition
                if (world_layer > prev_layer) {
                    face = 1; // bottom face
                    p_gcol = gcol; p_grow = grow; p_layer = world_layer - 1;
                } else {
                    face = 0; // top face
                    p_gcol = gcol; p_grow = grow; p_layer = world_layer + 1;
                }
            } else if (gcol != prev_gcol || grow != prev_grow) {
                // Horizontal transition: side face
                face = 2; // generic side
                p_gcol = prev_gcol; p_grow = prev_grow; p_layer = world_layer;
            } else {
                // First step hit (no previous)
                face = 0;
                p_gcol = gcol; p_grow = grow; p_layer = world_layer + 1;
            }

            // HIT accepted
            int cx = (int)floorf((float)gcol / HEX_CHUNK_SIZE);
            int cz = (int)floorf((float)grow / HEX_CHUNK_SIZE);
            int local_col = gcol - cx * HEX_CHUNK_SIZE;
            int local_row = grow - cz * HEX_CHUNK_SIZE;
            if (local_col < 0) { local_col += HEX_CHUNK_SIZE; cx--; }
            if (local_row < 0) { local_row += HEX_CHUNK_SIZE; cz--; }

            result.valid = true;
            result.chunk_index = find_chunk_const(ht, cx, cz);
            result.col = local_col;
            result.row = local_row;
            result.gcol = gcol;
            result.grow = grow;
            result.layer = world_layer;
            result.type = vtype;
            result.face = face;
            result.place_gcol = p_gcol;
            result.place_grow = p_grow;
            result.place_layer = p_layer;
            return result;
        }

        prev_gcol = gcol;
        prev_grow = grow;
        prev_layer = world_layer;
    }

    return result;
}

// ---- Break / Place (per-voxel) ----

bool hex_terrain_break(HexTerrain* ht, const HexHitResult* hit) {
    if (!hit->valid) return false;

    // Re-derive chunk from global coords (hit->chunk_index may be stale)
    int cx = (int)floorf((float)hit->gcol / HEX_CHUNK_SIZE);
    int cz = (int)floorf((float)hit->grow / HEX_CHUNK_SIZE);
    int local_col = hit->gcol - cx * HEX_CHUNK_SIZE;
    int local_row = hit->grow - cz * HEX_CHUNK_SIZE;
    if (local_col < 0) { local_col += HEX_CHUNK_SIZE; cx--; }
    if (local_row < 0) { local_row += HEX_CHUNK_SIZE; cz--; }

    int chunk_idx = find_chunk(ht, cx, cz);
    if (chunk_idx < 0) return false;

    HexChunk* chunk = &ht->chunks[chunk_idx];
    if (!chunk->active || !chunk->voxels) return false;

    int local_layer = hit->layer - chunk->base_layer;
    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) return false;

    // Bedrock is unbreakable
    uint8_t broken_type = HEX_VOXEL(chunk->voxels, local_col, local_row, local_layer);
    if (broken_type == VOXEL_BEDROCK) {
        float altitude_m = (float)hit->layer * HEX_HEIGHT;
        float base_alt_m = (float)chunk->base_layer * HEX_HEIGHT;
        hex_log("[BREAK] BLOCKED: bedrock at gcol=%d grow=%d world_layer=%d (alt=%.0fm) | "
                "chunk(%d,%d) local(%d,%d) local_layer=%d base_layer=%d (base_alt=%.0fm) | "
                "bedrock_range=[0..%d] col_min=%d col_max=%d\n",
                hit->gcol, hit->grow, hit->layer, altitude_m,
                cx, cz, local_col, local_row, local_layer,
                chunk->base_layer, base_alt_m,
                HEX_BEDROCK_LAYERS - 1,
                chunk->col_min_solid[local_col][local_row],
                chunk->col_max_solid[local_col][local_row]);
        return false;
    }

    hex_log("[BREAK] type=%d at gcol=%d grow=%d layer=%d | chunk(%d,%d) local(%d,%d)\n",
            broken_type, hit->gcol, hit->grow, hit->layer, cx, cz, local_col, local_row);

    HEX_VOXEL(chunk->voxels, local_col, local_row, local_layer) = VOXEL_AIR;
    update_column_cache(chunk, local_col, local_row);
    chunk->dirty = true;
    chunk->edit_dirty = true;

    // Persist edit to disk-backed sector cache
    {
        double ux, uy, uz;
        gcol_grow_to_unit(ht, hit->gcol, hit->grow, &ux, &uy, &uz);
        edit_cache_record(&ht->edits, ux, uy, uz, hit->layer, VOXEL_AIR);
    }

    // Always dirty all adjacent chunks: any block edit affects cross-chunk
    // BFS lighting (sky light range 32, torch range 15) and boundary geometry.
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (!ht->chunks[i].active || i == chunk_idx) continue;
        int dx = ht->chunks[i].cx - cx;
        int dz = ht->chunks[i].cz - cz;
        if (dx >= -1 && dx <= 1 && dz >= -1 && dz <= 1) {
            ht->chunks[i].dirty = true;
        }
    }
    return true;
}

bool hex_terrain_place(HexTerrain* ht, const HexHitResult* hit, uint8_t voxel_type) {
    if (!hit->valid) return false;

    int place_gcol = hit->place_gcol;
    int place_grow = hit->place_grow;
    int place_layer = hit->place_layer;

    // Find chunk for placement position
    int cx = (int)floorf((float)place_gcol / HEX_CHUNK_SIZE);
    int cz = (int)floorf((float)place_grow / HEX_CHUNK_SIZE);
    int local_col = place_gcol - cx * HEX_CHUNK_SIZE;
    int local_row = place_grow - cz * HEX_CHUNK_SIZE;
    if (local_col < 0) { local_col += HEX_CHUNK_SIZE; cx--; }
    if (local_row < 0) { local_row += HEX_CHUNK_SIZE; cz--; }

    int chunk_idx = find_chunk(ht, cx, cz);
    if (chunk_idx < 0) return false;

    HexChunk* chunk = &ht->chunks[chunk_idx];
    if (!chunk->active || !chunk->voxels) return false;

    int local_layer = place_layer - chunk->base_layer;
    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) return false;
    if (HEX_VOXEL(chunk->voxels, local_col, local_row, local_layer) != VOXEL_AIR) return false;

    HEX_VOXEL(chunk->voxels, local_col, local_row, local_layer) = voxel_type;
    update_column_cache(chunk, local_col, local_row);
    chunk->dirty = true;
    chunk->edit_dirty = true;

    hex_log("[PLACE] type=%d at gcol=%d grow=%d layer=%d | chunk(%d,%d) local(%d,%d)\n",
            voxel_type, place_gcol, place_grow, place_layer, cx, cz, local_col, local_row);

    // Persist edit to disk-backed sector cache
    {
        double ux, uy, uz;
        gcol_grow_to_unit(ht, place_gcol, place_grow, &ux, &uy, &uz);
        edit_cache_record(&ht->edits, ux, uy, uz, place_layer, voxel_type);
    }

    // Always dirty all adjacent chunks: any block edit affects cross-chunk
    // BFS lighting (sky light range 32, torch range 15) and boundary geometry.
    int dirtied_count = 0;
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (!ht->chunks[i].active) continue;
        if (ht->chunks[i].cx == cx && ht->chunks[i].cz == cz) continue;
        int dx = ht->chunks[i].cx - cx;
        int dz = ht->chunks[i].cz - cz;
        if (dx >= -1 && dx <= 1 && dz >= -1 && dz <= 1) {
            ht->chunks[i].dirty = true;
            dirtied_count++;
        }
    }
    hex_log("[PLACE] dirtied %d neighbor chunks\n", dirtied_count);
    return true;
}

// ---- Ctrl placement (inverted ray) ----

// Ctrl mode: find placement face by inverting the view direction.
// Picks the side face on the FAR side of the selected block (away from camera).
// If that side neighbor is solid, falls back to top or bottom face.
void hex_terrain_ctrl_placement(const HexTerrain* ht, HexHitResult* hit,
                                 HMM_Vec3 camera_pos) {
    if (!hit->valid || !ht->frame_valid) return;

    // Get hex center and camera in tangent plane
    float hx, hz;
    hex_local_pos(hit->gcol, hit->grow, &hx, &hz);

    float cx, cz;
    world_to_tangent(ht, camera_pos, &cx, &cz);

    // Inverted direction: from camera TOWARD hex (opposite of camera-facing)
    float dx = hx - cx;
    float dz = hz - cz;

    // Find which side face best matches the inverted direction
    int best_dir = 0;
    float best_dot = -1e30f;
    for (int dir = 0; dir < 6; dir++) {
        int edge = DIR_TO_EDGE[dir];
        float nx = (HEX_COS[edge] + HEX_COS[(edge + 1) % 6]) * 0.5f;
        float nz = (HEX_SIN[edge] + HEX_SIN[(edge + 1) % 6]) * 0.5f;
        float d = dx * nx + dz * nz;
        if (d > best_dot) {
            best_dot = d;
            best_dir = dir;
        }
    }

    // Check if the neighbor on that side is solid
    int ngcol, ngrow;
    hex_neighbor(hit->gcol, hit->grow, best_dir, &ngcol, &ngrow);
    uint8_t neighbor = hex_terrain_get_voxel(ht, ngcol, ngrow, hit->layer);

    if (neighbor == VOXEL_AIR) {
        // Side is free — place there
        hit->face = 2;
        hit->place_gcol = ngcol;
        hit->place_grow = ngrow;
        hit->place_layer = hit->layer;
    } else {
        // Side blocked — fall back to top or bottom
        uint8_t above = hex_terrain_get_voxel(ht, hit->gcol, hit->grow, hit->layer + 1);
        if (above == VOXEL_AIR) {
            hit->face = 0;  // top
            hit->place_gcol = hit->gcol;
            hit->place_grow = hit->grow;
            hit->place_layer = hit->layer + 1;
        } else {
            hit->face = 1;  // bottom
            hit->place_gcol = hit->gcol;
            hit->place_grow = hit->grow;
            hit->place_layer = hit->layer - 1;
        }
    }
}

// ---- Highlight ----

// Full hex prism outline: top hex + bottom hex + 6 vertical edges = 18 line segments = 36 verts = 108 floats.
// Slightly oversized (HEX_RADIUS * 1.04) so the wireframe wraps visibly around the target block.
bool hex_terrain_build_highlight(const HexTerrain* ht, const HexHitResult* hit,
                                  const double world_origin[3], float* out_verts) {
    if (!hit->valid || !ht->frame_valid) return false;

    float lx, lz;
    hex_local_pos(hit->gcol, hit->grow, &lx, &lz);

    // Slightly oversized radius so outline wraps around the block
    float overshoot = HEX_RADIUS * 1.04f;
    float margin = 0.05f;  // radial outset so lines sit just outside block faces
    float bot_r = ht->planet_radius + (float)hit->layer * HEX_HEIGHT + HEX_SURFACE_BIAS - margin;
    float top_r = ht->planet_radius + (float)(hit->layer + 1) * HEX_HEIGHT + HEX_SURFACE_BIAS + margin;

    // Compute 6 vertex directions (oversized hex)
    HMM_Vec3 dirs[6];
    for (int i = 0; i < 6; i++) {
        float vx = lx + overshoot * HEX_COS[i];
        float vz = lz + overshoot * HEX_SIN[i];
        dirs[i] = vec3_normalize(
            vec3_add(ht->tangent_origin,
                vec3_add(vec3_scale(ht->tangent_east, vx),
                         vec3_scale(ht->tangent_north, vz))));
    }

    // Build 12 vertices for top hex and 12 for bottom hex (6 edges each)
    HMM_Vec3 top_v[6], bot_v[6];
    for (int i = 0; i < 6; i++) {
        top_v[i] = vec3_scale(dirs[i], top_r);
        bot_v[i] = vec3_scale(dirs[i], bot_r);
        // Apply floating origin offset
        top_v[i].X -= (float)world_origin[0];
        top_v[i].Y -= (float)world_origin[1];
        top_v[i].Z -= (float)world_origin[2];
        bot_v[i].X -= (float)world_origin[0];
        bot_v[i].Y -= (float)world_origin[1];
        bot_v[i].Z -= (float)world_origin[2];
    }

    int idx = 0;
    // Top hexagon (6 edges)
    for (int i = 0; i < 6; i++) {
        int j = (i + 1) % 6;
        out_verts[idx++] = top_v[i].X; out_verts[idx++] = top_v[i].Y; out_verts[idx++] = top_v[i].Z;
        out_verts[idx++] = top_v[j].X; out_verts[idx++] = top_v[j].Y; out_verts[idx++] = top_v[j].Z;
    }
    // Bottom hexagon (6 edges)
    for (int i = 0; i < 6; i++) {
        int j = (i + 1) % 6;
        out_verts[idx++] = bot_v[i].X; out_verts[idx++] = bot_v[i].Y; out_verts[idx++] = bot_v[i].Z;
        out_verts[idx++] = bot_v[j].X; out_verts[idx++] = bot_v[j].Y; out_verts[idx++] = bot_v[j].Z;
    }
    // 6 vertical edges connecting top to bottom
    for (int i = 0; i < 6; i++) {
        out_verts[idx++] = top_v[i].X; out_verts[idx++] = top_v[i].Y; out_verts[idx++] = top_v[i].Z;
        out_verts[idx++] = bot_v[i].X; out_verts[idx++] = bot_v[i].Y; out_verts[idx++] = bot_v[i].Z;
    }
    // idx should be 108 (36 vertices * 3 floats)

    return true;
}

int hex_terrain_build_placement_face(const HexTerrain* ht, const HexHitResult* hit,
                                      const double world_origin[3], float* out_verts) {
    if (!hit->valid || !ht->frame_valid) return 0;

    float lx, lz;
    hex_local_pos(hit->gcol, hit->grow, &lx, &lz);

    // Use the same offset as the white wireframe (slightly oversized to avoid z-fighting)
    float bot_r = ht->planet_radius + (float)hit->layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
    float top_r = bot_r + HEX_HEIGHT;

    // Inset factor: shrink vertices toward face center (0.85 = 15% inset)
    float inset = 0.85f;

    // Compute 6 vertex directions (at full hex radius)
    HMM_Vec3 dirs[6];
    for (int i = 0; i < 6; i++) {
        float vx = lx + HEX_RADIUS * HEX_COS[i];
        float vz = lz + HEX_RADIUS * HEX_SIN[i];
        dirs[i] = vec3_normalize(
            vec3_add(ht->tangent_origin,
                vec3_add(vec3_scale(ht->tangent_east, vx),
                         vec3_scale(ht->tangent_north, vz))));
    }

    // Center direction (for inset lerp)
    HMM_Vec3 center_dir = vec3_normalize(
        vec3_add(ht->tangent_origin,
            vec3_add(vec3_scale(ht->tangent_east, lx),
                     vec3_scale(ht->tangent_north, lz))));

    // Compute inset vertices: lerp each dir toward center
    HMM_Vec3 top_v[6], bot_v[6];
    for (int i = 0; i < 6; i++) {
        HMM_Vec3 inset_dir = vec3_normalize(vec3_lerp(center_dir, dirs[i], inset));
        top_v[i] = vec3_scale(inset_dir, top_r);
        bot_v[i] = vec3_scale(inset_dir, bot_r);
        top_v[i].X -= (float)world_origin[0];
        top_v[i].Y -= (float)world_origin[1];
        top_v[i].Z -= (float)world_origin[2];
        bot_v[i].X -= (float)world_origin[0];
        bot_v[i].Y -= (float)world_origin[1];
        bot_v[i].Z -= (float)world_origin[2];
    }

    int idx = 0;

    if (hit->face == 0) {
        // Top face: 6 line segments forming inset hexagon
        for (int i = 0; i < 6; i++) {
            int j = (i + 1) % 6;
            out_verts[idx++] = top_v[i].X; out_verts[idx++] = top_v[i].Y; out_verts[idx++] = top_v[i].Z;
            out_verts[idx++] = top_v[j].X; out_verts[idx++] = top_v[j].Y; out_verts[idx++] = top_v[j].Z;
        }
    } else if (hit->face == 1) {
        // Bottom face: 6 line segments forming inset hexagon
        for (int i = 0; i < 6; i++) {
            int j = (i + 1) % 6;
            out_verts[idx++] = bot_v[i].X; out_verts[idx++] = bot_v[i].Y; out_verts[idx++] = bot_v[i].Z;
            out_verts[idx++] = bot_v[j].X; out_verts[idx++] = bot_v[j].Y; out_verts[idx++] = bot_v[j].Z;
        }
    } else {
        // Side face: 4 line segments forming inset rectangle
        int edge = 0;
        if (hit->place_gcol != hit->gcol || hit->place_grow != hit->grow) {
            for (int d = 0; d < 6; d++) {
                int ncol, nrow;
                hex_neighbor(hit->gcol, hit->grow, d, &ncol, &nrow);
                if (ncol == hit->place_gcol && nrow == hit->place_grow) {
                    edge = DIR_TO_EDGE[d];
                    break;
                }
            }
        }
        int vi0 = edge;
        int vi1 = (edge + 1) % 6;

        // Inset the quad: also shrink vertically toward face center
        HMM_Vec3 bl = bot_v[vi0], br = bot_v[vi1];
        HMM_Vec3 tl = top_v[vi0], tr = top_v[vi1];

        // Vertical inset: lerp top/bottom toward the vertical midpoint
        HMM_Vec3 mid_l = vec3_lerp(bl, tl, 0.5f);
        HMM_Vec3 mid_r = vec3_lerp(br, tr, 0.5f);
        bl = vec3_lerp(mid_l, bl, inset);
        tl = vec3_lerp(mid_l, tl, inset);
        br = vec3_lerp(mid_r, br, inset);
        tr = vec3_lerp(mid_r, tr, inset);

        // Horizontal inset: lerp left/right toward horizontal midpoint
        HMM_Vec3 mid_b = vec3_lerp(bl, br, 0.5f);
        HMM_Vec3 mid_t = vec3_lerp(tl, tr, 0.5f);
        bl = vec3_lerp(mid_b, bl, inset);
        br = vec3_lerp(mid_b, br, inset);
        tl = vec3_lerp(mid_t, tl, inset);
        tr = vec3_lerp(mid_t, tr, inset);

        // 4 line segments: BL-BR, BR-TR, TR-TL, TL-BL
        out_verts[idx++] = bl.X; out_verts[idx++] = bl.Y; out_verts[idx++] = bl.Z;
        out_verts[idx++] = br.X; out_verts[idx++] = br.Y; out_verts[idx++] = br.Z;

        out_verts[idx++] = br.X; out_verts[idx++] = br.Y; out_verts[idx++] = br.Z;
        out_verts[idx++] = tr.X; out_verts[idx++] = tr.Y; out_verts[idx++] = tr.Z;

        out_verts[idx++] = tr.X; out_verts[idx++] = tr.Y; out_verts[idx++] = tr.Z;
        out_verts[idx++] = tl.X; out_verts[idx++] = tl.Y; out_verts[idx++] = tl.Z;

        out_verts[idx++] = tl.X; out_verts[idx++] = tl.Y; out_verts[idx++] = tl.Z;
        out_verts[idx++] = bl.X; out_verts[idx++] = bl.Y; out_verts[idx++] = bl.Z;
    }

    return idx / 3; // number of vertices
}

// ---- Voxel query API (for collision) ----

uint8_t hex_terrain_get_voxel(const HexTerrain* ht, int gcol, int grow, int world_layer) {
    int cx = (int)floorf((float)gcol / HEX_CHUNK_SIZE);
    int cz = (int)floorf((float)grow / HEX_CHUNK_SIZE);
    int local_col = gcol - cx * HEX_CHUNK_SIZE;
    int local_row = grow - cz * HEX_CHUNK_SIZE;
    if (local_col < 0) { local_col += HEX_CHUNK_SIZE; cx--; }
    if (local_row < 0) { local_row += HEX_CHUNK_SIZE; cz--; }

    int idx = find_chunk_const(ht, cx, cz);
    if (idx < 0) return VOXEL_AIR;

    const HexChunk* chunk = &ht->chunks[idx];
    if (!chunk->voxels) return VOXEL_AIR;

    int local_layer = world_layer - chunk->base_layer;
    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) return VOXEL_AIR;

    return HEX_VOXEL(chunk->voxels, local_col, local_row, local_layer);
}

float hex_terrain_ground_height(const HexTerrain* ht, int gcol, int grow) {
    int cx = (int)floorf((float)gcol / HEX_CHUNK_SIZE);
    int cz = (int)floorf((float)grow / HEX_CHUNK_SIZE);
    int local_col = gcol - cx * HEX_CHUNK_SIZE;
    int local_row = grow - cz * HEX_CHUNK_SIZE;
    if (local_col < 0) { local_col += HEX_CHUNK_SIZE; cx--; }
    if (local_row < 0) { local_row += HEX_CHUNK_SIZE; cz--; }

    int idx = find_chunk_const(ht, cx, cz);
    if (idx < 0) return 0.0f;

    const HexChunk* chunk = &ht->chunks[idx];
    if (!chunk->voxels) return 0.0f;

    // Scan from max_solid downward to find topmost solid with air above
    int max_s = chunk->col_max_solid[local_col][local_row];
    if (max_s < 0) return ht->planet_radius + chunk->base_layer * HEX_HEIGHT;

    for (int l = max_s; l >= 0; l--) {
        if (!is_transparent(HEX_VOXEL(chunk->voxels, local_col, local_row, l))) {
            uint8_t above = (l + 1 < HEX_CHUNK_LAYERS)
                ? HEX_VOXEL(chunk->voxels, local_col, local_row, l + 1)
                : VOXEL_AIR;
            if (is_transparent(above)) {
                int world_layer = l + chunk->base_layer;
                return ht->planet_radius + (world_layer + 1) * HEX_HEIGHT;
            }
        }
    }

    return ht->planet_radius + chunk->base_layer * HEX_HEIGHT;
}

float hex_terrain_ground_height_below(const HexTerrain* ht, int gcol, int grow,
                                       int max_world_layer) {
    int cx = (int)floorf((float)gcol / HEX_CHUNK_SIZE);
    int cz = (int)floorf((float)grow / HEX_CHUNK_SIZE);
    int local_col = gcol - cx * HEX_CHUNK_SIZE;
    int local_row = grow - cz * HEX_CHUNK_SIZE;
    if (local_col < 0) { local_col += HEX_CHUNK_SIZE; cx--; }
    if (local_row < 0) { local_row += HEX_CHUNK_SIZE; cz--; }

    int idx = find_chunk_const(ht, cx, cz);
    if (idx < 0) return 0.0f;

    const HexChunk* chunk = &ht->chunks[idx];
    if (!chunk->voxels) return 0.0f;

    // Convert max_world_layer to local layer, clamp to chunk range
    int start_local = max_world_layer - chunk->base_layer;
    if (start_local >= HEX_CHUNK_LAYERS) start_local = HEX_CHUNK_LAYERS - 1;
    if (start_local < 0) return ht->planet_radius + chunk->base_layer * HEX_HEIGHT;

    // Scan downward from start_local to find topmost solid with transparent above
    for (int l = start_local; l >= 0; l--) {
        if (!is_transparent(HEX_VOXEL(chunk->voxels, local_col, local_row, l))) {
            uint8_t above = (l + 1 < HEX_CHUNK_LAYERS)
                ? HEX_VOXEL(chunk->voxels, local_col, local_row, l + 1)
                : VOXEL_AIR;
            if (is_transparent(above)) {
                int world_layer = l + chunk->base_layer;
                return ht->planet_radius + (world_layer + 1) * HEX_HEIGHT;
            }
        }
    }

    return ht->planet_radius + chunk->base_layer * HEX_HEIGHT;
}

bool hex_terrain_has_headroom(const HexTerrain* ht, int gcol, int grow,
                               int feet_layer, int clearance) {
    for (int l = feet_layer; l < feet_layer + clearance; l++) {
        if (hex_terrain_get_voxel(ht, gcol, grow, l) != VOXEL_AIR) {
            return false;
        }
    }
    return true;
}

void hex_terrain_world_to_hex(const HexTerrain* ht, HMM_Vec3 world_pos,
                               int* out_gcol, int* out_grow, int* out_layer) {
    float lx, lz;
    world_to_tangent(ht, world_pos, &lx, &lz);
    pixel_to_hex(lx, lz, out_gcol, out_grow);

    float point_r = sqrtf(vec3_dot(world_pos, world_pos));
    float altitude = point_r - ht->planet_radius;
    *out_layer = (int)floorf(altitude / HEX_HEIGHT);
}

void hex_terrain_hex_to_world(const HexTerrain* ht, int gcol, int grow, int layer,
                               float* out_x, float* out_y, float* out_z,
                               float* out_up_x, float* out_up_y, float* out_up_z) {
    if (!ht->frame_valid) {
        *out_x = *out_y = *out_z = 0.0f;
        *out_up_x = 0.0f; *out_up_y = 1.0f; *out_up_z = 0.0f;
        return;
    }
    float lx, lz;
    hex_local_pos(gcol, grow, &lx, &lz);
    float radius = ht->planet_radius + layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
    HMM_Vec3 pos = tangent_to_world(ht, lx, lz, radius);
    *out_x = pos.X;
    *out_y = pos.Y;
    *out_z = pos.Z;
    HMM_Vec3 up = vec3_normalize(pos);
    *out_up_x = up.X;
    *out_up_y = up.Y;
    *out_up_z = up.Z;
}

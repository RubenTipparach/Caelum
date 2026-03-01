#include "hex_terrain.h"
#include "hex_vertex.h"
#include "math_utils.h"
#include "planet.h"
#include "log_config.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FastNoiseLite.h"

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

// ---- Terrain noise (identical to lod.c) ----
#define TERRAIN_SEA_LEVEL_M   4000.0f
#define TERRAIN_AMPLITUDE_M   8000.0f
#define TERRAIN_MIN_M         500.0f

static float ht_smoothstepf(float edge0, float edge1, float x) {
    float t = (x - edge0) / (edge1 - edge0);
    t = fminf(1.0f, fmaxf(0.0f, t));
    return t * t * (3.0f - 2.0f * t);
}

static fnl_state ht_create_continental_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 3;
    noise.frequency = 0.6f;
    noise.seed = seed;
    return noise;
}

static fnl_state ht_create_mountain_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_RIDGED;
    noise.octaves = 5;
    noise.frequency = 1.5f;
    noise.seed = seed + 4000;
    return noise;
}

static fnl_state ht_create_warp_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 3;
    noise.frequency = 4.0f;
    noise.seed = seed + 1000;
    return noise;
}

static fnl_state ht_create_detail_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_RIDGED;
    noise.octaves = 3;
    noise.frequency = 16.0f;
    noise.seed = seed + 2000;
    return noise;
}

static float ht_sample_terrain_noise(fnl_state* continental, fnl_state* mountain,
                                      fnl_state* warp_noise, fnl_state* detail_noise,
                                      HMM_Vec3 unit_pos) {
    float scale = 3.0f;
    float px = unit_pos.X * scale;
    float py = unit_pos.Y * scale;
    float pz = unit_pos.Z * scale;

    float warp_strength = 0.5f;
    float wx = fnlGetNoise3D(warp_noise, px + 5.2f, py + 1.3f, pz + 3.7f);
    float wy = fnlGetNoise3D(warp_noise, px + 9.1f, py + 4.8f, pz + 7.2f);
    float wz = fnlGetNoise3D(warp_noise, px + 2.6f, py + 8.4f, pz + 0.9f);
    float wpx = px + wx * warp_strength;
    float wpy = py + wy * warp_strength;
    float wpz = pz + wz * warp_strength;

    float continent = fnlGetNoise3D(continental, wpx, wpy, wpz);

    float mountain_raw = fnlGetNoise3D(mountain, wpx, wpy, wpz);
    float mountain_val = (mountain_raw + 1.0f) * 0.5f;
    mountain_val *= mountain_val;
    float land_factor = ht_smoothstepf(-0.05f, 0.35f, continent);
    float mountain_height = mountain_val * land_factor;

    float detail = fnlGetNoise3D(detail_noise, px, py, pz);
    float detail_weight = 0.05f + land_factor * 0.10f;

    float height = continent * 0.55f - 0.22f;
    height += mountain_height * 0.45f;
    height += detail * detail_weight;

    if (height > 1.0f) height = 1.0f;
    if (height < -1.0f) height = -1.0f;
    return height;
}

static float ht_sample_height_m(fnl_state* continental, fnl_state* mountain,
                                  fnl_state* warp_noise, fnl_state* detail_noise,
                                  HMM_Vec3 unit_pos) {
    float n = ht_sample_terrain_noise(continental, mountain, warp_noise, detail_noise, unit_pos);
    float height = TERRAIN_MIN_M + (n + 1.0f) * 0.5f * TERRAIN_AMPLITUDE_M;
    if (height < 0.0f) height = 0.0f;
    return height;
}

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
    f.center = (HMM_Vec3){{0.0f, 1.0f, 0.0f}};
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

// Voxel type by height (matches lod.c compute_voxel_type logic)
static VoxelType ht_voxel_type(float height_m) {
    float rel = height_m - TERRAIN_SEA_LEVEL_M;
    if (height_m < TERRAIN_SEA_LEVEL_M) return VOXEL_WATER;
    if (rel < 200.0f) return VOXEL_SAND;
    if (rel < 1500.0f) return VOXEL_GRASS;
    if (rel < 3000.0f) return VOXEL_STONE;
    if (rel < 4500.0f) return VOXEL_STONE;
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
#define HEX_ATLAS_TILES 9
#define ATLAS_WATER      0
#define ATLAS_SAND       1
#define ATLAS_DIRT       2
#define ATLAS_GRASS      3
#define ATLAS_STONE      4
#define ATLAS_ICE        5
#define ATLAS_SNOW       6
#define ATLAS_DIRT_GRASS  7
#define ATLAS_DIRT_SNOW   8

static int voxel_top_atlas(VoxelType type) {
    switch (type) {
        case VOXEL_WATER: return ATLAS_WATER;
        case VOXEL_SAND:  return ATLAS_SAND;
        case VOXEL_DIRT:  return ATLAS_DIRT;
        case VOXEL_GRASS: return ATLAS_GRASS;
        case VOXEL_STONE: return ATLAS_STONE;
        case VOXEL_ICE:   return ATLAS_ICE;
        default:          return ATLAS_GRASS;
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

    int chunk_index;
    volatile int completed;
} HexMeshJob;

// ---- Emit triangle helper (textured) ----

typedef struct { float u, v; } HexUV;

static void hex_emit_tri(HexVertex** verts, int* count, int* cap,
                          HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                          HMM_Vec3 winding_normal,
                          HMM_Vec3 shading_normal,
                          HexUV uv0, HexUV uv1, HexUV uv2,
                          HMM_Vec3 vtx_color) {
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

// Get neighbor voxel, handling cross-chunk boundaries by re-sampling noise
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

    // Cross-chunk boundary: re-sample terrain from noise
    float nlx, nlz;
    hex_local_pos(ncol_g, nrow_g, &nlx, &nlz);
    HMM_Vec3 ndir = vec3_normalize(
        vec3_add(job->tangent_origin,
            vec3_add(vec3_scale(job->tangent_east, nlx),
                     vec3_scale(job->tangent_north, nlz))));
    float nh_m = ht_sample_height_m(continental, mountain, warp, detail, ndir);
    float nh_eff = fmaxf(nh_m, TERRAIN_SEA_LEVEL_M);
    int surface_layer = (int)ceilf(nh_eff / HEX_HEIGHT);

    int world_layer = layer + job->base_layer;

    // Check arch at this neighbor too
    float alx, alz;
    ht_arch_project(&job->arch_frame, vec3_normalize(ndir), job->planet_radius, &alx, &alz);
    int ab, at;
    ht_arch_query(alx, alz, job->arch_base, &ab, &at);
    if (ab != INT16_MIN && world_layer >= ab && world_layer <= at)
        return VOXEL_STONE;

    if (world_layer <= surface_layer) return VOXEL_STONE; // approximate as solid
    return VOXEL_AIR;
}

// ---- Core mesh generation for a chunk ----

static void generate_chunk_mesh(HexMeshJob* job) {
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
                VoxelType surface_type = ht_voxel_type(h_m);

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
    }

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

            // Compute per-column terrain color for vertex color fade
            float col_h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, center_dir);
            HMM_Vec3 col_color = ht_terrain_color(col_h_m);

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
                        local_up, local_up, top_uvs[0], top_uvs[i + 1], top_uvs[i + 2],
                        col_color);
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
                if (vtype == VOXEL_AIR) continue;

                int world_layer = l + job->base_layer;
                float bot_r = job->planet_radius + world_layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
                float top_r = bot_r + HEX_HEIGHT;

                // --- TOP FACE: emit if layer above is AIR ---
                uint8_t above = job_get_voxel(job, col, row, l + 1);
                if (above == VOXEL_AIR) {
                    int atlas = voxel_top_atlas((VoxelType)vtype);
                    HMM_Vec3 verts[6];
                    HexUV uvs[6];
                    for (int i = 0; i < 6; i++) {
                        verts[i] = vec3_scale(hex_dirs[i], top_r);
                        uvs[i] = hex_top_uv(hex_vx[i], hex_vz[i], lx, lz, atlas);
                    }
                    for (int i = 0; i < 4; i++) {
                        hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                            verts[0], verts[i + 1], verts[i + 2],
                            local_up, local_up, uvs[0], uvs[i + 1], uvs[i + 2],
                            col_color);
                    }
                }

                // --- BOTTOM FACE: emit if layer below is AIR ---
                uint8_t below = job_get_voxel(job, col, row, l - 1);
                if (below == VOXEL_AIR) {
                    int atlas = voxel_bottom_atlas((VoxelType)vtype);
                    HMM_Vec3 verts[6];
                    HexUV uvs[6];
                    for (int i = 0; i < 6; i++) {
                        verts[i] = vec3_scale(hex_dirs[i], bot_r);
                        uvs[i] = hex_top_uv(hex_vx[i], hex_vz[i], lx, lz, atlas);
                    }
                    // Reversed winding for bottom face
                    for (int i = 0; i < 4; i++) {
                        hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                            verts[0], verts[i + 2], verts[i + 1],
                            local_down, local_down, uvs[0], uvs[i + 2], uvs[i + 1],
                            col_color);
                    }
                }

                // --- SIDE FACES: 6 hex neighbors ---
                for (int dir = 0; dir < 6; dir++) {
                    uint8_t nv = job_get_neighbor_voxel(job, &continental, &mountain, &warp, &detail,
                                                         col, row, l, dir);
                    if (nv != VOXEL_AIR) continue;

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
                        col_color);
                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        p0b, p1t, p0t, wall_outs[dir], local_up, wuv00, wuv11, wuv01,
                        col_color);
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

    if (chunk->voxels) {
        free(chunk->voxels);
        chunk->voxels = NULL;
    }

    if (chunk->pending_job) {
        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        if (job->completed) {
            if (job->vertices) free(job->vertices);
            if (job->voxels) free(job->voxels);
            free(job);
        } else {
            orphan_job(job);
        }
        chunk->pending_job = NULL;
    }

    chunk->active = false;
    chunk->dirty = false;
    chunk->generating = false;
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
        free(job);
    }
    s_orphan_count = 0;
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
    }

    // Process completed mesh generation jobs
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (!chunk->active || !chunk->generating || !chunk->pending_job) continue;

        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        if (!job->completed) continue;

        // Transfer results: mesh
        chunk->cpu_vertices = job->vertices;
        chunk->cpu_vertex_count = job->vertex_count;
        job->vertices = NULL;

        // Transfer results: voxel data
        if (chunk->voxels) free(chunk->voxels);
        chunk->voxels = job->voxels;
        job->voxels = NULL;
        chunk->base_layer = job->base_layer;
        memcpy(chunk->col_min_solid, job->col_min_solid, sizeof(chunk->col_min_solid));
        memcpy(chunk->col_max_solid, job->col_max_solid, sizeof(chunk->col_max_solid));

        free(job);
        chunk->pending_job = NULL;
        chunk->generating = false;
        chunk->dirty = false;
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
            bool want_transition = (dist > HEX_INNER_RANGE);

            if (idx >= 0) {
                HexChunk* chunk = &ht->chunks[idx];
                if (chunk->gpu_vertex_count == 0 && !chunk->generating &&
                    chunk->cpu_vertex_count == 0) {
                    chunk->dirty = true;
                }
                if (!chunk->is_transition && dist > HEX_TRANSITION_ON) {
                    chunk->is_transition = true;
                    if (!chunk->generating) chunk->dirty = true;
                } else if (chunk->is_transition && dist < HEX_TRANSITION_OFF) {
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
                chunk->is_transition = want_transition;
                chunk->gpu_buffer.id = SG_INVALID_ID;
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
                        // Will be filled by generate_chunk_mesh
                    }

                    job->vertices = NULL;
                    job->vertex_count = 0;
                    job->completed = 0;

                    if (job_system_try_submit(ht->jobs, hex_mesh_worker, job)) {
                        chunk->pending_job = job;
                        chunk->generating = true;
                        chunk->dirty = false;
                        jobs_submitted++;
                    } else {
                        free(job->voxels);
                        free(job);
                    }
                }
            }
        }
        if (new_activations >= HEX_MAX_ACTIVATIONS || jobs_submitted >= 8) break;
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
    int meshed = 0, total = 0;
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (!ht->chunks[i].active) continue;
        if (ht->chunks[i].is_transition) continue;
        total++;
        if (ht->chunks[i].gpu_vertex_count > 0) meshed++;
    }
    if (total == 0) return 0.0f;
    float coverage = (float)meshed / (float)total;
    return (coverage >= 0.9f) ? HEX_INNER_RANGE : 0.0f;
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

        // Project to tangent plane
        HMM_Vec3 offset = vec3_sub(world_pt, ht->tangent_origin);
        float lx = vec3_dot(offset, ht->tangent_east);
        float lz = vec3_dot(offset, ht->tangent_north);

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
            // HIT
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

            // Determine hit face and placement position
            if (world_layer != prev_layer && gcol == prev_gcol && grow == prev_grow) {
                // Vertical transition
                if (world_layer > prev_layer) {
                    result.face = 1; // bottom face
                    result.place_gcol = gcol;
                    result.place_grow = grow;
                    result.place_layer = world_layer - 1;
                } else {
                    result.face = 0; // top face
                    result.place_gcol = gcol;
                    result.place_grow = grow;
                    result.place_layer = world_layer + 1;
                }
            } else if (gcol != prev_gcol || grow != prev_grow) {
                // Horizontal transition: side face
                result.face = 2; // generic side
                result.place_gcol = prev_gcol;
                result.place_grow = prev_grow;
                result.place_layer = world_layer;
            } else {
                // First step hit (no previous)
                result.face = 0;
                result.place_gcol = gcol;
                result.place_grow = grow;
                result.place_layer = world_layer + 1;
            }
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
    if (!hit->valid || hit->chunk_index < 0 || hit->chunk_index >= HEX_MAX_CHUNKS) return false;

    HexChunk* chunk = &ht->chunks[hit->chunk_index];
    if (!chunk->active || !chunk->voxels) return false;

    int local_layer = hit->layer - chunk->base_layer;
    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) return false;

    HEX_VOXEL(chunk->voxels, hit->col, hit->row, local_layer) = VOXEL_AIR;
    update_column_cache(chunk, hit->col, hit->row);
    chunk->dirty = true;

    // Dirty neighbor chunks if at chunk boundary
    if (hit->col == 0 || hit->col == HEX_CHUNK_SIZE - 1 ||
        hit->row == 0 || hit->row == HEX_CHUNK_SIZE - 1) {
        for (int dir = 0; dir < 6; dir++) {
            int ncol_g, nrow_g;
            hex_neighbor(hit->gcol, hit->grow, dir, &ncol_g, &nrow_g);
            int ncx = (int)floorf((float)ncol_g / HEX_CHUNK_SIZE);
            int ncz = (int)floorf((float)nrow_g / HEX_CHUNK_SIZE);
            if (ncx != chunk->cx || ncz != chunk->cz) {
                int ni = find_chunk(ht, ncx, ncz);
                if (ni >= 0 && !ht->chunks[ni].generating) {
                    ht->chunks[ni].dirty = true;
                }
            }
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

    // Dirty neighbor chunks if at boundary
    if (local_col == 0 || local_col == HEX_CHUNK_SIZE - 1 ||
        local_row == 0 || local_row == HEX_CHUNK_SIZE - 1) {
        for (int dir = 0; dir < 6; dir++) {
            int ncol_g, nrow_g;
            hex_neighbor(place_gcol, place_grow, dir, &ncol_g, &nrow_g);
            int ncx = (int)floorf((float)ncol_g / HEX_CHUNK_SIZE);
            int ncz = (int)floorf((float)nrow_g / HEX_CHUNK_SIZE);
            if (ncx != chunk->cx || ncz != chunk->cz) {
                int ni = find_chunk(ht, ncx, ncz);
                if (ni >= 0 && !ht->chunks[ni].generating) {
                    ht->chunks[ni].dirty = true;
                }
            }
        }
    }
    return true;
}

// ---- Highlight ----

bool hex_terrain_build_highlight(const HexTerrain* ht, const HexHitResult* hit,
                                  const double world_origin[3], float* out_verts) {
    if (!hit->valid || !ht->frame_valid) return false;

    float lx, lz;
    hex_local_pos(hit->gcol, hit->grow, &lx, &lz);

    // Highlight at the specific layer, not column top
    float top_r = ht->planet_radius + (float)(hit->layer + 1) * HEX_HEIGHT + HEX_SURFACE_BIAS + 0.05f;

    HMM_Vec3 hex_v[6];
    for (int i = 0; i < 6; i++) {
        float vx = lx + HEX_RADIUS * HEX_COS[i];
        float vz = lz + HEX_RADIUS * HEX_SIN[i];
        HMM_Vec3 vdir = vec3_normalize(
            vec3_add(ht->tangent_origin,
                vec3_add(vec3_scale(ht->tangent_east, vx),
                         vec3_scale(ht->tangent_north, vz))));
        hex_v[i] = vec3_scale(vdir, top_r);
    }

    for (int i = 0; i < 6; i++) {
        hex_v[i].X -= (float)world_origin[0];
        hex_v[i].Y -= (float)world_origin[1];
        hex_v[i].Z -= (float)world_origin[2];
    }

    for (int i = 0; i < 6; i++) {
        int j = (i + 1) % 6;
        out_verts[i * 6 + 0] = hex_v[i].X;
        out_verts[i * 6 + 1] = hex_v[i].Y;
        out_verts[i * 6 + 2] = hex_v[i].Z;
        out_verts[i * 6 + 3] = hex_v[j].X;
        out_verts[i * 6 + 4] = hex_v[j].Y;
        out_verts[i * 6 + 5] = hex_v[j].Z;
    }

    return true;
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
        if (HEX_VOXEL(chunk->voxels, local_col, local_row, l) != VOXEL_AIR) {
            uint8_t above = (l + 1 < HEX_CHUNK_LAYERS)
                ? HEX_VOXEL(chunk->voxels, local_col, local_row, l + 1)
                : VOXEL_AIR;
            if (above == VOXEL_AIR) {
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
    HMM_Vec3 offset = vec3_sub(world_pos, ht->tangent_origin);
    float lx = vec3_dot(offset, ht->tangent_east);
    float lz = vec3_dot(offset, ht->tangent_north);
    pixel_to_hex(lx, lz, out_gcol, out_grow);

    float point_r = sqrtf(vec3_dot(world_pos, world_pos));
    float altitude = point_r - ht->planet_radius;
    *out_layer = (int)floorf(altitude / HEX_HEIGHT);
}

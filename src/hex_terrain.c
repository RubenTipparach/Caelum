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
        case VOXEL_WATER:   return ATLAS_WATER;
        case VOXEL_SAND:    return ATLAS_SAND;
        case VOXEL_DIRT:    return ATLAS_DIRT;
        case VOXEL_GRASS:   return ATLAS_GRASS;
        case VOXEL_STONE:   return ATLAS_STONE;
        case VOXEL_ICE:     return ATLAS_ICE;
        case VOXEL_BEDROCK: return ATLAS_STONE;
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
    uint8_t* sky_map; // Per-voxel sky light (0..SKY_MAX), same layout as voxels (malloc'd)
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

    // Edit cache (read-only pointer for edit queries on worker thread)
    const EditCache* edit_cache;
} HexMeshJob;

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
                          HMM_Vec3 vtx_color, float sky_light) {
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

    // Be conservative near any surface boundary to prevent missing faces.
    // The noise approximation can disagree with actual voxel data by ±1 layer.
    // Returning VOXEL_AIR forces face emission (safe: worst case = hidden extra face).
    int margin = 2;

    // Near arch surfaces: conservative
    if (ab != INT16_MIN) {
        if (world_layer >= ab - margin && world_layer <= at + margin) {
            // Near arch boundary — could be inside or outside the shell
            if (world_layer >= ab && world_layer <= at)
                return VOXEL_STONE;  // definitely inside arch shell
            return VOXEL_AIR;  // near arch edge, be safe
        }
    }

    // Near terrain surface: conservative
    if (world_layer >= surface_layer - margin && world_layer <= surface_layer + margin)
        return VOXEL_AIR;

    // Well below surface: definitely solid
    if (world_layer < surface_layer - margin) return VOXEL_STONE;

    // Well above surface: definitely air
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
            // Scan from the TOP of the chunk downward until we hit a solid block
            for (int l = HEX_CHUNK_LAYERS - 1; l >= 0; l--) {
                if (HEX_VOXEL(job->voxels, col, row, l) != VOXEL_AIR) break;
                SKY_VOXEL(sky, col, row, l) = SKY_MAX;
            }
        }
    }

    // ---- Pass 2: Seed BFS queue ----
    // Any sky-lit air block adjacent to a non-sky-lit air block is a seed.
    // Also seed from sky-lit blocks adjacent to columns with higher terrain
    // (cross-chunk boundary awareness).
    int* queue = (int*)malloc(SKY_BFS_CAP * sizeof(int));
    int qhead = 0, qtail = 0;

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
                            HEX_VOXEL(job->voxels, nc, nr, l) == VOXEL_AIR &&
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
            if (HEX_VOXEL(job->voxels, c, r, nl) != VOXEL_AIR) continue;

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

        // 6 lateral hex neighbors (in-chunk only)
        int gcol = job->cx * HEX_CHUNK_SIZE + c;
        int grow = job->cz * HEX_CHUNK_SIZE + r;
        for (int dir = 0; dir < 6; dir++) {
            int ncol_g, nrow_g;
            hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);
            int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
            int nr_local = nrow_g - job->cz * HEX_CHUNK_SIZE;

            if (nc < 0 || nc >= HEX_CHUNK_SIZE || nr_local < 0 || nr_local >= HEX_CHUNK_SIZE)
                continue;
            if (HEX_VOXEL(job->voxels, nc, nr_local, l) != VOXEL_AIR) continue;

            uint8_t new_light = current - 1;
            if (new_light > SKY_VOXEL(sky, nc, nr_local, l)) {
                SKY_VOXEL(sky, nc, nr_local, l) = new_light;
                if (qtail < SKY_BFS_CAP) {
                    queue[qtail++] = nc * HEX_CHUNK_SIZE * HEX_CHUNK_LAYERS
                                   + nr_local * HEX_CHUNK_LAYERS + l;
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
                        col_color, 1.0f);
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
                    for (int i = 0; i < 4; i++) {
                        hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                            verts[0], verts[i + 1], verts[i + 2],
                            local_up, slope_normal, uvs[0], uvs[i + 1], uvs[i + 2],
                            col_color, face_sky);
                    }
                }

                // --- BOTTOM FACE: emit if layer below is AIR ---
                uint8_t below = job_get_voxel(job, col, row, l - 1);
                if (below == VOXEL_AIR) {
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
                    for (int i = 0; i < 4; i++) {
                        hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                            verts[0], verts[i + 2], verts[i + 1],
                            local_down, local_down, uvs[0], uvs[i + 2], uvs[i + 1],
                            col_color, face_sky);
                    }
                }

                // --- SIDE FACES: 6 hex neighbors ---
                for (int dir = 0; dir < 6; dir++) {
                    uint8_t nv = job_get_neighbor_voxel(job, &continental, &mountain, &warp, &detail,
                                                         col, row, l, dir);
                    if (nv != VOXEL_AIR) continue;

                    // Column sky light from the air block the side faces into
                    float face_sky = 1.0f;  // default for cross-chunk
                    {
                        int ncol_g, nrow_g;
                        hex_neighbor(job->cx * HEX_CHUNK_SIZE + col,
                                     job->cz * HEX_CHUNK_SIZE + row, dir, &ncol_g, &nrow_g);
                        int nc = ncol_g - job->cx * HEX_CHUNK_SIZE;
                        int nr_local = nrow_g - job->cz * HEX_CHUNK_SIZE;
                        if (nc >= 0 && nc < HEX_CHUNK_SIZE &&
                            nr_local >= 0 && nr_local < HEX_CHUNK_SIZE) {
                            face_sky = read_sky_light(job, nc, nr_local, l);
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
                        col_color, face_sky);
                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        p0b, p1t, p0t, wall_outs[dir], local_up, wuv00, wuv11, wuv01,
                        col_color, face_sky);
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
                if (vtype == VOXEL_AIR) continue;

                int world_layer = l + job->base_layer;
                float bot_r = job->planet_radius + world_layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
                float top_r = bot_r + HEX_HEIGHT;

                // Top face outline: hexagon at top_r
                uint8_t above = (l + 1 < HEX_CHUNK_LAYERS)
                    ? HEX_VOXEL(job->voxels, col, row, l + 1) : VOXEL_AIR;
                if (above == VOXEL_AIR) {
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
                if (below == VOXEL_AIR) {
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
                        exposed = (HEX_VOXEL(job->voxels, nc, nr, l) == VOXEL_AIR);
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

    // Free sky_map (only needed during mesh gen, not after)
    if (job->sky_map) {
        free(job->sky_map);
        job->sky_map = NULL;
    }

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
            if (job->wire_verts) free(job->wire_verts);
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

    if (chunk->pending_job) {
        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        if (job->completed) {
            if (job->vertices) free(job->vertices);
            if (job->voxels) free(job->voxels);
            if (job->sky_map) free(job->sky_map);
            if (job->wire_verts) free(job->wire_verts);
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
        if (job->wire_verts) free(job->wire_verts);
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
        }

        free(job);
        chunk->pending_job = NULL;
        chunk->generating = false;
        // Don't clear dirty — it was already cleared at submission time.
        // If it's true now, something re-dirtied it and we need another pass.
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
                        // Will be filled by generate_chunk_mesh
                    }

                    job->vertices = NULL;
                    job->vertex_count = 0;
                    job->completed = 0;

                    if (job_system_try_submit(ht->jobs, hex_mesh_worker, job)) {
                        chunk->pending_job = job;
                        chunk->generating = true;
                        chunk->dirty = false;
                        chunk->edit_dirty = false;
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
            }

            job->vertices = NULL;
            job->vertex_count = 0;
            job->completed = 0;

            if (job_system_try_submit(ht->jobs, hex_mesh_worker, job)) {
                chunk->pending_job = job;
                chunk->generating = true;
                chunk->dirty = false;
                chunk->edit_dirty = false;
                jobs_submitted++;
            } else {
                free(job->voxels);
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
    if (!hit->valid || hit->chunk_index < 0 || hit->chunk_index >= HEX_MAX_CHUNKS) return false;

    HexChunk* chunk = &ht->chunks[hit->chunk_index];
    if (!chunk->active || !chunk->voxels) return false;

    int local_layer = hit->layer - chunk->base_layer;
    if (local_layer < 0 || local_layer >= HEX_CHUNK_LAYERS) return false;

    // Bedrock is unbreakable
    if (HEX_VOXEL(chunk->voxels, hit->col, hit->row, local_layer) == VOXEL_BEDROCK) return false;

    HEX_VOXEL(chunk->voxels, hit->col, hit->row, local_layer) = VOXEL_AIR;
    update_column_cache(chunk, hit->col, hit->row);
    chunk->dirty = true;
    chunk->edit_dirty = true;

    // Persist edit to disk-backed sector cache
    {
        double ux, uy, uz;
        gcol_grow_to_unit(ht, hit->gcol, hit->grow, &ux, &uy, &uz);
        edit_cache_record(&ht->edits, ux, uy, uz, hit->layer, VOXEL_AIR);
    }

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
                if (ni >= 0) {
                    ht->chunks[ni].dirty = true;
                    ht->chunks[ni].edit_dirty = true;
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
    chunk->edit_dirty = true;

    // Persist edit to disk-backed sector cache
    {
        double ux, uy, uz;
        gcol_grow_to_unit(ht, place_gcol, place_grow, &ux, &uy, &uz);
        edit_cache_record(&ht->edits, ux, uy, uz, place_layer, voxel_type);
    }

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
                if (ni >= 0) {
                    ht->chunks[ni].dirty = true;
                    ht->chunks[ni].edit_dirty = true;
                }
            }
        }
    }
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

    // Scan downward from start_local to find topmost solid with air above
    for (int l = start_local; l >= 0; l--) {
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
    float lx, lz;
    world_to_tangent(ht, world_pos, &lx, &lz);
    pixel_to_hex(lx, lz, out_gcol, out_grow);

    float point_r = sqrtf(vec3_dot(world_pos, world_pos));
    float altitude = point_r - ht->planet_radius;
    *out_layer = (int)floorf(altitude / HEX_HEIGHT);
}

#include "hex_terrain.h"
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
// Direction 0: NE (+1, 0), Direction 1: SE (+1,-1), Direction 2: S (0,-1)
// Direction 3: SW (-1, 0), Direction 4: NW (-1,+1), Direction 5: N (0,+1)
static const int AX_NEIGHBORS[6][2] = {
    { +1,  0 },   // direction 0: NE
    { +1, -1 },   // direction 1: SE
    {  0, -1 },   // direction 2: S
    { -1,  0 },   // direction 3: SW
    { -1, +1 },   // direction 4: NW
    {  0, +1 },   // direction 5: N
};

// Map neighbor direction to the hex edge index that faces that neighbor.
// Edge i = the edge between vertex[i] and vertex[(i+1)%6].
// Flat-topped hex vertices: v0=E, v1=NE, v2=NW, v3=W, v4=SW, v5=SE
// Edge 0 (v0-v1) faces NE, Edge 1 (v1-v2) faces N, Edge 2 (v2-v3) faces NW
// Edge 3 (v3-v4) faces SW, Edge 4 (v4-v5) faces S, Edge 5 (v5-v0) faces SE
static const int DIR_TO_EDGE[6] = {
    0,  // dir 0 (NE neighbor) → edge 0 (v0-v1, faces NE)
    5,  // dir 1 (SE neighbor) → edge 5 (v5-v0, faces SE)
    4,  // dir 2 (S  neighbor) → edge 4 (v4-v5, faces S)
    3,  // dir 3 (SW neighbor) → edge 3 (v3-v4, faces SW)
    2,  // dir 4 (NW neighbor) → edge 2 (v2-v3, faces NW)
    1,  // dir 5 (N  neighbor) → edge 1 (v1-v2, faces N)
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

// ---- Terrain noise (same as lod.c — shared constants) ----
#define TERRAIN_SEA_LEVEL_M   4000.0f
#define TERRAIN_AMPLITUDE_M   8000.0f
#define TERRAIN_MIN_M         500.0f

static fnl_state ht_create_terrain_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 4;
    noise.frequency = 2.0f;
    noise.seed = seed;
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

static float ht_sample_terrain_noise(fnl_state* base, fnl_state* warp,
                                      fnl_state* detail, HMM_Vec3 unit_pos) {
    float scale = 3.0f;
    float px = unit_pos.X * scale;
    float py = unit_pos.Y * scale;
    float pz = unit_pos.Z * scale;

    float warp_strength = 0.3f;
    float wx = fnlGetNoise3D(warp, px + 5.2f, py + 1.3f, pz + 3.7f);
    float wy = fnlGetNoise3D(warp, px + 9.1f, py + 4.8f, pz + 7.2f);
    float wz = fnlGetNoise3D(warp, px + 2.6f, py + 8.4f, pz + 0.9f);
    px += wx * warp_strength;
    py += wy * warp_strength;
    pz += wz * warp_strength;

    float n = fnlGetNoise3D(base, px, py, pz);
    float d = fnlGetNoise3D(detail,
        unit_pos.X * scale, unit_pos.Y * scale, unit_pos.Z * scale);
    float detail_weight = fmaxf(0.0f, n * 0.5f + 0.5f) * 0.3f;
    n += d * detail_weight;
    return n;
}

static float ht_sample_height_m(fnl_state* base, fnl_state* warp,
                                 fnl_state* detail, HMM_Vec3 unit_pos) {
    float n = ht_sample_terrain_noise(base, warp, detail, unit_pos);
    float height = TERRAIN_MIN_M + (n + 1.0f) * 0.5f * TERRAIN_AMPLITUDE_M;
    if (height < 0.0f) height = 0.0f;
    return height;
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

static HMM_Vec3 ht_voxel_color(VoxelType type) {
    switch (type) {
        case VOXEL_WATER: return (HMM_Vec3){{0.15f, 0.35f, 0.65f}};
        case VOXEL_SAND:  return (HMM_Vec3){{0.76f, 0.70f, 0.40f}};
        case VOXEL_DIRT:  return (HMM_Vec3){{0.45f, 0.32f, 0.18f}};
        case VOXEL_GRASS: return (HMM_Vec3){{0.20f, 0.60f, 0.15f}};
        case VOXEL_STONE: return (HMM_Vec3){{0.50f, 0.50f, 0.50f}};
        case VOXEL_ICE:   return (HMM_Vec3){{0.90f, 0.95f, 1.00f}};
        default:          return (HMM_Vec3){{1.0f, 0.0f, 1.0f}};
    }
}

// Color for exposed side walls (dirt/stone below surface)
static HMM_Vec3 ht_wall_color(float height_m, float wall_y) {
    float depth_below_surface = height_m - wall_y;
    if (depth_below_surface < 1.0f) {
        return ht_voxel_color(ht_voxel_type(height_m));
    } else if (depth_below_surface < 3.0f) {
        return ht_voxel_color(VOXEL_DIRT);
    } else {
        return ht_voxel_color(VOXEL_STONE);
    }
}

// ---- Tangent frame computation ----

static void compute_tangent_frame(HexTerrain* ht, HMM_Vec3 camera_pos) {
    // Camera radial direction = local up
    float cam_r = sqrtf(vec3_dot(camera_pos, camera_pos));
    if (cam_r < 1.0f) cam_r = 1.0f;
    HMM_Vec3 cam_dir = vec3_scale(camera_pos, 1.0f / cam_r);

    ht->tangent_up = cam_dir;

    // Ground point: project camera down to terrain surface
    // (Use planet_radius + sea_level as approximate ground)
    float ground_r = ht->planet_radius + TERRAIN_SEA_LEVEL_M;
    ht->tangent_origin = vec3_scale(cam_dir, ground_r);

    // Build tangent frame: choose an arbitrary "north" that isn't parallel to up
    HMM_Vec3 world_y = (HMM_Vec3){{0.0f, 1.0f, 0.0f}};
    if (fabsf(vec3_dot(cam_dir, world_y)) > 0.99f) {
        world_y = (HMM_Vec3){{1.0f, 0.0f, 0.0f}};
    }

    ht->tangent_east = vec3_normalize(vec3_cross(world_y, cam_dir));
    ht->tangent_north = vec3_normalize(vec3_cross(cam_dir, ht->tangent_east));
}

// ---- Hex position computation ----

// Get the local tangent-plane 2D position of a hex at global offset (gcol, grow)
static void hex_local_pos(int gcol, int grow, float* out_x, float* out_z) {
    *out_x = gcol * HEX_COL_SPACING;
    *out_z = ((gcol & 1) ? (grow + 0.5f) : (float)grow) * HEX_ROW_SPACING;
}

// Convert a local tangent-plane position to world position on the sphere at given radius
static HMM_Vec3 tangent_to_world(const HexTerrain* ht, float lx, float lz, float radius) {
    HMM_Vec3 world_dir = vec3_add(ht->tangent_origin,
        vec3_add(vec3_scale(ht->tangent_east, lx),
                 vec3_scale(ht->tangent_north, lz)));
    return vec3_scale(vec3_normalize(world_dir), radius);
}

// ---- Mesh generation job ----

typedef struct HexMeshJob {
    // Chunk grid position
    int cx, cz;

    // Planet params
    float planet_radius;
    int seed;

    // Tangent frame (copied for thread safety)
    HMM_Vec3 tangent_origin;
    HMM_Vec3 tangent_up;
    HMM_Vec3 tangent_east;
    HMM_Vec3 tangent_north;

    // Floating origin
    float origin[3];

    // Output
    int16_t heights[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE];
    uint8_t types[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE];
    LodVertex* vertices;
    int vertex_count;

    // Linkage
    int chunk_index;
    volatile int completed;
} HexMeshJob;

// ---- Emit triangle helper ----

static void hex_emit_tri(LodVertex** verts, int* count, int* cap,
                          HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                          HMM_Vec3 normal, HMM_Vec3 color) {
    // Ensure correct winding (CCW from normal side)
    HMM_Vec3 e1 = vec3_sub(p1, p0);
    HMM_Vec3 e2 = vec3_sub(p2, p0);
    HMM_Vec3 face_n = vec3_cross(e1, e2);
    if (vec3_dot(face_n, normal) < 0.0f) {
        HMM_Vec3 tmp = p1; p1 = p2; p2 = tmp;
    }

    if (*count + 3 > *cap) {
        *cap = (*cap) * 2;
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

// ---- Core mesh generation for a chunk ----

static void generate_chunk_mesh(HexMeshJob* job) {
    fnl_state noise = ht_create_terrain_noise(job->seed);
    fnl_state warp = ht_create_warp_noise(job->seed);
    fnl_state detail = ht_create_detail_noise(job->seed);

    // Phase 1: Sample terrain heights for all columns in this chunk
    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            int gcol = job->cx * HEX_CHUNK_SIZE + col;
            int grow = job->cz * HEX_CHUNK_SIZE + row;

            float lx, lz;
            hex_local_pos(gcol, grow, &lx, &lz);

            // World direction for this hex center
            HMM_Vec3 world_dir = vec3_add(job->tangent_origin,
                vec3_add(vec3_scale(job->tangent_east, lx),
                         vec3_scale(job->tangent_north, lz)));
            HMM_Vec3 unit = vec3_normalize(world_dir);

            float h_m = ht_sample_height_m(&noise, &warp, &detail, unit);
            float effective_h = fmaxf(h_m, TERRAIN_SEA_LEVEL_M);

            // Discretize to integer layers (1m each).
            // Use ceil so hex prism tops are always at or above the smooth LOD surface.
            int h_layers = (int)ceilf(effective_h / HEX_HEIGHT);
            if (h_layers < 0) h_layers = 0;
            if (h_layers >= HEX_MAX_COLUMN_H) h_layers = HEX_MAX_COLUMN_H - 1;

            job->heights[col][row] = (int16_t)h_layers;
            job->types[col][row] = (uint8_t)ht_voxel_type(h_m);
        }
    }

    // Phase 2: Also sample a 1-cell border ring for neighbor wall checks
    // We store border heights in temporary arrays
    // border[dir][i] = height of neighbor in direction 'dir' for boundary cells
    // For simplicity, we'll re-sample neighbors on the fly during mesh gen

    // Phase 3: Generate hex prism geometry
    // Estimate max vertices: each hex can have top (4 tris=12 verts) + 6 walls (2 tris=12 verts each)
    // Worst case: 1024 * (12 + 72) = 86K verts. Typical much less.
    int cap = HEX_CHUNK_SIZE * HEX_CHUNK_SIZE * 24; // conservative initial
    job->vertices = (LodVertex*)malloc(cap * sizeof(LodVertex));
    job->vertex_count = 0;

    for (int col = 0; col < HEX_CHUNK_SIZE; col++) {
        for (int row = 0; row < HEX_CHUNK_SIZE; row++) {
            int gcol = job->cx * HEX_CHUNK_SIZE + col;
            int grow = job->cz * HEX_CHUNK_SIZE + row;
            int h = job->heights[col][row];

            float lx, lz;
            hex_local_pos(gcol, grow, &lx, &lz);

            // Compute hex center world direction
            HMM_Vec3 center_dir = vec3_normalize(
                vec3_add(job->tangent_origin,
                    vec3_add(vec3_scale(job->tangent_east, lx),
                             vec3_scale(job->tangent_north, lz))));

            float top_r = job->planet_radius + h * HEX_HEIGHT + HEX_SURFACE_BIAS;

            // Local up at this hex (radial direction)
            HMM_Vec3 local_up = center_dir;

            // Compute 6 vertex world positions for the top face
            HMM_Vec3 hex_verts_top[6];
            for (int i = 0; i < 6; i++) {
                float vx = lx + HEX_RADIUS * HEX_COS[i];
                float vz = lz + HEX_RADIUS * HEX_SIN[i];

                HMM_Vec3 vdir = vec3_normalize(
                    vec3_add(job->tangent_origin,
                        vec3_add(vec3_scale(job->tangent_east, vx),
                                 vec3_scale(job->tangent_north, vz))));
                hex_verts_top[i] = vec3_scale(vdir, top_r);
            }

            // -- Top face: 4 triangles (fan from vertex 0) --
            VoxelType surface_type = (VoxelType)job->types[col][row];
            HMM_Vec3 top_color = ht_voxel_color(surface_type);

            // Normal = local up (pointing outward from planet)
            hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                hex_verts_top[0], hex_verts_top[1], hex_verts_top[2],
                local_up, top_color);
            hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                hex_verts_top[0], hex_verts_top[2], hex_verts_top[3],
                local_up, top_color);
            hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                hex_verts_top[0], hex_verts_top[3], hex_verts_top[4],
                local_up, top_color);
            hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                hex_verts_top[0], hex_verts_top[4], hex_verts_top[5],
                local_up, top_color);

            // -- Side walls: check each of 6 neighbors --
            for (int dir = 0; dir < 6; dir++) {
                int ncol_g, nrow_g;
                hex_neighbor(gcol, grow, dir, &ncol_g, &nrow_g);

                // Get neighbor height
                int ncol_local = ncol_g - job->cx * HEX_CHUNK_SIZE;
                int nrow_local = nrow_g - job->cz * HEX_CHUNK_SIZE;

                int nh;
                if (ncol_local >= 0 && ncol_local < HEX_CHUNK_SIZE &&
                    nrow_local >= 0 && nrow_local < HEX_CHUNK_SIZE) {
                    nh = job->heights[ncol_local][nrow_local];
                } else {
                    // Neighbor is in another chunk — sample terrain directly
                    float nlx, nlz;
                    hex_local_pos(ncol_g, nrow_g, &nlx, &nlz);

                    HMM_Vec3 ndir = vec3_normalize(
                        vec3_add(job->tangent_origin,
                            vec3_add(vec3_scale(job->tangent_east, nlx),
                                     vec3_scale(job->tangent_north, nlz))));
                    float nh_m = ht_sample_height_m(&noise, &warp, &detail, ndir);
                    float nh_eff = fmaxf(nh_m, TERRAIN_SEA_LEVEL_M);
                    nh = (int)ceilf(nh_eff / HEX_HEIGHT);
                    if (nh < 0) nh = 0;
                }

                // Only emit wall if we're taller than the neighbor
                if (h <= nh) continue;

                // Wall edge: use DIR_TO_EDGE to get the hex edge facing this neighbor
                int edge = DIR_TO_EDGE[dir];
                int vi0 = edge;
                int vi1 = (edge + 1) % 6;

                // For each exposed layer, emit a wall quad
                int wall_bottom = nh;
                int wall_top = h;

                // Limit wall height to avoid excessive geometry
                if (wall_top - wall_bottom > 32) wall_bottom = wall_top - 32;

                for (int layer = wall_bottom; layer < wall_top; layer++) {
                    float bot_r = job->planet_radius + layer * HEX_HEIGHT + HEX_SURFACE_BIAS;
                    float top_r_wall = job->planet_radius + (layer + 1) * HEX_HEIGHT + HEX_SURFACE_BIAS;

                    // Recompute vertex positions at edge for this layer
                    float vx0 = lx + HEX_RADIUS * HEX_COS[vi0];
                    float vz0 = lz + HEX_RADIUS * HEX_SIN[vi0];
                    float vx1 = lx + HEX_RADIUS * HEX_COS[vi1];
                    float vz1 = lz + HEX_RADIUS * HEX_SIN[vi1];

                    HMM_Vec3 d0 = vec3_normalize(
                        vec3_add(job->tangent_origin,
                            vec3_add(vec3_scale(job->tangent_east, vx0),
                                     vec3_scale(job->tangent_north, vz0))));
                    HMM_Vec3 d1 = vec3_normalize(
                        vec3_add(job->tangent_origin,
                            vec3_add(vec3_scale(job->tangent_east, vx1),
                                     vec3_scale(job->tangent_north, vz1))));

                    HMM_Vec3 p0_bot = vec3_scale(d0, bot_r);
                    HMM_Vec3 p1_bot = vec3_scale(d1, bot_r);
                    HMM_Vec3 p0_top = vec3_scale(d0, top_r_wall);
                    HMM_Vec3 p1_top = vec3_scale(d1, top_r_wall);

                    // Wall normal: outward from hex center (perpendicular to edge, in tangent plane)
                    float edge_mid_x = (vx0 + vx1) * 0.5f;
                    float edge_mid_z = (vz0 + vz1) * 0.5f;
                    HMM_Vec3 wall_outward = vec3_normalize(
                        vec3_sub(
                            vec3_add(job->tangent_origin,
                                vec3_add(vec3_scale(job->tangent_east, edge_mid_x),
                                         vec3_scale(job->tangent_north, edge_mid_z))),
                            vec3_add(job->tangent_origin,
                                vec3_add(vec3_scale(job->tangent_east, lx),
                                         vec3_scale(job->tangent_north, lz)))
                        ));

                    float wall_y_m = (layer + 0.5f) * HEX_HEIGHT;
                    float surface_h_m = h * HEX_HEIGHT;
                    HMM_Vec3 wall_col = ht_wall_color(surface_h_m, wall_y_m);

                    // Two triangles for the wall quad
                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        p0_bot, p1_bot, p1_top, wall_outward, wall_col);
                    hex_emit_tri(&job->vertices, &job->vertex_count, &cap,
                        p0_bot, p1_top, p0_top, wall_outward, wall_col);
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

// ---- Chunk management ----

static int find_chunk(HexTerrain* ht, int cx, int cz) {
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

    if (chunk->pending_job) {
        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        // If still running, we can't free the job memory.
        // Mark chunk inactive; the job will be cleaned up when it completes.
        if (job->completed) {
            if (job->vertices) free(job->vertices);
            free(job);
        }
        // If not completed, leak the job (it'll finish and the memory will be orphaned).
        // In practice, chunk deactivation happens when walking away, and jobs complete quickly.
        chunk->pending_job = NULL;
    }

    chunk->active = false;
    chunk->dirty = false;
    chunk->generating = false;
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
    }

    printf("[HEX] Terrain system initialized: radius=%.0f, range=%.0fm, chunk=%dx%d (%.0fm x %.0fm)\n",
        planet_radius, HEX_RANGE, HEX_CHUNK_SIZE, HEX_CHUNK_SIZE,
        HEX_CHUNK_WIDTH, HEX_CHUNK_DEPTH);
    fflush(stdout);
}

void hex_terrain_destroy(HexTerrain* ht) {
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        if (ht->chunks[i].active) {
            free_chunk(&ht->chunks[i]);
        }
    }
}

// Threshold for re-anchoring the tangent frame (in meters).
// On an 800km planet, tangent-plane error at 10km is only ~6cm — negligible
// for 1m hex prisms. Set high so reanchor is rare; by the time it triggers,
// old chunks are out of activation range and already cleaned up.
#define FRAME_REANCHOR_THRESHOLD 10000.0f

void hex_terrain_update(HexTerrain* ht, HMM_Vec3 camera_pos,
                        const double world_origin[3]) {
    ht->camera_pos = camera_pos;
    ht->world_origin[0] = world_origin[0];
    ht->world_origin[1] = world_origin[1];
    ht->world_origin[2] = world_origin[2];

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
        // Recompute tangent frame centered on new camera position
        compute_tangent_frame(ht, camera_pos);
        ht->frame_anchor = camera_pos;
        ht->frame_valid = true;

        // Invalidate all existing chunks — they were in the old coordinate system
        for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
            if (ht->chunks[i].active) {
                free_chunk(&ht->chunks[i]);
            }
        }
    }

    // Process completed mesh generation jobs
    for (int i = 0; i < HEX_MAX_CHUNKS; i++) {
        HexChunk* chunk = &ht->chunks[i];
        if (!chunk->active || !chunk->generating || !chunk->pending_job) continue;

        HexMeshJob* job = (HexMeshJob*)chunk->pending_job;
        if (!job->completed) continue;

        // Transfer results
        chunk->cpu_vertices = job->vertices;
        chunk->cpu_vertex_count = job->vertex_count;
        memcpy(chunk->heights, job->heights, sizeof(chunk->heights));
        memcpy(chunk->types, job->types, sizeof(chunk->types));
        job->vertices = NULL;

        free(job);
        chunk->pending_job = NULL;
        chunk->generating = false;
        chunk->dirty = false;
    }

    // Determine chunk loading range.
    // Camera offset from tangent origin in the tangent plane:
    HMM_Vec3 cam_offset = vec3_sub(camera_pos, ht->tangent_origin);
    float cam_lx = vec3_dot(cam_offset, ht->tangent_east);
    float cam_lz = vec3_dot(cam_offset, ht->tangent_north);

    // Camera chunk coordinates
    int cam_cx = (int)floorf(cam_lx / HEX_CHUNK_WIDTH);
    int cam_cz = (int)floorf(cam_lz / HEX_CHUNK_DEPTH);

    // Number of chunks in each direction to cover HEX_RANGE
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

    // Activate/generate chunks within range (limited per frame to avoid stalls)
    int jobs_submitted = 0;
    int new_activations = 0;
    for (int cx = cam_cx - range_cx; cx <= cam_cx + range_cx; cx++) {
        for (int cz = cam_cz - range_cz; cz <= cam_cz + range_cz; cz++) {
            // Budget check: stop activating new chunks this frame if we hit the limit
            if (new_activations >= HEX_MAX_ACTIVATIONS || jobs_submitted >= 8) break;

            // Check if chunk center is within range of camera
            float chunk_center_x = (cx * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE / 2) * HEX_COL_SPACING;
            float chunk_center_z = (cz * HEX_CHUNK_SIZE + HEX_CHUNK_SIZE / 2) * HEX_ROW_SPACING;
            float dx = chunk_center_x - cam_lx;
            float dz = chunk_center_z - cam_lz;
            float dist = sqrtf(dx * dx + dz * dz);
            if (dist > HEX_RANGE + HEX_CHUNK_WIDTH) continue;

            int idx = find_chunk(ht, cx, cz);

            if (idx >= 0) {
                // Chunk exists — check if it needs meshing
                HexChunk* chunk = &ht->chunks[idx];
                if (chunk->gpu_vertex_count == 0 && !chunk->generating &&
                    chunk->cpu_vertex_count == 0) {
                    chunk->dirty = true;
                }
            } else {
                // New chunk needed — respect per-frame activation limit
                if (new_activations >= HEX_MAX_ACTIVATIONS) continue;

                idx = find_free_chunk(ht);
                if (idx < 0) continue;

                HexChunk* chunk = &ht->chunks[idx];
                memset(chunk, 0, sizeof(HexChunk));
                chunk->cx = cx;
                chunk->cz = cz;
                chunk->active = true;
                chunk->dirty = true;
                chunk->gpu_buffer.id = SG_INVALID_ID;
                new_activations++;
            }

            // Submit mesh gen job for dirty chunks
            if (idx >= 0 && idx < HEX_MAX_CHUNKS && jobs_submitted < 8) {
                HexChunk* chunk = &ht->chunks[idx];
                if (chunk->dirty && !chunk->generating) {
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
                    job->vertices = NULL;
                    job->vertex_count = 0;
                    job->completed = 0;

                    if (job_system_try_submit(ht->jobs, hex_mesh_worker, job)) {
                        chunk->pending_job = job;
                        chunk->generating = true;
                        chunk->dirty = false;
                        jobs_submitted++;
                    } else {
                        free(job);
                    }
                }
            }
        }
        if (new_activations >= HEX_MAX_ACTIVATIONS || jobs_submitted >= 8) break;
    }

    // Count active chunks
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

        // Destroy old buffer
        if (chunk->gpu_buffer.id != SG_INVALID_ID) {
            sg_destroy_buffer(chunk->gpu_buffer);
        }

        chunk->gpu_buffer = sg_make_buffer(&(sg_buffer_desc){
            .data = (sg_range){
                chunk->cpu_vertices,
                (size_t)chunk->cpu_vertex_count * sizeof(LodVertex)
            },
            .label = "hex-chunk-vertices",
        });

        if (chunk->gpu_buffer.id == SG_INVALID_ID) {
            continue; // Retry next frame
        }

        chunk->gpu_vertex_count = chunk->cpu_vertex_count;
        uploads++;

        // Free CPU data
        free(chunk->cpu_vertices);
        chunk->cpu_vertices = NULL;
        chunk->cpu_vertex_count = 0;
    }
}

void hex_terrain_render(HexTerrain* ht, sg_pipeline pip) {
    (void)pip; // Pipeline already applied by caller

    sg_bindings bind = {0};
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
    // Check if the position is within HEX_RANGE of the camera
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
        total++;
        if (ht->chunks[i].gpu_vertex_count > 0) meshed++;
    }
    if (total == 0) return 0.0f;
    float coverage = (float)meshed / (float)total;
    // Only suppress LOD when we have good hex terrain coverage.
    // Below 80%, LOD renders underneath to prevent black holes.
    return (coverage >= 0.8f) ? HEX_RANGE : 0.0f;
}

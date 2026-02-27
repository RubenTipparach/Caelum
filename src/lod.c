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

    // Output (written by worker thread)
    LodVertex* vertices;
    int vertex_count;

    // Linkage
    int node_index;
    volatile int completed;  // 0 = pending, 1 = done
} MeshGenJob;

// ---- Terrain noise (same config as planet.c) ----

static fnl_state create_terrain_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 4;
    noise.frequency = 2.0f;
    noise.seed = seed;
    return noise;
}

// Sample terrain height at a point on the unit sphere
// Returns height in layers (same scale as planet.c)
static int sample_terrain_height(fnl_state* noise, HMM_Vec3 unit_pos, int sea_level) {
    float n = fnlGetNoise3D(noise,
        unit_pos.X * 3.0f, unit_pos.Y * 3.0f, unit_pos.Z * 3.0f);
    int height = 8 + (int)((n + 1.0f) * 0.5f * 40.0f);
    if (height < 1) height = 1;
    if (height >= 64) height = 63;
    (void)sea_level;
    return height;
}

// Get terrain color based on height relative to sea level (matches planet.c)
static HMM_Vec3 terrain_color(int height, int sea_level) {
    int rel = height - sea_level;
    if (height < sea_level) {
        // Underwater: show water color
        return (HMM_Vec3){{0.15f, 0.35f, 0.65f}};
    } else if (rel <= 1) {
        return (HMM_Vec3){{0.76f, 0.70f, 0.40f}};  // Sand
    } else if (rel <= 6) {
        return (HMM_Vec3){{0.20f, 0.60f, 0.15f}};  // Grass
    } else if (rel <= 12) {
        return (HMM_Vec3){{0.50f, 0.50f, 0.50f}};  // Stone
    } else {
        return (HMM_Vec3){{0.90f, 0.95f, 1.00f}};  // Ice
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

// Get color for a voxel type (thread-safe, matches planet.c voxel_color)
static HMM_Vec3 voxel_type_color(VoxelType type) {
    switch (type) {
        case VOXEL_WATER: return (HMM_Vec3){{0.15f, 0.35f, 0.65f}};
        case VOXEL_SAND:  return (HMM_Vec3){{0.76f, 0.70f, 0.40f}};
        case VOXEL_DIRT:  return (HMM_Vec3){{0.45f, 0.32f, 0.18f}};
        case VOXEL_GRASS: return (HMM_Vec3){{0.20f, 0.60f, 0.15f}};
        case VOXEL_STONE: return (HMM_Vec3){{0.50f, 0.50f, 0.50f}};
        case VOXEL_ICE:   return (HMM_Vec3){{0.90f, 0.95f, 1.00f}};
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
    if (tree->node_count >= tree->node_capacity) {
        // Grow the array
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
    tree->active_leaf_count = 0;
    tree->total_vertex_count = 0;
    tree->stats_frame_counter = 0;

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

// ---- Mesh generation ----

// Helper: emit a triangle into the vertex array
static void lod_emit_tri(LodVertex** verts, int* count, int* cap,
                          HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                          HMM_Vec3 normal, HMM_Vec3 color) {
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

// Generate mesh for a coarse node (tessellated spherical triangle).
// Thread-safe: takes all inputs by value, writes to output pointers.
static void generate_coarse_mesh_params(
    const SphericalTriangle* tri, int depth,
    float planet_radius, float layer_thickness, int sea_level, int seed,
    LodVertex** out_vertices, int* out_count)
{
    // Tessellation resolution depends on depth
    int tess = 4;
    if (depth >= 4) tess = 6;
    if (depth >= 8) tess = 8;

    fnl_state noise = create_terrain_noise(seed);

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

            int h = sample_terrain_height(&noise, p, sea_level);
            int effective_h = h;
            if (effective_h < sea_level) effective_h = sea_level;

            float radius = planet_radius +
                           (effective_h + 1) * layer_thickness;

            points[idx] = vec3_scale(p, radius);
            heights[idx] = (float)h;
            colors[idx] = terrain_color(h, sea_level);
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
    // Extend each boundary edge downward by a few layers to cover seams
    float skirt_drop = 3.0f * layer_thickness;

    // Helper macro to emit a skirt wall quad between two boundary vertices
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

    // Edge 1: v0 -> v1 (bottom row, row=0)
    for (int col = 0; col < tess; col++) {
        EMIT_SKIRT(row_offsets[0] + col, row_offsets[0] + col + 1);
    }

    // Edge 2: v2 -> v0 (left side, col=0 in each row)
    for (int row = 0; row < tess; row++) {
        EMIT_SKIRT(row_offsets[row + 1], row_offsets[row]);
    }

    // Edge 3: v1 -> v2 (hypotenuse, last vertex in each row)
    for (int row = 0; row < tess; row++) {
        EMIT_SKIRT(row_offsets[row] + (tess - row),
                   row_offsets[row + 1] + (tess - row - 1));
    }

    #undef EMIT_SKIRT

    free(points);
    free(heights);
    free(colors);
    free(row_offsets);
}

// ---- Voxel mesh generation for fine LOD nodes (depth >= LOD_VOXEL_DEPTH) ----
// Each tessellation triangle becomes a discrete-height "cell" with wall faces
// between adjacent cells at different heights, giving a Minecraft-like stepped look.

// Edge info for tracking adjacency between tessellation triangles
typedef struct {
    int va, vb;     // Sorted vertex indices (va < vb)
    int tri[2];     // Up to 2 triangles sharing this edge
    int count;      // Number of triangles (1 = boundary, 2 = interior)
} VoxelEdge;

// Find or create an edge entry in the flat array
static int find_or_add_edge(VoxelEdge* edges, int* edge_count, int cap, int va, int vb) {
    int a = va < vb ? va : vb;
    int b = va < vb ? vb : va;
    for (int i = 0; i < *edge_count; i++) {
        if (edges[i].va == a && edges[i].vb == b) {
            return i;
        }
    }
    if (*edge_count >= cap) return -1;
    int idx = (*edge_count)++;
    edges[idx].va = a;
    edges[idx].vb = b;
    edges[idx].count = 0;
    edges[idx].tri[0] = -1;
    edges[idx].tri[1] = -1;
    return idx;
}

static void generate_voxel_mesh_params(
    const SphericalTriangle* stri, int depth,
    float planet_radius, float layer_thickness, int sea_level, int seed,
    LodVertex** out_vertices, int* out_count)
{
    (void)depth;
    int tess = 8;  // Fine tessellation for voxel detail
    fnl_state noise = create_terrain_noise(seed);

    // 1. Generate grid vertices on unit sphere
    int num_pts = (tess + 1) * (tess + 2) / 2;
    HMM_Vec3* unit_pts = (HMM_Vec3*)malloc(num_pts * sizeof(HMM_Vec3));

    int pidx = 0;
    for (int row = 0; row <= tess; row++) {
        for (int col = 0; col <= tess - row; col++) {
            float u = (float)col / (float)tess;
            float v = (float)row / (float)tess;
            float w = 1.0f - u - v;
            HMM_Vec3 p = vec3_add(
                vec3_add(vec3_scale(stri->v0, w), vec3_scale(stri->v1, u)),
                vec3_scale(stri->v2, v));
            unit_pts[pidx++] = vec3_normalize(p);
        }
    }

    // 2. Build row offsets for grid indexing
    int* row_off = (int*)malloc((tess + 1) * sizeof(int));
    row_off[0] = 0;
    for (int r = 1; r <= tess; r++) {
        row_off[r] = row_off[r - 1] + (tess - (r - 1) + 1);
    }

    // 3. Build triangle list with terrain heights at centroids
    int num_tris = tess * tess;
    int* tri_v = (int*)malloc(num_tris * 3 * sizeof(int));
    int* tri_h = (int*)malloc(num_tris * sizeof(int));   // raw terrain height
    int* tri_eh = (int*)malloc(num_tris * sizeof(int));  // effective height (clamped to sea level)

    int tidx = 0;
    for (int row = 0; row < tess; row++) {
        for (int col = 0; col < tess - row; col++) {
            // Upward triangle: V(row,col), V(row,col+1), V(row+1,col)
            int i0 = row_off[row] + col;
            int i1 = row_off[row] + col + 1;
            int i2 = row_off[row + 1] + col;
            tri_v[tidx * 3 + 0] = i0;
            tri_v[tidx * 3 + 1] = i1;
            tri_v[tidx * 3 + 2] = i2;

            HMM_Vec3 c = vec3_normalize(vec3_scale(
                vec3_add(vec3_add(unit_pts[i0], unit_pts[i1]), unit_pts[i2]),
                1.0f / 3.0f));
            int h = sample_terrain_height(&noise, c, sea_level);
            tri_h[tidx] = h;
            tri_eh[tidx] = h < sea_level ? sea_level : h;
            tidx++;

            // Downward triangle: V(row+1,col), V(row,col+1), V(row+1,col+1)
            if (col < tess - row - 1) {
                int j0 = row_off[row + 1] + col;
                int j1 = row_off[row] + col + 1;
                int j2 = row_off[row + 1] + col + 1;
                tri_v[tidx * 3 + 0] = j0;
                tri_v[tidx * 3 + 1] = j1;
                tri_v[tidx * 3 + 2] = j2;

                c = vec3_normalize(vec3_scale(
                    vec3_add(vec3_add(unit_pts[j0], unit_pts[j1]), unit_pts[j2]),
                    1.0f / 3.0f));
                h = sample_terrain_height(&noise, c, sea_level);
                tri_h[tidx] = h;
                tri_eh[tidx] = h < sea_level ? sea_level : h;
                tidx++;
            }
        }
    }

    // 4. Build edge adjacency
    int edge_cap = num_tris * 2;
    VoxelEdge* edges = (VoxelEdge*)calloc(edge_cap, sizeof(VoxelEdge));
    int edge_count = 0;

    for (int t = 0; t < tidx; t++) {
        for (int e = 0; e < 3; e++) {
            int va = tri_v[t * 3 + e];
            int vb = tri_v[t * 3 + (e + 1) % 3];
            int ei = find_or_add_edge(edges, &edge_count, edge_cap, va, vb);
            if (ei >= 0 && edges[ei].count < 2) {
                edges[ei].tri[edges[ei].count++] = t;
            }
        }
    }

    // 5. Allocate output
    int max_verts = num_tris * 3 * 20;  // Extra room for wall quads and skirts
    *out_vertices = (LodVertex*)malloc(max_verts * sizeof(LodVertex));
    *out_count = 0;
    int cap = max_verts;

    // 6. Render top faces (each triangle at its discrete height)
    for (int t = 0; t < tidx; t++) {
        int eh = tri_eh[t];
        float r_top = planet_radius + (eh + 1) * layer_thickness;

        HMM_Vec3 p0 = vec3_scale(unit_pts[tri_v[t * 3 + 0]], r_top);
        HMM_Vec3 p1 = vec3_scale(unit_pts[tri_v[t * 3 + 1]], r_top);
        HMM_Vec3 p2 = vec3_scale(unit_pts[tri_v[t * 3 + 2]], r_top);

        HMM_Vec3 edge1 = vec3_sub(p1, p0);
        HMM_Vec3 edge2 = vec3_sub(p2, p0);
        HMM_Vec3 normal = vec3_normalize(vec3_cross(edge1, edge2));

        // Ensure normal points outward (away from planet center)
        HMM_Vec3 fc = vec3_scale(vec3_add(vec3_add(p0, p1), p2), 1.0f / 3.0f);
        if (vec3_dot(normal, fc) < 0) {
            normal = vec3_scale(normal, -1.0f);
            HMM_Vec3 tmp = p1; p1 = p2; p2 = tmp;
        }

        VoxelType vt = compute_voxel_type(eh, tri_h[t], sea_level);
        HMM_Vec3 color = voxel_type_color(vt);

        lod_emit_tri(out_vertices, out_count, &cap, p0, p1, p2, normal, color);
    }

    // 7. Render walls between adjacent triangles and boundary skirts
    for (int i = 0; i < edge_count; i++) {
        VoxelEdge* e = &edges[i];
        HMM_Vec3 ua = unit_pts[e->va];
        HMM_Vec3 ub = unit_pts[e->vb];

        if (e->count == 1) {
            // Boundary edge: render skirt extending downward
            int t = e->tri[0];
            int eh = tri_eh[t];
            int th = tri_h[t];
            int skirt_bottom = eh - 3;
            if (skirt_bottom < 0) skirt_bottom = 0;

            // Normal points outward from the patch (away from triangle centroid)
            HMM_Vec3 tri_c = vec3_normalize(vec3_scale(
                vec3_add(vec3_add(
                    unit_pts[tri_v[t * 3 + 0]],
                    unit_pts[tri_v[t * 3 + 1]]),
                    unit_pts[tri_v[t * 3 + 2]]), 1.0f / 3.0f));
            HMM_Vec3 edge_mid = vec3_normalize(vec3_scale(vec3_add(ua, ub), 0.5f));
            HMM_Vec3 outward = vec3_normalize(vec3_sub(edge_mid, tri_c));

            for (int layer = skirt_bottom; layer <= eh; layer++) {
                float r_bot = planet_radius + layer * layer_thickness;
                float r_top = planet_radius + (layer + 1) * layer_thickness;

                HMM_Vec3 bl = vec3_scale(ua, r_bot);
                HMM_Vec3 br = vec3_scale(ub, r_bot);
                HMM_Vec3 tl = vec3_scale(ua, r_top);
                HMM_Vec3 tr = vec3_scale(ub, r_top);

                VoxelType vt = compute_voxel_type(layer, th, sea_level);
                HMM_Vec3 color = voxel_type_color(vt);

                lod_emit_tri(out_vertices, out_count, &cap, bl, tr, br, outward, color);
                lod_emit_tri(out_vertices, out_count, &cap, bl, tl, tr, outward, color);
            }
        } else {
            // Interior edge: render wall if heights differ
            int t0 = e->tri[0];
            int t1 = e->tri[1];
            int h0 = tri_eh[t0];
            int h1 = tri_eh[t1];

            if (h0 == h1) continue;  // Same height, no wall needed

            int tall = h0 > h1 ? t0 : t1;
            int shor = h0 > h1 ? t1 : t0;
            int high = h0 > h1 ? h0 : h1;
            int low = h0 > h1 ? h1 : h0;
            int tall_th = h0 > h1 ? tri_h[t0] : tri_h[t1];

            // Normal from tall toward short (faces the exposed cliff)
            HMM_Vec3 tall_c = vec3_normalize(vec3_scale(
                vec3_add(vec3_add(
                    unit_pts[tri_v[tall * 3 + 0]],
                    unit_pts[tri_v[tall * 3 + 1]]),
                    unit_pts[tri_v[tall * 3 + 2]]), 1.0f / 3.0f));
            HMM_Vec3 short_c = vec3_normalize(vec3_scale(
                vec3_add(vec3_add(
                    unit_pts[tri_v[shor * 3 + 0]],
                    unit_pts[tri_v[shor * 3 + 1]]),
                    unit_pts[tri_v[shor * 3 + 2]]), 1.0f / 3.0f));
            HMM_Vec3 wall_normal = vec3_normalize(vec3_sub(short_c, tall_c));

            for (int layer = low + 1; layer <= high; layer++) {
                float r_bot = planet_radius + layer * layer_thickness;
                float r_top = planet_radius + (layer + 1) * layer_thickness;

                HMM_Vec3 bl = vec3_scale(ua, r_bot);
                HMM_Vec3 br = vec3_scale(ub, r_bot);
                HMM_Vec3 tl = vec3_scale(ua, r_top);
                HMM_Vec3 tr = vec3_scale(ub, r_top);

                VoxelType vt = compute_voxel_type(layer, tall_th, sea_level);
                HMM_Vec3 color = voxel_type_color(vt);

                lod_emit_tri(out_vertices, out_count, &cap, bl, tr, br, wall_normal, color);
                lod_emit_tri(out_vertices, out_count, &cap, bl, tl, tr, wall_normal, color);
            }
        }
    }

    // Cleanup
    free(unit_pts);
    free(row_off);
    free(tri_v);
    free(tri_h);
    free(tri_eh);
    free(edges);
}

// Worker thread callback for async mesh generation
static void mesh_gen_worker(void* data) {
    MeshGenJob* job = (MeshGenJob*)data;

    if (job->depth >= LOD_VOXEL_DEPTH) {
        generate_voxel_mesh_params(
            &job->tri, job->depth,
            job->planet_radius, job->layer_thickness, job->sea_level, job->seed,
            &job->vertices, &job->vertex_count);
    } else {
        generate_coarse_mesh_params(
            &job->tri, job->depth,
            job->planet_radius, job->layer_thickness, job->sea_level, job->seed,
            &job->vertices, &job->vertex_count);
    }

    // Signal completion (volatile write - safe on x86 with store ordering)
    job->completed = 1;
}

// Synchronous mesh generation (fallback when no job system)
static void generate_node_mesh(LodTree* tree, LodNode* node) {
    if (node->depth >= LOD_VOXEL_DEPTH) {
        generate_voxel_mesh_params(
            &node->tri, node->depth,
            tree->planet_radius, tree->layer_thickness, tree->sea_level, tree->seed,
            &node->cpu_vertices, &node->cpu_vertex_count);
    } else {
        generate_coarse_mesh_params(
            &node->tri, node->depth,
            tree->planet_radius, tree->layer_thickness, tree->sea_level, tree->seed,
            &node->cpu_vertices, &node->cpu_vertex_count);
    }
    node->state = LOD_READY;
}

// Submit async mesh generation job for a node
static void submit_mesh_job(LodTree* tree, int node_idx) {
    LodNode* node = &tree->nodes[node_idx];

    MeshGenJob* job = (MeshGenJob*)calloc(1, sizeof(MeshGenJob));
    job->tri = node->tri;
    job->depth = node->depth;
    job->planet_radius = tree->planet_radius;
    job->layer_thickness = tree->layer_thickness;
    job->sea_level = tree->sea_level;
    job->seed = tree->seed;
    job->node_index = node_idx;
    job->vertices = NULL;
    job->vertex_count = 0;
    job->completed = 0;

    node->pending_job = job;
    node->state = LOD_GENERATING;

    job_system_submit(tree->jobs, mesh_gen_worker, job);
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

        // Transfer mesh data from job to node
        node->cpu_vertices = job->vertices;
        node->cpu_vertex_count = job->vertex_count;
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

    // Free children's resources
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

        // Note: we don't compact the array (would invalidate indices).
        // Dead nodes stay with state=LOD_UNLOADED.
        node->children[i] = -1;
    }

    node->is_leaf = true;
}

// ---- Screen-space metric ----

// Compute how large this node appears on screen.
// Uses distance to NEAREST EDGE of the patch, not center.
// This is critical: when the camera is inside a patch (standing on it),
// the center may be hundreds of km away, but the patch extends beneath us.
static float compute_screen_metric(const LodTree* tree, const LodNode* node) {
    float patch_size = node->tri.angular_radius * tree->planet_radius;

    // Angular distance from camera direction to patch center on the unit sphere
    float cam_r = sqrtf(vec3_dot(tree->camera_pos, tree->camera_pos));
    if (cam_r < 1.0f) cam_r = 1.0f;
    HMM_Vec3 cam_dir = vec3_scale(tree->camera_pos, 1.0f / cam_r);
    float cos_angle = vec3_dot(cam_dir, node->tri.center);
    if (cos_angle > 1.0f) cos_angle = 1.0f;
    if (cos_angle < -1.0f) cos_angle = -1.0f;
    float angle_to_center = acosf(cos_angle);

    // Distance to nearest edge of the patch (0 if camera is above the patch)
    float edge_angle = angle_to_center - node->tri.angular_radius;
    if (edge_angle < 0.0f) edge_angle = 0.0f;

    // Camera altitude above the planet surface
    float altitude = cam_r - tree->planet_radius;
    if (altitude < 1.0f) altitude = 1.0f;

    // 3D distance to nearest patch edge: Pythagorean with surface arc + altitude
    float surface_dist = edge_angle * tree->planet_radius;
    float dist = sqrtf(surface_dist * surface_dist + altitude * altitude);

    return patch_size / dist;
}

// ---- Update (split/merge based on camera) ----

// Recursive traversal: decide split/merge for each node
static void update_node(LodTree* tree, int node_idx) {
    LodNode* node = &tree->nodes[node_idx];
    float metric = compute_screen_metric(tree, node);

    if (node->is_leaf) {
        // Submit async mesh generation if needed
        if (node->state == LOD_UNLOADED) {
            if (tree->jobs) {
                submit_mesh_job(tree, node_idx);
            } else {
                // Fallback: synchronous generation
                generate_node_mesh(tree, node);
            }
        }

        // Only split nodes with GPU buffers (ensures parent can be fallback)
        if (node->state == LOD_ACTIVE) {
            if (metric > LOD_SPLIT_THRESHOLD && node->depth < LOD_MAX_DEPTH) {
                if (tree->active_leaf_count < LOD_TARGET_LEAVES * 2 &&
                    tree->splits_this_frame < LOD_MAX_SPLITS_PER_FRAME) {
                    split_node(tree, node_idx);
                    tree->splits_this_frame++;
                    // Recurse into new children
                    // Must re-fetch node pointer EACH iteration — recursive update_node
                    // may trigger split_node → allocate_node → realloc, invalidating pointers.
                    if (!tree->nodes[node_idx].is_leaf) {
                        for (int i = 0; i < LOD_CHILDREN; i++) {
                            int ci = tree->nodes[node_idx].children[i];
                            if (ci >= 0) {
                                update_node(tree, ci);
                            }
                        }
                    }
                    node = &tree->nodes[node_idx];
                }
            }
        }
    } else {
        // Internal node: check if we should merge children
        bool all_children_want_merge = true;
        for (int i = 0; i < LOD_CHILDREN; i++) {
            int child_idx = node->children[i];
            if (child_idx < 0) continue;

            LodNode* child = &tree->nodes[child_idx];
            if (!child->is_leaf) {
                all_children_want_merge = false;
                break;
            }
            // Don't merge if any child is still generating
            if (child->state == LOD_GENERATING) {
                all_children_want_merge = false;
                break;
            }
            float child_metric = compute_screen_metric(tree, child);
            if (child_metric > LOD_MERGE_THRESHOLD) {
                all_children_want_merge = false;
                break;
            }
        }

        if (all_children_want_merge && metric < LOD_MERGE_THRESHOLD) {
            merge_node(tree, node_idx);
            // Parent is now a leaf; if it still has a valid mesh, use it
            node = &tree->nodes[node_idx];
            if (node->state == LOD_UNLOADED) {
                if (tree->jobs) {
                    submit_mesh_job(tree, node_idx);
                } else {
                    generate_node_mesh(tree, node);
                }
            }
        } else {
            // Recurse into children
            // Must re-fetch via tree->nodes[] each iteration — recursive calls may realloc.
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
    (void)view_proj;  // May use for frustum culling later
    tree->camera_pos = camera_pos;
    tree->active_leaf_count = 0;
    tree->splits_this_frame = 0;

    // Process completed async jobs first (transfers mesh data to nodes)
    process_completed_jobs(tree);

    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        update_node(tree, tree->root_nodes[i]);
    }

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

        // Compute distance from camera to patch center
        HMM_Vec3 world_center = HMM_MulV3F(node->tri.center, tree->planet_radius);
        HMM_Vec3 diff = HMM_SubV3(world_center, tree->camera_pos);
        float dist = HMM_LenV3(diff);

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

        node->gpu_buffer = sg_make_buffer(&(sg_buffer_desc){
            .data = (sg_range){
                node->cpu_vertices,
                (size_t)node->cpu_vertex_count * sizeof(LodVertex)
            },
            .label = "lod-patch-vertices",
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

// Recursive render: traverses the tree from each node downward.
// An internal node renders itself as fallback if any child is not yet ACTIVE,
// preventing temporal popping (parent stays visible until all children are ready).
static void render_node(const LodTree* tree, int node_idx,
                         sg_bindings* bind, int* total_verts,
                         LodPreDrawFunc pre_draw, void* user_data) {
    const LodNode* node = &tree->nodes[node_idx];

    if (node->is_leaf) {
        // Leaf: render if active
        if (node->state == LOD_ACTIVE && node->gpu_vertex_count > 0 &&
            node->gpu_buffer.id != SG_INVALID_ID) {
            if (pre_draw) pre_draw(node->depth, user_data);
            bind->vertex_buffers[0] = node->gpu_buffer;
            sg_apply_bindings(bind);
            sg_draw(0, node->gpu_vertex_count, 1);
            *total_verts += node->gpu_vertex_count;
        }
        return;
    }

    // Internal node: check if all children have ACTIVE meshes
    bool all_children_active = true;
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int ci = node->children[i];
        if (ci < 0 || tree->nodes[ci].state != LOD_ACTIVE) {
            all_children_active = false;
            break;
        }
    }

    if (all_children_active) {
        // All children ready: recurse into them (don't render parent)
        for (int i = 0; i < LOD_CHILDREN; i++) {
            if (node->children[i] >= 0) {
                render_node(tree, node->children[i], bind, total_verts,
                           pre_draw, user_data);
            }
        }
    } else {
        // Not all children ready: render this node as fallback
        // (parent's mesh covers all children's regions until they're ready)
        if (node->state == LOD_ACTIVE && node->gpu_vertex_count > 0 &&
            node->gpu_buffer.id != SG_INVALID_ID) {
            if (pre_draw) pre_draw(node->depth, user_data);
            bind->vertex_buffers[0] = node->gpu_buffer;
            sg_apply_bindings(bind);
            sg_draw(0, node->gpu_vertex_count, 1);
            *total_verts += node->gpu_vertex_count;
        }
    }
}

void lod_tree_render(LodTree* tree, sg_pipeline pip, HMM_Mat4 vp,
                     LodPreDrawFunc pre_draw, void* user_data) {
    (void)pip;
    (void)vp;

    sg_bindings bind = {0};
    int total_verts = 0;

    // Traverse from each root node recursively
    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        render_node(tree, tree->root_nodes[i], &bind, &total_verts,
                   pre_draw, user_data);
    }

    tree->total_vertex_count = total_verts;
}

// ---- Spatial queries ----

static int find_node_recursive(const LodTree* tree, int node_idx, HMM_Vec3 unit_pos) {
    const LodNode* node = &tree->nodes[node_idx];

    if (node->is_leaf) {
        return node_idx;
    }

    // Find which child contains the point
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int child_idx = node->children[i];
        if (child_idx < 0) continue;
        if (point_in_spherical_triangle(&tree->nodes[child_idx].tri, unit_pos)) {
            return find_node_recursive(tree, child_idx, unit_pos);
        }
    }

    // Fallback: closest child center
    float best_dot = -2.0f;
    int best_child = node_idx;
    for (int i = 0; i < LOD_CHILDREN; i++) {
        int child_idx = node->children[i];
        if (child_idx < 0) continue;
        float d = vec3_dot(tree->nodes[child_idx].tri.center, unit_pos);
        if (d > best_dot) {
            best_dot = d;
            best_child = child_idx;
        }
    }
    return find_node_recursive(tree, best_child, unit_pos);
}

int lod_tree_find_node(const LodTree* tree, HMM_Vec3 world_pos) {
    HMM_Vec3 unit = vec3_normalize(world_pos);

    // Find which root triangle contains this point
    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        int root = tree->root_nodes[i];
        if (point_in_spherical_triangle(&tree->nodes[root].tri, unit)) {
            return find_node_recursive(tree, root, unit);
        }
    }

    // Fallback: closest root center
    float best_dot = -2.0f;
    int best_root = 0;
    for (int i = 0; i < LOD_ROOT_COUNT; i++) {
        int root = tree->root_nodes[i];
        float d = vec3_dot(tree->nodes[root].tri.center, unit);
        if (d > best_dot) {
            best_dot = d;
            best_root = root;
        }
    }
    return find_node_recursive(tree, best_root, unit);
}

float lod_tree_terrain_height(const LodTree* tree, HMM_Vec3 world_pos) {
    HMM_Vec3 unit = vec3_normalize(world_pos);
    fnl_state noise = create_terrain_noise(tree->seed);
    int h = sample_terrain_height(&noise, unit, tree->sea_level);
    int effective_h = h;
    if (effective_h < tree->sea_level) effective_h = tree->sea_level;
    return tree->planet_radius + (effective_h + 1) * tree->layer_thickness;
}

int lod_tree_active_leaves(const LodTree* tree) {
    return tree->active_leaf_count;
}

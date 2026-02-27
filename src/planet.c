#include "planet.h"
#include "math_utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#define FNL_IMPL
#include "FastNoiseLite.h"

// ---- Geodesic sphere subdivision + dual construction ----

typedef struct {
    HMM_Vec3* verts;
    int vert_count;
    int vert_cap;
    int* tris;       // 3 indices per triangle
    int tri_count;
    int tri_cap;
} GeoMesh;

static void geo_init(GeoMesh* g, int vert_cap, int tri_cap) {
    g->verts = (HMM_Vec3*)calloc(vert_cap, sizeof(HMM_Vec3));
    g->vert_count = 0;
    g->vert_cap = vert_cap;
    g->tris = (int*)calloc(tri_cap * 3, sizeof(int));
    g->tri_count = 0;
    g->tri_cap = tri_cap;
}

static void geo_free(GeoMesh* g) {
    free(g->verts);
    free(g->tris);
}

static int geo_add_vert(GeoMesh* g, HMM_Vec3 v) {
    if (g->vert_count >= g->vert_cap) {
        g->vert_cap *= 2;
        g->verts = (HMM_Vec3*)realloc(g->verts, g->vert_cap * sizeof(HMM_Vec3));
    }
    g->verts[g->vert_count] = v;
    return g->vert_count++;
}

static void geo_add_tri(GeoMesh* g, int a, int b, int c) {
    if (g->tri_count >= g->tri_cap) {
        g->tri_cap *= 2;
        g->tris = (int*)realloc(g->tris, g->tri_cap * 3 * sizeof(int));
    }
    int i = g->tri_count * 3;
    g->tris[i] = a;
    g->tris[i + 1] = b;
    g->tris[i + 2] = c;
    g->tri_count++;
}

// ---- Vertex deduplication via hash map ----

typedef struct {
    uint64_t key;
    int value;
} EdgeEntry;

typedef struct {
    EdgeEntry* entries;
    int capacity;
} EdgeMap;

static void edgemap_init(EdgeMap* m, int capacity) {
    m->capacity = capacity;
    m->entries = (EdgeEntry*)calloc(capacity, sizeof(EdgeEntry));
}

static void edgemap_free(EdgeMap* m) {
    free(m->entries);
}

static uint64_t edge_key(int a, int b) {
    if (a > b) { int t = a; a = b; b = t; }
    return ((uint64_t)(a + 1) << 32) | (uint64_t)(b + 1);
}

static int edgemap_get(EdgeMap* m, uint64_t key) {
    int idx = (int)(key % (uint64_t)m->capacity);
    for (int i = 0; i < m->capacity; i++) {
        int slot = (idx + i) % m->capacity;
        if (m->entries[slot].key == 0) return -1;
        if (m->entries[slot].key == key) return m->entries[slot].value;
    }
    return -1;
}

static void edgemap_put(EdgeMap* m, uint64_t key, int value) {
    int idx = (int)(key % (uint64_t)m->capacity);
    for (int i = 0; i < m->capacity; i++) {
        int slot = (idx + i) % m->capacity;
        if (m->entries[slot].key == 0 || m->entries[slot].key == key) {
            m->entries[slot].key = key;
            m->entries[slot].value = value;
            return;
        }
    }
}

// ---- Subdivide the icosahedron ----

static int get_midpoint(GeoMesh* g, EdgeMap* em, int a, int b) {
    uint64_t key = edge_key(a, b);
    int existing = edgemap_get(em, key);
    if (existing >= 0) return existing;

    HMM_Vec3 mid = vec3_normalize(vec3_scale(
        vec3_add(g->verts[a], g->verts[b]), 0.5f));
    int idx = geo_add_vert(g, mid);
    edgemap_put(em, key, idx);
    return idx;
}

static void build_geodesic(GeoMesh* g, int subdivisions) {
    for (int i = 0; i < ICO_VERTEX_COUNT; i++) {
        geo_add_vert(g, vec3_normalize(ICO_VERTICES[i]));
    }
    for (int i = 0; i < ICO_FACE_COUNT; i++) {
        geo_add_tri(g, ICO_FACES[i][0], ICO_FACES[i][1], ICO_FACES[i][2]);
    }

    int iters = 0;
    int m = subdivisions;
    while (m > 1) { m /= 2; iters++; }

    for (int iter = 0; iter < iters; iter++) {
        int old_tri_count = g->tri_count;
        EdgeMap em;
        edgemap_init(&em, old_tri_count * 4);

        for (int t = 0; t < old_tri_count; t++) {
            int a = g->tris[t * 3];
            int b = g->tris[t * 3 + 1];
            int c = g->tris[t * 3 + 2];

            int ab = get_midpoint(g, &em, a, b);
            int bc = get_midpoint(g, &em, b, c);
            int ca = get_midpoint(g, &em, c, a);

            g->tris[t * 3] = a;
            g->tris[t * 3 + 1] = ab;
            g->tris[t * 3 + 2] = ca;

            geo_add_tri(g, ab, b, bc);
            geo_add_tri(g, ca, bc, c);
            geo_add_tri(g, ab, bc, ca);
        }

        edgemap_free(&em);
    }
}

// ---- Build the dual (hex/pentagon cells) ----

typedef struct {
    int* lists;
    int* offsets;
    int* counts;
} VertexTriMap;

static VertexTriMap build_vertex_tri_map(const GeoMesh* g) {
    VertexTriMap m;
    m.counts = (int*)calloc(g->vert_count, sizeof(int));
    m.offsets = (int*)calloc(g->vert_count, sizeof(int));

    for (int t = 0; t < g->tri_count; t++) {
        m.counts[g->tris[t * 3]]++;
        m.counts[g->tris[t * 3 + 1]]++;
        m.counts[g->tris[t * 3 + 2]]++;
    }

    int total = 0;
    for (int i = 0; i < g->vert_count; i++) {
        m.offsets[i] = total;
        total += m.counts[i];
    }

    m.lists = (int*)calloc(total, sizeof(int));
    int* cursor = (int*)calloc(g->vert_count, sizeof(int));

    for (int t = 0; t < g->tri_count; t++) {
        for (int k = 0; k < 3; k++) {
            int vi = g->tris[t * 3 + k];
            m.lists[m.offsets[vi] + cursor[vi]] = t;
            cursor[vi]++;
        }
    }

    free(cursor);
    return m;
}

static void free_vertex_tri_map(VertexTriMap* m) {
    free(m->lists);
    free(m->offsets);
    free(m->counts);
}

static void order_ring(const GeoMesh* g, int vertex_idx,
                       const int* tri_list, int tri_count,
                       int* ordered) {
    if (tri_count == 0) return;

    bool* used = (bool*)calloc(tri_count, sizeof(bool));
    ordered[0] = tri_list[0];
    used[0] = true;

    for (int i = 1; i < tri_count; i++) {
        int prev_tri = ordered[i - 1];
        int pa = g->tris[prev_tri * 3];
        int pb = g->tris[prev_tri * 3 + 1];
        int pc = g->tris[prev_tri * 3 + 2];

        int other[2];
        if (pa == vertex_idx) { other[0] = pb; other[1] = pc; }
        else if (pb == vertex_idx) { other[0] = pc; other[1] = pa; }
        else { other[0] = pa; other[1] = pb; }

        int shared_v = other[1];
        bool found = false;
        for (int j = 0; j < tri_count; j++) {
            if (used[j]) continue;
            int ta = g->tris[tri_list[j] * 3];
            int tb = g->tris[tri_list[j] * 3 + 1];
            int tc = g->tris[tri_list[j] * 3 + 2];
            if ((ta == vertex_idx || tb == vertex_idx || tc == vertex_idx) &&
                (ta == shared_v || tb == shared_v || tc == shared_v)) {
                ordered[i] = tri_list[j];
                used[j] = true;
                found = true;
                break;
            }
        }

        if (!found) {
            for (int j = 0; j < tri_count; j++) {
                if (!used[j]) {
                    ordered[i] = tri_list[j];
                    used[j] = true;
                    break;
                }
            }
        }
    }

    free(used);
}

static HMM_Vec3 tri_centroid(const GeoMesh* g, int tri_idx) {
    int a = g->tris[tri_idx * 3];
    int b = g->tris[tri_idx * 3 + 1];
    int c = g->tris[tri_idx * 3 + 2];
    HMM_Vec3 center = vec3_scale(
        vec3_add(vec3_add(g->verts[a], g->verts[b]), g->verts[c]),
        1.0f / 3.0f);
    return vec3_normalize(center);
}

// ---- Neighbor building ----
// Two dual cells (geodesic vertices) are neighbors if they share a geodesic edge.
// We extract all unique edges from the geodesic triangles and use them to populate
// the neighbor arrays. The neighbor order matches the cell vertex order: neighbor[k]
// shares the edge between vertices[k] and vertices[(k+1)%n].

static void build_neighbors(Planet* planet, const GeoMesh* geo) {
    // First, initialize all neighbors to -1
    for (int i = 0; i < planet->cell_count; i++) {
        for (int k = 0; k < MAX_CELL_VERTICES; k++) {
            planet->cells[i].neighbors[k] = -1;
        }
        planet->cells[i].neighbor_count = 0;
    }

    // Collect unique edges using EdgeMap
    // Each geodesic edge connects two vertices = two dual cells that are neighbors
    int edge_cap = geo->tri_count * 4;
    EdgeMap seen;
    edgemap_init(&seen, edge_cap);

    // For each cell, we need to know which other cells share each of its edges.
    // A dual cell edge (between vertices[k] and vertices[k+1]) corresponds to the
    // geodesic edge between the two geodesic vertices whose triangle rings produce
    // those centroids. More directly: two geodesic vertices connected by a geodesic
    // edge are always neighbors in the dual.

    // Collect all geodesic edges and mark neighbor pairs
    for (int t = 0; t < geo->tri_count; t++) {
        int v[3] = { geo->tris[t * 3], geo->tris[t * 3 + 1], geo->tris[t * 3 + 2] };
        for (int e = 0; e < 3; e++) {
            int a = v[e];
            int b = v[(e + 1) % 3];
            uint64_t key = edge_key(a, b);
            if (edgemap_get(&seen, key) >= 0) continue;
            edgemap_put(&seen, key, 1);

            // a and b are geodesic vertex indices = cell indices
            // They are neighbors. Find the correct slot in each cell's neighbor array.
            HexCell* ca = &planet->cells[a];
            HexCell* cb = &planet->cells[b];

            // For cell a: find which edge slot corresponds to neighbor b
            // The shared dual edge is between two consecutive cell vertices that
            // are centroids of triangles containing both a and b.
            // Simple approach: just append to neighbor list
            if (ca->neighbor_count < ca->vertex_count) {
                ca->neighbors[ca->neighbor_count++] = b;
            }
            if (cb->neighbor_count < cb->vertex_count) {
                cb->neighbors[cb->neighbor_count++] = a;
            }
        }
    }

    edgemap_free(&seen);

    // Now reorder neighbors to match vertex edge order.
    // neighbor[k] should be the cell that shares the edge between vertices[k] and vertices[(k+1)%n].
    // The shared edge's midpoint should be closest to the midpoint of (vertices[k], vertices[(k+1)%n]).
    for (int i = 0; i < planet->cell_count; i++) {
        HexCell* cell = &planet->cells[i];
        int n = cell->vertex_count;
        int* reordered = (int*)calloc(n, sizeof(int));
        for (int k = 0; k < n; k++) reordered[k] = -1;

        for (int k = 0; k < n; k++) {
            // Midpoint of edge k
            int next = (k + 1) % n;
            HMM_Vec3 mid = vec3_scale(vec3_add(cell->vertices[k], cell->vertices[next]), 0.5f);

            // Find the neighbor whose center is closest to this edge midpoint
            float best_dist = 1e30f;
            int best_nb = -1;
            for (int j = 0; j < cell->neighbor_count; j++) {
                int nb = cell->neighbors[j];
                // Check this neighbor hasn't already been assigned
                bool already_used = false;
                for (int kk = 0; kk < k; kk++) {
                    if (reordered[kk] == nb) { already_used = true; break; }
                }
                if (already_used) continue;

                HMM_Vec3 diff = vec3_sub(planet->cells[nb].center, mid);
                float dist = vec3_dot(diff, diff);
                if (dist < best_dist) {
                    best_dist = dist;
                    best_nb = nb;
                }
            }
            reordered[k] = best_nb;
        }

        memcpy(cell->neighbors, reordered, n * sizeof(int));
        free(reordered);
    }
}

// ---- Terrain generation with FastNoiseLite ----

static void generate_terrain(Planet* planet) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 4;
    noise.frequency = 2.0f;
    noise.seed = 42;

    planet->sea_level = 24;
    planet->layer_thickness = 1.0f;  // Each layer is 1 meter thick

    for (int i = 0; i < planet->cell_count; i++) {
        HexCell* cell = &planet->cells[i];
        HMM_Vec3 c = cell->center;

        // Sample 3D noise at the cell center on the unit sphere
        float n = fnlGetNoise3D(&noise, c.X * 3.0f, c.Y * 3.0f, c.Z * 3.0f);
        // n is in [-1, 1], map to height range [8..48]
        int height = 8 + (int)((n + 1.0f) * 0.5f * 40.0f);
        if (height < 1) height = 1;
        if (height >= MAX_VOXEL_HEIGHT) height = MAX_VOXEL_HEIGHT - 1;

        cell->terrain_height = height;

        // Fill voxel column
        memset(cell->voxels, VOXEL_AIR, MAX_VOXEL_HEIGHT);

        for (int layer = 0; layer <= height; layer++) {
            int rel = layer - planet->sea_level;  // Relative to sea level

            if (rel < -3) {
                cell->voxels[layer] = VOXEL_STONE;
            } else if (rel < 0) {
                cell->voxels[layer] = VOXEL_DIRT;
            } else if (rel == 0) {
                // At sea level: sand beach if terrain is near sea level
                if (height <= planet->sea_level + 1) {
                    cell->voxels[layer] = VOXEL_SAND;
                } else {
                    cell->voxels[layer] = VOXEL_DIRT;
                }
            } else if (rel == 1 && height <= planet->sea_level + 2) {
                cell->voxels[layer] = VOXEL_SAND;
            } else if (rel <= 6) {
                cell->voxels[layer] = VOXEL_GRASS;
            } else if (rel <= 12) {
                cell->voxels[layer] = VOXEL_STONE;
            } else {
                cell->voxels[layer] = VOXEL_ICE;
            }
        }

        // Fill water above terrain up to sea level
        if (height < planet->sea_level) {
            for (int layer = height + 1; layer <= planet->sea_level; layer++) {
                cell->voxels[layer] = VOXEL_WATER;
            }
            // Update terrain_height to include water surface for rendering
            // but keep terrain_height as the solid ground for collision
        }
    }
}

// ---- Color lookup for voxel types ----

HMM_Vec3 voxel_color(VoxelType type, int layer, int sea_level) {
    (void)layer;
    (void)sea_level;
    switch (type) {
        case VOXEL_WATER:
            return (HMM_Vec3){{0.15f, 0.35f, 0.65f}};
        case VOXEL_SAND:
            return (HMM_Vec3){{0.76f, 0.70f, 0.40f}};
        case VOXEL_DIRT:
            return (HMM_Vec3){{0.45f, 0.32f, 0.18f}};
        case VOXEL_GRASS:
            return (HMM_Vec3){{0.20f, 0.60f, 0.15f}};
        case VOXEL_STONE:
            return (HMM_Vec3){{0.50f, 0.50f, 0.50f}};
        case VOXEL_ICE:
            return (HMM_Vec3){{0.90f, 0.95f, 1.00f}};
        default:
            return (HMM_Vec3){{1.0f, 0.0f, 1.0f}};  // Magenta = error
    }
}

// ---- Init / Destroy ----

void planet_init(Planet* planet, int subdivision) {
    planet->subdivision = subdivision;
    planet->radius = 796000.0f;
    planet->dirty = false;
    planet->last_found_cell = -1;

    int est_verts = 10 * subdivision * subdivision + 100;
    int est_tris = 20 * subdivision * subdivision + 100;

    GeoMesh geo;
    geo_init(&geo, est_verts, est_tris);
    build_geodesic(&geo, subdivision);

    VertexTriMap vtmap = build_vertex_tri_map(&geo);

    planet->cell_count = geo.vert_count;
    planet->cells = (HexCell*)calloc(planet->cell_count, sizeof(HexCell));

    int* ordered_ring = (int*)calloc(12, sizeof(int));

    for (int vi = 0; vi < geo.vert_count; vi++) {
        HexCell* cell = &planet->cells[vi];
        int ntris = vtmap.counts[vi];
        const int* tri_list = &vtmap.lists[vtmap.offsets[vi]];

        order_ring(&geo, vi, tri_list, ntris, ordered_ring);

        cell->center = geo.verts[vi];
        cell->normal = vec3_normalize(geo.verts[vi]);
        cell->vertex_count = ntris;

        for (int t = 0; t < ntris && t < MAX_CELL_VERTICES; t++) {
            cell->vertices[t] = tri_centroid(&geo, ordered_ring[t]);
        }

        // Verify winding: CCW when viewed from outside
        if (ntris >= 3) {
            HMM_Vec3 edge1 = vec3_sub(cell->vertices[1], cell->vertices[0]);
            HMM_Vec3 edge2 = vec3_sub(cell->vertices[2], cell->vertices[0]);
            HMM_Vec3 face_normal = vec3_cross(edge1, edge2);
            if (vec3_dot(face_normal, cell->normal) < 0) {
                for (int j = 0; j < ntris / 2; j++) {
                    HMM_Vec3 tmp = cell->vertices[j];
                    cell->vertices[j] = cell->vertices[ntris - 1 - j];
                    cell->vertices[ntris - 1 - j] = tmp;
                }
            }
        }
    }

    free(ordered_ring);

    // Build neighbor adjacency from geodesic edges
    build_neighbors(planet, &geo);

    free_vertex_tri_map(&vtmap);
    geo_free(&geo);

    // Generate terrain with noise
    generate_terrain(planet);

    planet->dirty = true;
}

void planet_destroy(Planet* planet) {
    free(planet->cells);
    planet->cells = NULL;
    planet->cell_count = 0;
}

// ---- Spatial queries ----

// Brute-force O(N) cell search (used as fallback)
static int planet_find_cell_brute(const Planet* planet, HMM_Vec3 pos) {
    HMM_Vec3 dir = vec3_normalize(pos);
    int best = -1;
    float best_dot = -2.0f;

    for (int i = 0; i < planet->cell_count; i++) {
        float d = vec3_dot(dir, planet->cells[i].normal);
        if (d > best_dot) {
            best_dot = d;
            best = i;
        }
    }
    return best;
}

int planet_find_cell_from_hint(const Planet* planet, HMM_Vec3 pos, int hint) {
    if (hint < 0 || hint >= planet->cell_count) {
        return planet_find_cell_brute(planet, pos);
    }

    HMM_Vec3 dir = vec3_normalize(pos);
    int current = hint;

    // Greedy walk: move to whichever neighbor has higher dot product
    for (int iter = 0; iter < planet->cell_count; iter++) {
        float current_dot = vec3_dot(dir, planet->cells[current].normal);
        int best_neighbor = -1;
        float best_dot = current_dot;

        const HexCell* cell = &planet->cells[current];
        for (int i = 0; i < cell->neighbor_count; i++) {
            int nb = cell->neighbors[i];
            if (nb < 0) continue;
            float d = vec3_dot(dir, planet->cells[nb].normal);
            if (d > best_dot) {
                best_dot = d;
                best_neighbor = nb;
            }
        }

        if (best_neighbor < 0) {
            // 1-ring converged. Do a 2-ring check to handle near-Voronoi edges.
            for (int i = 0; i < cell->neighbor_count; i++) {
                int nb = cell->neighbors[i];
                if (nb < 0) continue;
                const HexCell* nb_cell = &planet->cells[nb];
                for (int j = 0; j < nb_cell->neighbor_count; j++) {
                    int nb2 = nb_cell->neighbors[j];
                    if (nb2 < 0 || nb2 == current) continue;
                    float d = vec3_dot(dir, planet->cells[nb2].normal);
                    if (d > best_dot) {
                        best_dot = d;
                        best_neighbor = nb2;
                    }
                }
            }
            if (best_neighbor < 0) {
                return current;
            }
            // Found a better cell in 2-ring, continue walking from there
        }
        current = best_neighbor;
    }

    return current;
}

int planet_find_cell(Planet* planet, HMM_Vec3 pos) {
    int result = planet_find_cell_from_hint(planet, pos, planet->last_found_cell);
    planet->last_found_cell = result;
    return result;
}

int planet_cell_distance(const Planet* planet, int a, int b, int max_dist) {
    if (a == b) return 0;
    if (a < 0 || b < 0 || a >= planet->cell_count || b >= planet->cell_count) return -1;

    // BFS from a to b, bounded by max_dist
    int* queue = (int*)malloc(planet->cell_count * sizeof(int));
    int* dist = (int*)malloc(planet->cell_count * sizeof(int));
    for (int i = 0; i < planet->cell_count; i++) dist[i] = -1;

    dist[a] = 0;
    queue[0] = a;
    int head = 0, tail = 1;

    while (head < tail) {
        int cur = queue[head++];
        int cur_dist = dist[cur];
        if (cur_dist >= max_dist) continue;

        const HexCell* cell = &planet->cells[cur];
        for (int i = 0; i < cell->neighbor_count; i++) {
            int nb = cell->neighbors[i];
            if (nb < 0 || dist[nb] >= 0) continue;
            dist[nb] = cur_dist + 1;
            if (nb == b) {
                int result = dist[nb];
                free(queue);
                free(dist);
                return result;
            }
            queue[tail++] = nb;
        }
    }

    free(queue);
    free(dist);
    return -1;
}

float planet_terrain_radius(const Planet* planet, int cell_index) {
    if (cell_index < 0 || cell_index >= planet->cell_count) return planet->radius;
    const HexCell* cell = &planet->cells[cell_index];
    // Top of the terrain (or sea level if underwater)
    int top = cell->terrain_height;
    if (top < planet->sea_level) top = planet->sea_level;
    return planet->radius + (top + 1) * planet->layer_thickness;
}

// ---- Raycast ----

static bool ray_tri_intersect(HMM_Vec3 orig, HMM_Vec3 dir,
                               HMM_Vec3 v0, HMM_Vec3 v1, HMM_Vec3 v2,
                               float* t_out) {
    HMM_Vec3 e1 = vec3_sub(v1, v0);
    HMM_Vec3 e2 = vec3_sub(v2, v0);
    HMM_Vec3 h = vec3_cross(dir, e2);
    float a = vec3_dot(e1, h);
    if (a > -1e-6f && a < 1e-6f) return false;

    float f = 1.0f / a;
    HMM_Vec3 s = vec3_sub(orig, v0);
    float u = f * vec3_dot(s, h);
    if (u < 0.0f || u > 1.0f) return false;

    HMM_Vec3 q = vec3_cross(s, e1);
    float v = f * vec3_dot(dir, q);
    if (v < 0.0f || u + v > 1.0f) return false;

    float t = f * vec3_dot(e2, q);
    if (t > 1e-6f) {
        *t_out = t;
        return true;
    }
    return false;
}

int planet_raycast(const Planet* planet, HMM_Vec3 ray_origin, HMM_Vec3 ray_dir, float max_dist) {
    // Simplified raycast: find the cell whose top surface is hit
    // For each cell, test ray against the top face triangle fan at terrain radius
    int best_cell = -1;
    float best_t = max_dist;

    for (int i = 0; i < planet->cell_count; i++) {
        const HexCell* cell = &planet->cells[i];
        int top_layer = cell->terrain_height;
        if (top_layer < planet->sea_level) top_layer = planet->sea_level;
        if (cell->voxels[top_layer] == VOXEL_AIR) continue;

        float r_top = planet->radius + (top_layer + 1) * planet->layer_thickness;

        // Quick distance check
        HMM_Vec3 scaled_center = vec3_scale(cell->center, r_top);
        HMM_Vec3 oc = vec3_sub(scaled_center, ray_origin);
        float proj = vec3_dot(oc, ray_dir);
        if (proj < 0 || proj > best_t) continue;

        // Test ray against top face triangle fan
        HMM_Vec3 center_pos = vec3_scale(cell->center, r_top);
        for (int v = 0; v < cell->vertex_count; v++) {
            int next = (v + 1) % cell->vertex_count;
            HMM_Vec3 v0 = center_pos;
            HMM_Vec3 v1 = vec3_scale(cell->vertices[next], r_top);
            HMM_Vec3 v2 = vec3_scale(cell->vertices[v], r_top);
            float t;
            if (ray_tri_intersect(ray_origin, ray_dir, v0, v1, v2, &t)) {
                if (t < best_t) {
                    best_t = t;
                    best_cell = i;
                }
            }
        }
    }

    return best_cell;
}

// ---- Block interactions ----

bool planet_break_cell(Planet* planet, int cell_index) {
    if (cell_index < 0 || cell_index >= planet->cell_count) return false;
    HexCell* cell = &planet->cells[cell_index];
    if (cell->terrain_height < 0) return false;

    // Remove the top layer
    int top = cell->terrain_height;
    cell->voxels[top] = VOXEL_AIR;

    // Find new terrain height
    while (cell->terrain_height > 0 && cell->voxels[cell->terrain_height] == VOXEL_AIR) {
        cell->terrain_height--;
    }
    if (cell->voxels[cell->terrain_height] == VOXEL_AIR) {
        cell->terrain_height = -1;  // Column is completely empty
    }

    planet->dirty = true;
    return true;
}

int planet_place_cell(Planet* planet, int adjacent_to, HMM_Vec3 hit_point) {
    (void)hit_point;
    if (adjacent_to < 0 || adjacent_to >= planet->cell_count) return -1;
    HexCell* cell = &planet->cells[adjacent_to];

    // Place a stone block on top of the terrain
    int new_layer = cell->terrain_height + 1;
    if (new_layer >= MAX_VOXEL_HEIGHT) return -1;

    cell->voxels[new_layer] = VOXEL_STONE;
    cell->terrain_height = new_layer;
    planet->dirty = true;
    return adjacent_to;
}

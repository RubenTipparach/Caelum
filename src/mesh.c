#include "mesh.h"
#include "math_utils.h"
#include <stdlib.h>
#include <string.h>

// Forward declaration from planet.c
extern HMM_Vec3 voxel_color(VoxelType type, int layer, int sea_level);

// Helper: emit a single triangle (3 vertices) into the mesh
static void emit_tri(PlanetMesh* mesh, int* count, int* max_verts,
                     HMM_Vec3 p0, HMM_Vec3 p1, HMM_Vec3 p2,
                     HMM_Vec3 normal, HMM_Vec3 color) {
    if (*count + 3 > *max_verts) {
        *max_verts *= 2;
        mesh->vertices = (Vertex*)realloc(mesh->vertices, *max_verts * sizeof(Vertex));
    }
    Vertex* vp = &mesh->vertices[*count];

    vp[0].pos[0] = p0.X; vp[0].pos[1] = p0.Y; vp[0].pos[2] = p0.Z;
    vp[0].normal[0] = normal.X; vp[0].normal[1] = normal.Y; vp[0].normal[2] = normal.Z;
    vp[0].color[0] = color.X; vp[0].color[1] = color.Y; vp[0].color[2] = color.Z;

    vp[1].pos[0] = p1.X; vp[1].pos[1] = p1.Y; vp[1].pos[2] = p1.Z;
    vp[1].normal[0] = normal.X; vp[1].normal[1] = normal.Y; vp[1].normal[2] = normal.Z;
    vp[1].color[0] = color.X; vp[1].color[1] = color.Y; vp[1].color[2] = color.Z;

    vp[2].pos[0] = p2.X; vp[2].pos[1] = p2.Y; vp[2].pos[2] = p2.Z;
    vp[2].normal[0] = normal.X; vp[2].normal[1] = normal.Y; vp[2].normal[2] = normal.Z;
    vp[2].color[0] = color.X; vp[2].color[1] = color.Y; vp[2].color[2] = color.Z;

    *count += 3;
}

// Generate mesh geometry for a single cell.
// visible_set: if non-NULL, only treat neighbors in the set as potentially solid.
//              Neighbors outside the set are treated as AIR (renders boundary walls).
// visible_set is a boolean array of size planet->cell_count.
static void mesh_generate_cell(const Planet* planet, int cell_idx,
                                const bool* visible_set,
                                PlanetMesh* mesh, int* count, int* max_verts) {
    const HexCell* cell = &planet->cells[cell_idx];
    int n = cell->vertex_count;

    for (int layer = 0; layer < MAX_VOXEL_HEIGHT; layer++) {
        VoxelType vtype = (VoxelType)cell->voxels[layer];
        if (vtype == VOXEL_AIR) continue;

        HMM_Vec3 color = voxel_color(vtype, layer, planet->sea_level);

        float r_bot = planet->radius + layer * planet->layer_thickness;
        float r_top = r_bot + planet->layer_thickness;

        // ---- Top face: render if layer above is AIR ----
        bool render_top = (layer == MAX_VOXEL_HEIGHT - 1) ||
                          (cell->voxels[layer + 1] == VOXEL_AIR);
        if (render_top) {
            HMM_Vec3 center_top = vec3_scale(cell->center, r_top);
            HMM_Vec3 top_normal = cell->normal;

            for (int k = 0; k < n; k++) {
                int next = (k + 1) % n;
                HMM_Vec3 v0 = center_top;
                HMM_Vec3 v1 = vec3_scale(cell->vertices[next], r_top);
                HMM_Vec3 v2 = vec3_scale(cell->vertices[k], r_top);
                emit_tri(mesh, count, max_verts, v0, v1, v2, top_normal, color);
            }
        }

        // ---- Bottom face: render if layer below is AIR (or layer 0) ----
        bool render_bot = (layer == 0) ||
                          (cell->voxels[layer - 1] == VOXEL_AIR);
        if (render_bot) {
            HMM_Vec3 center_bot = vec3_scale(cell->center, r_bot);
            HMM_Vec3 bot_normal = vec3_scale(cell->normal, -1.0f);

            for (int k = 0; k < n; k++) {
                int next = (k + 1) % n;
                HMM_Vec3 v0 = center_bot;
                HMM_Vec3 v1 = vec3_scale(cell->vertices[k], r_bot);
                HMM_Vec3 v2 = vec3_scale(cell->vertices[next], r_bot);
                emit_tri(mesh, count, max_verts, v0, v1, v2, bot_normal, color);
            }
        }

        // ---- Side walls: render if neighbor at same layer is AIR ----
        for (int k = 0; k < n; k++) {
            int nb_idx = cell->neighbors[k];
            bool render_side = false;

            if (nb_idx < 0) {
                render_side = true;
            } else if (visible_set && !visible_set[nb_idx]) {
                // Neighbor is outside visible set — treat as AIR (boundary wall)
                render_side = true;
            } else {
                render_side = (planet->cells[nb_idx].voxels[layer] == VOXEL_AIR);
            }

            if (render_side) {
                int next = (k + 1) % n;
                HMM_Vec3 bl = vec3_scale(cell->vertices[k], r_bot);
                HMM_Vec3 br = vec3_scale(cell->vertices[next], r_bot);
                HMM_Vec3 tl = vec3_scale(cell->vertices[k], r_top);
                HMM_Vec3 tr = vec3_scale(cell->vertices[next], r_top);

                HMM_Vec3 edge_mid = vec3_scale(vec3_add(cell->vertices[k], cell->vertices[next]), 0.5f);
                HMM_Vec3 side_normal = vec3_normalize(vec3_sub(edge_mid, cell->center));

                emit_tri(mesh, count, max_verts, bl, tr, br, side_normal, color);
                emit_tri(mesh, count, max_verts, bl, tl, tr, side_normal, color);
            }
        }
    }
}

PlanetMesh planet_mesh_generate(const Planet* planet) {
    PlanetMesh mesh;
    int max_verts = 4000000;
    mesh.vertices = (Vertex*)malloc(max_verts * sizeof(Vertex));
    mesh.vertex_count = 0;
    int count = 0;

    for (int i = 0; i < planet->cell_count; i++) {
        mesh_generate_cell(planet, i, NULL, &mesh, &count, &max_verts);
    }

    mesh.vertex_count = count;
    if (count > 0) {
        mesh.vertices = (Vertex*)realloc(mesh.vertices, count * sizeof(Vertex));
    }
    return mesh;
}

PlanetMesh planet_mesh_generate_nearby(const Planet* planet, int center_cell, int render_distance) {
    if (center_cell < 0 || center_cell >= planet->cell_count) {
        return (PlanetMesh){ NULL, 0 };
    }

    // BFS to collect visible cells within render_distance hops
    bool* visible = (bool*)calloc(planet->cell_count, sizeof(bool));
    int* queue = (int*)malloc(planet->cell_count * sizeof(int));
    int* dist = (int*)malloc(planet->cell_count * sizeof(int));
    for (int i = 0; i < planet->cell_count; i++) dist[i] = -1;

    visible[center_cell] = true;
    dist[center_cell] = 0;
    queue[0] = center_cell;
    int head = 0, tail = 1;
    int visible_count = 1;

    while (head < tail) {
        int cur = queue[head++];
        int cur_dist = dist[cur];
        if (cur_dist >= render_distance) continue;

        const HexCell* cell = &planet->cells[cur];
        for (int i = 0; i < cell->neighbor_count; i++) {
            int nb = cell->neighbors[i];
            if (nb < 0 || dist[nb] >= 0) continue;
            dist[nb] = cur_dist + 1;
            visible[nb] = true;
            queue[tail++] = nb;
            visible_count++;
        }
    }

    // Generate mesh only for visible cells
    PlanetMesh mesh;
    int max_verts = visible_count * 200;  // ~200 verts per cell estimate
    if (max_verts < 10000) max_verts = 10000;
    mesh.vertices = (Vertex*)malloc(max_verts * sizeof(Vertex));
    mesh.vertex_count = 0;
    int count = 0;

    for (int i = 0; i < tail; i++) {
        mesh_generate_cell(planet, queue[i], visible, &mesh, &count, &max_verts);
    }

    mesh.vertex_count = count;
    if (count > 0) {
        mesh.vertices = (Vertex*)realloc(mesh.vertices, count * sizeof(Vertex));
    } else {
        free(mesh.vertices);
        mesh.vertices = NULL;
    }

    free(visible);
    free(queue);
    free(dist);
    return mesh;
}

void planet_mesh_destroy(PlanetMesh* mesh) {
    free(mesh->vertices);
    mesh->vertices = NULL;
    mesh->vertex_count = 0;
}

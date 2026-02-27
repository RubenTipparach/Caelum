#ifndef MESH_H
#define MESH_H

#include "planet.h"

typedef struct Vertex {
    float pos[3];
    float normal[3];
    float color[3];
} Vertex;

typedef struct PlanetMesh {
    Vertex* vertices;
    int vertex_count;
} PlanetMesh;

// Generate a triangle mesh from ALL planet hex cells.
// Caller must free the mesh with planet_mesh_destroy().
PlanetMesh planet_mesh_generate(const Planet* planet);

// Generate mesh only for cells within render_distance hops of center_cell.
// Uses BFS through neighbor graph. Boundary neighbors treated as AIR.
PlanetMesh planet_mesh_generate_nearby(const Planet* planet, int center_cell, int render_distance);

void planet_mesh_destroy(PlanetMesh* mesh);

#endif

#ifndef PLANET_H
#define PLANET_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"

#define MAX_CELL_VERTICES 6
#define MAX_VOXEL_HEIGHT 64

typedef enum {
    VOXEL_AIR = 0,
    VOXEL_WATER,
    VOXEL_SAND,
    VOXEL_DIRT,
    VOXEL_GRASS,
    VOXEL_STONE,
    VOXEL_ICE,
} VoxelType;

typedef struct HexCell {
    HMM_Vec3 center;                         // Position on unit sphere surface
    HMM_Vec3 normal;                         // Surface normal (= normalized center)
    HMM_Vec3 vertices[MAX_CELL_VERTICES];    // Corner vertices on unit sphere
    int vertex_count;                         // 6 for hex, 5 for pentagon
    int neighbors[MAX_CELL_VERTICES];         // Adjacent cell indices (-1 = none)
    int neighbor_count;
    uint8_t voxels[MAX_VOXEL_HEIGHT];        // VoxelType per radial layer (0=innermost)
    int terrain_height;                       // Index of topmost solid non-water layer
} HexCell;

typedef struct Planet {
    HexCell* cells;
    int cell_count;
    float radius;            // Base radius
    int subdivision;         // The 'm' in GP(m,0)
    int sea_level;           // Layer index for water surface
    float layer_thickness;   // Radial thickness per voxel layer (meters)
    bool dirty;              // Mesh needs regeneration
    int last_found_cell;     // Cached hint for greedy neighbor walk
} Planet;

void planet_init(Planet* planet, int subdivision);
void planet_destroy(Planet* planet);

// Find the cell index nearest to a world position (uses cached hint for speed)
int planet_find_cell(Planet* planet, HMM_Vec3 pos);

// Find the cell index using a greedy walk from a hint cell
int planet_find_cell_from_hint(const Planet* planet, HMM_Vec3 pos, int hint);

// BFS distance between two cells (-1 if unreachable within max_dist)
int planet_cell_distance(const Planet* planet, int a, int b, int max_dist);

// Returns the world radius at the top of the terrain column for a cell
float planet_terrain_radius(const Planet* planet, int cell_index);

// Raycast against planet cells. Returns cell index or -1 if no hit.
int planet_raycast(const Planet* planet, HMM_Vec3 ray_origin, HMM_Vec3 ray_dir, float max_dist);

// Break (remove) a voxel. Returns true if state changed.
bool planet_break_cell(Planet* planet, int cell_index);

// Place a voxel on top of existing terrain. Returns the placed cell index or -1.
int planet_place_cell(Planet* planet, int adjacent_to, HMM_Vec3 hit_point);

#endif

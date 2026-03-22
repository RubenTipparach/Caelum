#ifndef AI_PATHFIND_H
#define AI_PATHFIND_H

#include <stdbool.h>
#include "hex_terrain.h"

// --- Hex Grid A* Pathfinding ---
// Finds a walkable path between two hex grid positions.
// Uses hex_terrain for ground/obstacle queries.

#define AI_PATH_MAX_STEPS 128

typedef struct {
    int q, r;       // hex grid coordinates for each step
} AiPathStep;

typedef struct {
    AiPathStep steps[AI_PATH_MAX_STEPS];
    int step_count;
    int current_step;   // index of next step to walk to
    bool valid;
} AiPath;

// Find a path from (start_q, start_r) to (end_q, end_r).
// Returns true if a path was found. Uses A* with hex neighbors.
// The agent needs 2 layers of headroom to walk through.
bool ai_pathfind(AiPath* path, const HexTerrain* ht,
                 int start_q, int start_r,
                 int end_q, int end_r);

// Get the next step position. Returns false if path is complete.
bool ai_path_next(const AiPath* path, int* out_q, int* out_r);

// Advance to the next step.
void ai_path_advance(AiPath* path);

// Check if path is complete.
bool ai_path_done(const AiPath* path);

#endif

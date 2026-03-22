#ifndef AI_SENSORS_H
#define AI_SENSORS_H

#include "hex_terrain.h"
#include "HandmadeMath.h"

// --- AI Sensor System ---
// Provides terrain awareness for AI agents.
// Sensors scan the hex terrain and produce structured text descriptions
// that can be injected into LLM context for informed decision-making.

// Maximum size of a sensor report text buffer
#define AI_SENSOR_REPORT_MAX 2048

// A single hex cell observation
typedef struct {
    int q, r;                // global hex coords
    int ground_layer;        // topmost solid layer
    float ground_r;          // ground height (radius)
    uint8_t surface_type;    // VoxelType of topmost block
    float height_diff;       // relative to agent's ground height
    float distance;          // hex distance from agent
    bool walkable;           // has ground + 2 layers headroom
} AiHexObs;

// Scan result: nearby terrain around a position
#define AI_SCAN_MAX_CELLS 128
typedef struct {
    AiHexObs cells[AI_SCAN_MAX_CELLS];
    int count;
    int agent_q, agent_r;        // agent's current hex
    int agent_ground_layer;      // agent's ground layer
    float agent_ground_r;        // agent's ground radius
} AiScanResult;

// Scan terrain in a radius around global hex coords.
// radius: in hex cells (e.g. 5 = scan 5 cells out = ~61 cells)
void ai_sensors_scan(const HexTerrain* ht, int center_q, int center_r,
                     int radius, AiScanResult* out);

// Generate a text report from a scan result for LLM context injection.
// Includes: surface composition, elevation profile, walkable paths, obstacles.
// Returns the number of chars written.
int ai_sensors_report(const AiScanResult* scan, char* buf, int buf_size);

// Scan a vertical column at (q, r) and describe the voxel stack.
// Useful for understanding what's below/above a specific cell.
int ai_sensors_column_report(const HexTerrain* ht, int q, int r,
                             char* buf, int buf_size);

// Get a directional description: what the agent sees in each of 6 hex directions.
// Reports terrain height changes and block types in each direction up to 'depth' cells.
int ai_sensors_directional_report(const HexTerrain* ht,
                                  int center_q, int center_r,
                                  int depth, char* buf, int buf_size);

#endif

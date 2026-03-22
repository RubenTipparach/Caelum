#ifndef AI_ORDERS_H
#define AI_ORDERS_H

#include "ai_script.h"
#include "ai_memory.h"
#include "hex_terrain.h"
#include "HandmadeMath.h"

// --- AI Orders ---
// Parses natural language commands into executable scripts.
// Orders are saved to the agent's learning folder for future reference.
//
// Supported commands:
//   "follow me"              — follow the player continuously
//   "stay" / "stop"          — stop and idle
//   "go north 10"            — walk ~10m north (N/S/E/W/NE/NW/SE/SW)
//   "build wall north 5"     — build a 5-block wall going north
//   "build stairs north 3"   — build 3 steps going north
//   "clear area 3"           — break all blocks in a 3-cell radius above ground
//   "clear ahead 5"          — break blocks in the forward direction for 5 cells

typedef enum {
    ORDER_NONE = 0,
    ORDER_FOLLOW,
    ORDER_STAY,
    ORDER_GO,
    ORDER_BUILD_WALL,
    ORDER_BUILD_STAIRS,
    ORDER_CLEAR_AREA,
    ORDER_CLEAR_AHEAD,
} OrderType;

typedef struct {
    OrderType type;
    int direction;      // hex direction 0-5 (E=0, NE=1, NW=2, W=3, SW=4, SE=5), -1=forward
    int count;          // distance in cells, or block count, or radius
    char block[16];     // block type for building (default "stone")
} AiOrder;

// Parse a text command into an order. Returns true if recognized.
// Uses the agent's current forward direction to resolve "forward"/"ahead".
bool ai_order_parse(const char* text, AiOrder* out);

// Generate a script from an order, using the agent's current position.
// The script is ready to execute via ai_script_run().
bool ai_order_to_script(const AiOrder* order, const HexTerrain* ht,
                        int agent_q, int agent_r, int agent_ground_layer,
                        AiScript* out);

// Save an order (as its generated script) to the agent's learning folder.
void ai_order_save(const AiOrder* order, const AiScript* script,
                   AiMemory* mem);

#endif

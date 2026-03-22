#ifndef AI_SCRIPT_H
#define AI_SCRIPT_H

#include <stdbool.h>
#include <stdint.h>
#include "hex_terrain.h"

// --- AI Script System ---
// Scripts are generated from player-agent conversations in "engineer mode".
// They describe sequences of actions Walter can perform autonomously.
//
// Script format (JSON):
// {
//   "name": "build_staircase",
//   "description": "Build a 5-step staircase going NE",
//   "steps": [
//     { "action": "scan", "radius": 5 },
//     { "action": "move_to", "q": 10, "r": 5 },
//     { "action": "place", "q": 10, "r": 5, "layer": 5450, "block": "stone" },
//     { "action": "place", "q": 11, "r": 5, "layer": 5451, "block": "stone" },
//     { "action": "wait", "seconds": 1.0 },
//     { "action": "say", "text": "Done building!" },
//     { "action": "loop", "from": 2, "count": 5 },
//     { "action": "scan_and_report" }
//   ]
// }

#define AI_SCRIPT_MAX_STEPS  256
#define AI_SCRIPT_NAME_MAX   64
#define AI_SCRIPT_DESC_MAX   256
#define AI_SCRIPT_TEXT_MAX   256

typedef enum {
    SCRIPT_ACT_NONE = 0,
    SCRIPT_ACT_MOVE_TO,       // pathfind + walk to (q, r)
    SCRIPT_ACT_PLACE,         // place block at (q, r, layer) with type
    SCRIPT_ACT_BREAK,         // break block at (q, r, layer)
    SCRIPT_ACT_SAY,           // display text in chat
    SCRIPT_ACT_WAIT,          // pause for N seconds
    SCRIPT_ACT_SCAN,          // run sensor scan, store result
    SCRIPT_ACT_SCAN_REPORT,   // scan + inject report into LLM context
    SCRIPT_ACT_LOOK_AT,       // face toward (q, r)
    SCRIPT_ACT_JUMP,          // jump in place
    SCRIPT_ACT_LOOP,          // repeat from step 'from_step' N times
    SCRIPT_ACT_MOVE_REL,      // move relative to anchor (dq, dr)
    SCRIPT_ACT_PLACE_REL,     // place relative to anchor (dq, dr, dlayer)
    SCRIPT_ACT_BREAK_REL,     // break relative to anchor (dq, dr, dlayer)
    SCRIPT_ACT_SET_ANCHOR,    // set anchor to current position (or explicit q,r,layer)
} ScriptActionType;

typedef struct {
    ScriptActionType action;
    int q, r, layer;            // hex coords — absolute (move/place/break)
    int dq, dr, dlayer;         // hex coords — relative to anchor (move_rel/place_rel/break_rel)
    char block[16];             // block type name (place)
    char text[AI_SCRIPT_TEXT_MAX]; // say text
    float seconds;              // wait duration
    int radius;                 // scan radius
    int from_step;              // loop: go back to this step index
    int loop_count;             // loop: repeat N times
    int _loop_remaining;        // runtime: loops left (decremented during execution)
} AiScriptStep;

typedef struct {
    char name[AI_SCRIPT_NAME_MAX];
    char description[AI_SCRIPT_DESC_MAX];
    AiScriptStep steps[AI_SCRIPT_MAX_STEPS];
    int step_count;
} AiScript;

typedef enum {
    SCRIPT_IDLE,
    SCRIPT_RUNNING,
    SCRIPT_WAITING,     // paused on a wait step
    SCRIPT_MOVING,      // walking to a move_to target
    SCRIPT_BUILDING,    // facing target + cooldown before place/break
    SCRIPT_DONE,
    SCRIPT_ERROR,
} AiScriptState;

typedef struct {
    AiScript script;
    AiScriptState state;
    int current_step;
    float wait_timer;
    char error[128];

    // Anchor point for relative coordinates (set_anchor or auto-set on run)
    int anchor_q, anchor_r, anchor_layer;
    bool anchor_set;
} AiScriptRunner;

// Parse a script from JSON text. Returns true on success.
bool ai_script_parse(const char* json, AiScript* out);

// Load a script from a file. Returns true on success.
bool ai_script_load(const char* path, AiScript* out);

// Save a script to a file.
bool ai_script_save(const AiScript* script, const char* path);

// Initialize the script runner.
void ai_script_runner_init(AiScriptRunner* runner);

// Start executing a script.
void ai_script_run(AiScriptRunner* runner, const AiScript* script);

// Stop the current script.
void ai_script_stop(AiScriptRunner* runner);

// Update the script runner each frame. Returns the current step action
// if one needs external handling (move_to, place, break), or SCRIPT_ACT_NONE.
// The caller (ai_agent) handles pathfinding, block ops, etc.
ScriptActionType ai_script_update(AiScriptRunner* runner, float dt);

// Get the current step (for reading params). NULL if not running.
const AiScriptStep* ai_script_current_step(const AiScriptRunner* runner);

// Advance to the next step (called by agent after completing a move/place/break).
void ai_script_advance(AiScriptRunner* runner);

// Is the script done or idle?
bool ai_script_idle(const AiScriptRunner* runner);

#endif

#ifndef AI_ACTIONS_H
#define AI_ACTIONS_H

#include "ai_npc.h"
#include "hex_terrain.h"

// --- AI Action Executor ---
// Processes parsed AI actions and applies them to the hex terrain.
// Actions are queued and executed one at a time with rate limiting.

#define AI_ACTION_QUEUE_SIZE 64

typedef struct {
    AiAction queue[AI_ACTION_QUEUE_SIZE];
    int      queue_count;
    int      queue_head;        // next action to execute
    float    action_timer;      // countdown until next action fires
    float    action_interval;   // seconds between actions (default 0.5 = 2 blocks/sec)
    int      actions_executed;  // lifetime counter
    int      actions_failed;    // lifetime counter
} AiActionExecutor;

// Initialize the action executor.
void ai_actions_init(AiActionExecutor* exec);

// Enqueue actions from an AI response.
void ai_actions_enqueue(AiActionExecutor* exec, const AiAction* actions, int count);

// Update: execute queued actions at rate limit. Call each frame.
// ht: hex terrain to modify
// Returns true if an action was executed this frame.
bool ai_actions_update(AiActionExecutor* exec, HexTerrain* ht, float dt);

// Clear all pending actions.
void ai_actions_clear(AiActionExecutor* exec);

// Are there pending actions?
bool ai_actions_busy(const AiActionExecutor* exec);

#endif

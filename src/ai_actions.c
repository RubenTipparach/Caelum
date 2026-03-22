#include "ai_actions.h"
#include "planet.h"

#include <stdio.h>
#include <string.h>

// Map block name string to VoxelType
static uint8_t block_name_to_type(const char* name) {
    if (strcmp(name, "stone") == 0) return VOXEL_STONE;
    if (strcmp(name, "dirt") == 0)  return VOXEL_DIRT;
    if (strcmp(name, "grass") == 0) return VOXEL_GRASS;
    if (strcmp(name, "sand") == 0)  return VOXEL_SAND;
    if (strcmp(name, "water") == 0) return VOXEL_WATER;
    if (strcmp(name, "ice") == 0)   return VOXEL_ICE;
    if (strcmp(name, "snow") == 0)  return VOXEL_ICE; // snow uses ice texture
    if (strcmp(name, "torch") == 0) return VOXEL_TORCH;
    return VOXEL_STONE; // fallback
}

void ai_actions_init(AiActionExecutor* exec) {
    memset(exec, 0, sizeof(AiActionExecutor));
    exec->action_interval = 0.5f; // 2 blocks per second
}

void ai_actions_enqueue(AiActionExecutor* exec, const AiAction* actions, int count) {
    for (int i = 0; i < count && exec->queue_count < AI_ACTION_QUEUE_SIZE; i++) {
        exec->queue[exec->queue_count++] = actions[i];
    }
    if (count > 0) {
        printf("[AI-EXEC] Enqueued %d actions (%d total in queue)\n",
               count, exec->queue_count - exec->queue_head);
        fflush(stdout);
    }
}

bool ai_actions_update(AiActionExecutor* exec, HexTerrain* ht, float dt) {
    if (exec->queue_head >= exec->queue_count) return false;

    exec->action_timer -= dt;
    if (exec->action_timer > 0.0f) return false;

    exec->action_timer = exec->action_interval;

    AiAction* act = &exec->queue[exec->queue_head++];

    switch (act->type) {
        case AI_ACTION_PLACE_BLOCK: {
            // Construct a synthetic HexHitResult for placement
            HexHitResult hit;
            memset(&hit, 0, sizeof(hit));
            hit.valid = true;
            hit.place_gcol = act->q;
            hit.place_grow = act->r;
            hit.place_layer = act->layer;

            uint8_t voxel = block_name_to_type(act->block);
            if (hex_terrain_place(ht, &hit, voxel)) {
                exec->actions_executed++;
                printf("[AI-EXEC] Placed %s at q=%d r=%d layer=%d\n",
                       act->block, act->q, act->r, act->layer);
            } else {
                exec->actions_failed++;
                printf("[AI-EXEC] FAILED to place %s at q=%d r=%d layer=%d\n",
                       act->block, act->q, act->r, act->layer);
            }
            fflush(stdout);
            return true;
        }

        case AI_ACTION_BREAK_BLOCK: {
            HexHitResult hit;
            memset(&hit, 0, sizeof(hit));
            hit.valid = true;
            hit.gcol = act->q;
            hit.grow = act->r;
            hit.layer = act->layer;

            if (hex_terrain_break(ht, &hit)) {
                exec->actions_executed++;
                printf("[AI-EXEC] Broke block at q=%d r=%d layer=%d\n",
                       act->q, act->r, act->layer);
            } else {
                exec->actions_failed++;
                printf("[AI-EXEC] FAILED to break block at q=%d r=%d layer=%d\n",
                       act->q, act->r, act->layer);
            }
            fflush(stdout);
            return true;
        }

        case AI_ACTION_SAY:
            printf("[AI-NPC] \"%s\"\n", act->text);
            fflush(stdout);
            exec->actions_executed++;
            return true;

        case AI_ACTION_MOVE_TO:
            // TODO: pathfinding (Phase 2)
            printf("[AI-EXEC] Move to q=%d r=%d (pathfinding not yet implemented)\n",
                   act->q, act->r);
            fflush(stdout);
            exec->actions_executed++;
            return true;

        case AI_ACTION_LOOK_AT:
            // TODO: agent facing (Phase 2)
            printf("[AI-EXEC] Look at q=%d r=%d (not yet implemented)\n",
                   act->q, act->r);
            fflush(stdout);
            exec->actions_executed++;
            return true;
    }

    return false;
}

void ai_actions_clear(AiActionExecutor* exec) {
    exec->queue_count = 0;
    exec->queue_head = 0;
    exec->action_timer = 0.0f;
}

bool ai_actions_busy(const AiActionExecutor* exec) {
    return exec->queue_head < exec->queue_count;
}

#ifndef AI_AGENT_H
#define AI_AGENT_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "sokol_gfx.h"
#include "hex_terrain.h"
#include "ai_npc.h"
#include "ai_actions.h"
#include "ai_pathfind.h"

// --- AI Agent ---
// A physical NPC entity that exists in the game world.
// Mirrors the player controller: gravity, collision, ground detection.

#define AI_AGENT_MAX 8
#define AI_AGENT_MESH_MAX_VERTS 2048

typedef enum {
    AGENT_STATE_IDLE,
    AGENT_STATE_WALKING,
    AGENT_STATE_WORKING,
    AGENT_STATE_TALKING,
    AGENT_STATE_SLEEPING,
} AgentState;

typedef struct {
    // Identity
    char id[64];
    char name[32];
    char directive[512];
    float creativity;
    float chattiness;
    float diligence;

    // Appearance (from JSON)
    float height;
    float stockiness;
    float head_scale;
    float arm_length;
    float torso_width;
    float leg_length;
    float primary_color[3];
    float accent_color[3];

    // Physics (mirrors Camera)
    double pos_d[3];           // double-precision world position
    HMM_Vec3 position;        // float copy for rendering
    HMM_Vec3 local_up;        // radial up (normalized position on sphere)
    HMM_Vec3 forward;         // facing direction in tangent plane
    HMM_Vec3 velocity;
    float yaw;                // facing angle in tangent frame
    bool on_ground;
    int gravity_body;          // -1 = planet, 0-9 = moon

    // State machine
    AgentState state;
    float state_timer;         // time in current state

    // Navigation
    AiPath path;
    int target_q, target_r;    // hex grid target for pathfinding
    bool has_target;
    float move_speed;          // meters per second

    // Rendering — CPU-side model-space verts (9 floats each: pos, normal, color)
    float* mesh_verts;        // malloc'd array of all verts (primary + accent + eyes)
    int vert_count_total;
    bool mesh_valid;

    // AI connection
    AiNpc ai;
    AiActionExecutor executor;

    // Sleep range
    float sleep_range;         // distance from player to sleep (default 200m)
    bool active;               // false = not spawned or sleeping
} AiAgent;

typedef struct {
    AiAgent agents[AI_AGENT_MAX];
    int count;
    float planet_radius;
} AiAgentSystem;

// Initialize the agent system
void ai_agent_system_init(AiAgentSystem* sys, float planet_radius);

// Load an agent from a JSON file (cache/agents/*.json).
// Returns the agent index, or -1 on failure.
int ai_agent_load(AiAgentSystem* sys, const char* json_path);

// Spawn an agent near a world position.
// Places it on the ground at the nearest walkable hex.
void ai_agent_spawn(AiAgentSystem* sys, int idx, HexTerrain* ht,
                    double spawn_x, double spawn_y, double spawn_z);

// Update all agents: gravity, movement, state machine, AI polling.
// Call once per frame.
void ai_agent_system_update(AiAgentSystem* sys, HexTerrain* ht,
                            const double player_pos[3], float dt);

// Render all active agents. Call during the render pass.
// vp: view-projection matrix, world_origin: camera-relative origin
void ai_agent_system_render(AiAgentSystem* sys, HMM_Mat4 vp,
                            const double world_origin[3],
                            HMM_Vec3 sun_dir);

// Shutdown: free GPU resources
void ai_agent_system_shutdown(AiAgentSystem* sys);

// Find the nearest agent to a world position within max_dist.
// Returns agent index or -1.
int ai_agent_nearest(const AiAgentSystem* sys, const double pos[3], float max_dist);

// Save all agent positions to a file (e.g. cache/worlds/{id}/agents.dat)
void ai_agent_save_positions(const AiAgentSystem* sys, const char* path);

// Load agent positions from file. Matches agents by ID.
void ai_agent_load_positions(AiAgentSystem* sys, const char* path);

#endif

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
#include "ai_sensors.h"
#include "ai_script.h"
#include "ai_memory.h"
#include "ai_emotions.h"
#include "ai_vision.h"

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
    AGENT_STATE_JUMPING,
    AGENT_STATE_PICKUP,
    AGENT_STATE_PLACING,
    AGENT_STATE_WAVING,
    AGENT_STATE_CELEBRATING,
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

    // Current hex position (updated each frame)
    int hex_q, hex_r;          // global hex grid coords
    int hex_layer;             // current ground layer
    int hex_ground_layer;      // topmost solid layer at this column

    // Navigation
    AiPath path;
    int target_q, target_r;    // hex grid target for pathfinding
    bool has_target;
    float move_speed;          // meters per second
    float stuck_timer;         // time spent not making progress
    float last_dist_to_step;   // previous frame distance to next path step
    float last_ground_r;       // previous frame ground height (for bob detection)
    int   bob_count;           // consecutive frames where ground_r flipped
    float bob_locked_r;        // locked ground height when bobbing detected

    // Conversation
    float convo_timer;         // time since last player interaction (0 = not in convo)
    bool in_conversation;      // player is actively chatting

    // Rendering — CPU-side model-space verts (9 floats each: pos, normal, color)
    float* mesh_verts;
    int vert_count_total;
    bool mesh_valid;

    // Body part vertex ranges (for animation)
    int vr_torso_start, vr_torso_count;       // torso + shoulders + pelvis
    int vr_left_leg_start, vr_left_leg_count;
    int vr_right_leg_start, vr_right_leg_count;
    int vr_head_start, vr_head_count;          // neck + head + eyes
    int vr_left_arm_start, vr_left_arm_count;
    int vr_right_arm_start, vr_right_arm_count;

    // Pivot points (model-space Y)
    float pivot_hip_y;     // leg rotation pivot
    float pivot_shoulder_y; // arm rotation pivot
    float pivot_neck_y;    // head rotation pivot
    float pivot_arm_x;     // shoulder X offset from center

    // Animation
    float anim_time;

    // AI connection
    AiNpc ai;
    AiActionExecutor executor;

    // Scripting, sensors & memory
    AiScriptRunner script_runner;
    AiScanResult last_scan;
    char sensor_report[AI_SENSOR_REPORT_MAX];
    AiMemory memory;
    AiEmotions emotions;
    AiVision vision;

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

// Raycast against agent bounding boxes. Returns agent index or -1.
// ray_origin/ray_dir in world space (float). max_dist in meters.
int ai_agent_raycast(const AiAgentSystem* sys, HMM_Vec3 ray_origin,
                     HMM_Vec3 ray_dir, float max_dist);

// Save all agent positions to a file (e.g. cache/worlds/{id}/agents.dat)
void ai_agent_save_positions(const AiAgentSystem* sys, const char* path);

// Load agent positions from file. Matches agents by ID.
void ai_agent_load_positions(AiAgentSystem* sys, const char* path);

#endif

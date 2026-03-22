#include "ai_agent.h"
#include "planet.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ---- Minimal JSON parser (reuse pattern from lobby.c) ----

static bool json_str(const char* json, const char* key, char* out, int out_size) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    if (*p != '"') return false;
    p++;
    int i = 0;
    while (*p && *p != '"' && i < out_size - 1) {
        if (*p == '\\' && *(p+1)) { p++; }
        out[i++] = *p++;
    }
    out[i] = '\0';
    return true;
}

static bool json_float(const char* json, const char* key, float* out) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    *out = (float)atof(p);
    return true;
}

// Parse a JSON array of 3 floats: [r, g, b]
static bool json_float3(const char* json, const char* key, float out[3]) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p = strchr(p, '[');
    if (!p) return false;
    p++;
    out[0] = (float)atof(p);
    p = strchr(p, ','); if (!p) return false; p++;
    out[1] = (float)atof(p);
    p = strchr(p, ','); if (!p) return false; p++;
    out[2] = (float)atof(p);
    return true;
}

// ---- Hex prism mesh generation (C port of HexBodyGenerator) ----

typedef struct { float x, y, z, nx, ny, nz, r, g, b; } AgentVertex; // 36 bytes = LodVertex layout

static int g_vert_count;
static AgentVertex g_verts[AI_AGENT_MESH_MAX_VERTS];

static void emit_tri(float x0,float y0,float z0,
                     float x1,float y1,float z1,
                     float x2,float y2,float z2,
                     float cr, float cg, float cb) {
    if (g_vert_count + 3 > AI_AGENT_MESH_MAX_VERTS) return;
    // Compute flat normal
    float ux = x1-x0, uy = y1-y0, uz = z1-z0;
    float vx = x2-x0, vy = y2-y0, vz = z2-z0;
    float nx = uy*vz - uz*vy;
    float ny = uz*vx - ux*vz;
    float nz = ux*vy - uy*vx;
    float len = sqrtf(nx*nx + ny*ny + nz*nz);
    if (len > 0.0001f) { nx /= len; ny /= len; nz /= len; }

    AgentVertex* v = &g_verts[g_vert_count];
    v[0] = (AgentVertex){x0,y0,z0, nx,ny,nz, cr,cg,cb};
    v[1] = (AgentVertex){x1,y1,z1, nx,ny,nz, cr,cg,cb};
    v[2] = (AgentVertex){x2,y2,z2, nx,ny,nz, cr,cg,cb};
    g_vert_count += 3;
}

static void hex_prism(float cx, float cy, float cz,
                      float radius, float height,
                      float wscale, float dscale,
                      float cr, float cg, float cb) {
    float halfH = height * 0.5f;
    float topX[6], topY[6], topZ[6];
    float botX[6], botY[6], botZ[6];

    for (int i = 0; i < 6; i++) {
        float angle = 3.14159265f / 3.0f * i + 3.14159265f / 6.0f;
        float x = cx + radius * cosf(angle) * wscale;
        float z = cz + radius * sinf(angle) * dscale;
        topX[i] = x; topY[i] = cy + halfH; topZ[i] = z;
        botX[i] = x; botY[i] = cy - halfH; botZ[i] = z;
    }

    // Top face (CCW from above)
    for (int i = 0; i < 6; i++) {
        int n = (i + 1) % 6;
        emit_tri(cx, cy+halfH, cz, topX[n], topY[n], topZ[n], topX[i], topY[i], topZ[i], cr,cg,cb);
    }
    // Bottom face
    for (int i = 0; i < 6; i++) {
        int n = (i + 1) % 6;
        emit_tri(cx, cy-halfH, cz, botX[i], botY[i], botZ[i], botX[n], botY[n], botZ[n], cr,cg,cb);
    }
    // Sides
    for (int i = 0; i < 6; i++) {
        int n = (i + 1) % 6;
        emit_tri(topX[i],topY[i],topZ[i], topX[n],topY[n],topZ[n], botX[n],botY[n],botZ[n], cr,cg,cb);
        emit_tri(topX[i],topY[i],topZ[i], botX[n],botY[n],botZ[n], botX[i],botY[i],botZ[i], cr,cg,cb);
    }
}

// Z-facing hex prism (for eyes)
static void hex_prism_z(float cx, float cy, float cz,
                        float radius, float depth,
                        float cr, float cg, float cb) {
    float halfD = depth * 0.5f;
    float fX[6], fY[6], bX[6], bY[6];

    for (int i = 0; i < 6; i++) {
        float angle = 3.14159265f / 3.0f * i + 3.14159265f / 6.0f;
        float x = cx + radius * cosf(angle);
        float y = cy + radius * sinf(angle);
        fX[i] = x; fY[i] = y;
        bX[i] = x; bY[i] = y;
    }
    float fz = cz + halfD, bz = cz - halfD;

    // Front face
    for (int i = 0; i < 6; i++) {
        int n = (i + 1) % 6;
        emit_tri(cx,cy,fz, fX[i],fY[i],fz, fX[n],fY[n],fz, cr,cg,cb);
    }
    // Back face
    for (int i = 0; i < 6; i++) {
        int n = (i + 1) % 6;
        emit_tri(cx,cy,bz, bX[n],bY[n],bz, bX[i],bY[i],bz, cr,cg,cb);
    }
    // Sides
    for (int i = 0; i < 6; i++) {
        int n = (i + 1) % 6;
        emit_tri(fX[n],fY[n],fz, fX[i],fY[i],fz, bX[i],bY[i],bz, cr,cg,cb);
        emit_tri(fX[n],fY[n],fz, bX[i],bY[i],bz, bX[n],bY[n],bz, cr,cg,cb);
    }
}

static void generate_agent_mesh(AiAgent* agent) {
    float s = agent->height;
    float stocky = 0.5f + agent->stockiness * 0.5f;
    float hexR = 0.15f;

    float torsoW = stocky * agent->torso_width;
    float torsoR = hexR * (0.8f + torsoW * 0.8f);
    float torsoH = 0.4f * s;
    float torsoY = 0.25f * s;

    float pelvisR = torsoR * 1.05f;
    float pelvisH = 0.06f * s;
    float pelvisY = torsoY - torsoH/2 - pelvisH/2;

    float neckR = hexR * 0.5f;
    float neckH = 0.06f * s;
    float neckY = torsoY + torsoH/2 + neckH/2;

    float headR = hexR * 1.1f * agent->head_scale;
    float headH = 0.22f * s * agent->head_scale;
    float headY = neckY + neckH/2 + headH/2;

    float padR = hexR * 0.65f;
    float padH = 0.1f * s;
    float padSpacing = torsoR * 0.85f;
    float padY = torsoY + torsoH/2 - padH * 0.3f;

    float armR = hexR * 0.55f;
    float armH = 0.38f * s * agent->arm_length;
    float armSpacing = padSpacing + padR * 0.3f;
    float armY = padY - padH/2 - armH/2;

    float legR = hexR * 0.55f;
    float legH = 0.4f * s * agent->leg_length;
    float legSpacing = hexR * stocky * 0.75f;
    float legY = pelvisY - pelvisH/2 - legH/2;

    // Anchor feet at Y=0
    float feetY = legY - legH/2;
    float yOff = -feetY;
    torsoY += yOff; pelvisY += yOff; neckY += yOff; headY += yOff;
    padY += yOff; armY += yOff; legY += yOff;

    float* pc = agent->primary_color;
    float* ac = agent->accent_color;

    // Generate all parts, tracking vertex ranges per body part
    g_vert_count = 0;
    int mark;

    // Torso + shoulders + pelvis
    mark = g_vert_count;
    hex_prism(0, torsoY, 0, torsoR, torsoH, stocky, 0.8f, pc[0], pc[1], pc[2]);
    hex_prism(-padSpacing, padY, 0, padR, padH, 0.8f, 0.8f, pc[0], pc[1], pc[2]);
    hex_prism(padSpacing, padY, 0, padR, padH, 0.8f, 0.8f, pc[0], pc[1], pc[2]);
    hex_prism(0, pelvisY, 0, pelvisR, pelvisH, stocky, 0.7f, pc[0], pc[1], pc[2]);
    agent->vr_torso_start = mark; agent->vr_torso_count = g_vert_count - mark;

    // Left leg
    mark = g_vert_count;
    hex_prism(-legSpacing, legY, 0, legR, legH, stocky*0.65f, 0.65f, pc[0], pc[1], pc[2]);
    agent->vr_left_leg_start = mark; agent->vr_left_leg_count = g_vert_count - mark;

    // Right leg
    mark = g_vert_count;
    hex_prism(legSpacing, legY, 0, legR, legH, stocky*0.65f, 0.65f, pc[0], pc[1], pc[2]);
    agent->vr_right_leg_start = mark; agent->vr_right_leg_count = g_vert_count - mark;

    // Left arm
    mark = g_vert_count;
    hex_prism(-armSpacing, armY, 0, armR, armH, 0.55f, 0.55f, ac[0], ac[1], ac[2]);
    agent->vr_left_arm_start = mark; agent->vr_left_arm_count = g_vert_count - mark;

    // Right arm
    mark = g_vert_count;
    hex_prism(armSpacing, armY, 0, armR, armH, 0.55f, 0.55f, ac[0], ac[1], ac[2]);
    agent->vr_right_arm_start = mark; agent->vr_right_arm_count = g_vert_count - mark;

    // Head (neck + head + eyes)
    mark = g_vert_count;
    hex_prism(0, neckY, 0, neckR, neckH, 1.0f, 1.0f, ac[0], ac[1], ac[2]);
    hex_prism(0, headY, 0, headR, headH, 1.0f, 1.0f, ac[0], ac[1], ac[2]);
    float eyeR = headR * 0.22f;
    float eyeDepth = headR * 0.3f;
    float eyeSpacing = headR * 0.45f;
    float eyeY = headY + headH * 0.08f;
    float eyeZ = headR * 0.7f;
    hex_prism_z(-eyeSpacing, eyeY, eyeZ, eyeR, eyeDepth, 1,1,1);
    hex_prism_z(eyeSpacing, eyeY, eyeZ, eyeR, eyeDepth, 1,1,1);
    agent->vr_head_start = mark; agent->vr_head_count = g_vert_count - mark;

    // Store pivot points
    agent->pivot_hip_y = pelvisY - pelvisH * 0.5f;
    agent->pivot_shoulder_y = padY;
    agent->pivot_neck_y = neckY - neckH * 0.5f;

    // Store on CPU
    agent->vert_count_total = g_vert_count;
    agent->mesh_verts = (float*)malloc(g_vert_count * 9 * sizeof(float));
    memcpy(agent->mesh_verts, g_verts, g_vert_count * 9 * sizeof(float));

    agent->mesh_valid = true;
    printf("[AGENT] Generated mesh for '%s': %d verts\n", agent->name, g_vert_count);
    fflush(stdout);
}

// ---- Public API ----

void ai_agent_system_init(AiAgentSystem* sys, float planet_radius) {
    memset(sys, 0, sizeof(AiAgentSystem));
    sys->planet_radius = planet_radius;
}

int ai_agent_load(AiAgentSystem* sys, const char* json_path) {
    if (sys->count >= AI_AGENT_MAX) return -1;

    FILE* f = fopen(json_path, "r");
    if (!f) {
        printf("[AGENT] Failed to open: %s\n", json_path);
        return -1;
    }
    char json[4096];
    int len = (int)fread(json, 1, sizeof(json) - 1, f);
    json[len] = '\0';
    fclose(f);

    int idx = sys->count++;
    AiAgent* a = &sys->agents[idx];
    memset(a, 0, sizeof(AiAgent));

    // Parse identity
    json_str(json, "id", a->id, sizeof(a->id));
    json_str(json, "name", a->name, sizeof(a->name));
    json_str(json, "directive", a->directive, sizeof(a->directive));

    // Parse appearance
    a->height = 1.0f; a->stockiness = 0.5f; a->head_scale = 1.0f;
    a->arm_length = 1.0f; a->torso_width = 1.0f; a->leg_length = 1.0f;
    json_float(json, "height", &a->height);
    json_float(json, "stockiness", &a->stockiness);
    json_float(json, "head_scale", &a->head_scale);
    json_float(json, "arm_length", &a->arm_length);
    json_float(json, "torso_width", &a->torso_width);
    json_float(json, "leg_length", &a->leg_length);
    json_float3(json, "primary_color", a->primary_color);
    json_float3(json, "accent_color", a->accent_color);

    // Parse personality
    a->creativity = 0.5f; a->chattiness = 0.5f; a->diligence = 0.7f;
    json_float(json, "creativity", &a->creativity);
    json_float(json, "chattiness", &a->chattiness);
    json_float(json, "diligence", &a->diligence);

    // Defaults
    a->move_speed = 2.0f;
    a->sleep_range = 200.0f;
    a->gravity_body = -1;
    a->state = AGENT_STATE_IDLE;

    // Initialize AI subsystems — agent shares the already-running llama-server
    ai_actions_init(&a->executor);
    ai_script_runner_init(&a->script_runner);
    ai_emotions_init(&a->emotions);

    // Set up agent's AI — load provider config (local or Claude)
    memset(&a->ai, 0, sizeof(AiNpc));
    a->ai.server_port = 8080;
    ai_load_config(&a->ai);
    if (a->ai.provider == AI_PROVIDER_LOCAL) {
        a->ai.server_running = true; // assume server launched by main ai_npc_init
    } else {
        a->ai.server_running = true; // Claude is always "running"
    }

    // Initialize persistent memory
    {
        char safe_name[64];
        snprintf(safe_name, sizeof(safe_name), "%s", a->name);
        // Lowercase for filesystem
        for (int ci = 0; safe_name[ci]; ci++)
            if (safe_name[ci] >= 'A' && safe_name[ci] <= 'Z')
                safe_name[ci] += 32;

        char agent_dir[256];
        snprintf(agent_dir, sizeof(agent_dir), "cache/agents/%s", safe_name);
        ai_memory_init(&a->memory, agent_dir);
        ai_memory_load(&a->memory);
        ai_memory_load_chat_history(&a->memory, &a->ai);
    }

    // Build system prompt from directives + personality + memory
    {
        char directives[512] = "";
        FILE* df = fopen("config/ai_directives.txt", "r");
        if (df) {
            int dlen = (int)fread(directives, 1, sizeof(directives) - 1, df);
            directives[dlen] = '\0';
            fclose(df);
        }

        // Build memory context
        char memory_ctx[AI_MEMORY_MAX + AI_JOURNAL_MAX + 256];
        ai_memory_build_context(&a->memory, memory_ctx, sizeof(memory_ctx));

        if (a->ai.provider == AI_PROVIDER_CLAUDE) {
            snprintf(a->ai.system_prompt, sizeof(a->ai.system_prompt),
                "%s\n\nYour name is %s. %s\n%s",
                directives, a->name, a->directive, memory_ctx);
        } else {
            snprintf(a->ai.system_prompt, sizeof(a->ai.system_prompt),
                "%s\n\nYour name is %s. %s\n%s\n"
                "Respond with dialogue and a list of actions in JSON format.",
                directives, a->name, a->directive, memory_ctx);
        }
    }

    printf("[AGENT] Loaded '%s' (id=%s) from %s\n", a->name, a->id, json_path);
    fflush(stdout);

    return idx;
}

void ai_agent_spawn(AiAgentSystem* sys, int idx, HexTerrain* ht,
                    double spawn_x, double spawn_y, double spawn_z) {
    if (idx < 0 || idx >= sys->count) return;
    AiAgent* a = &sys->agents[idx];

    // Place agent at spawn position
    a->pos_d[0] = spawn_x;
    a->pos_d[1] = spawn_y;
    a->pos_d[2] = spawn_z;
    a->position = HMM_V3((float)spawn_x, (float)spawn_y, (float)spawn_z);

    // Compute local up (radial direction from planet center)
    float len = HMM_LenV3(a->position);
    a->local_up = (len > 0.001f) ? HMM_MulV3F(a->position, 1.0f / len) : HMM_V3(0,1,0);

    // Snap to ground
    int gcol, grow, layer;
    hex_terrain_world_to_hex(ht, a->position, &gcol, &grow, &layer);
    float ground_r = hex_terrain_ground_height(ht, gcol, grow);
    if (ground_r > 0.0f) {
        // Place feet on ground + eye_height to stand on surface
        float stand_r = ground_r + 1.7f;
        a->pos_d[0] = (double)a->local_up.X * stand_r;
        a->pos_d[1] = (double)a->local_up.Y * stand_r;
        a->pos_d[2] = (double)a->local_up.Z * stand_r;
        a->position = HMM_V3((float)a->pos_d[0], (float)a->pos_d[1], (float)a->pos_d[2]);
        printf("[AGENT] Ground snap: r=%.1f -> stand_r=%.1f\n", ground_r, stand_r);
    } else {
        printf("[AGENT] WARNING: no ground found at gcol=%d grow=%d\n", gcol, grow);
    }
    fflush(stdout);

    // Initial facing: tangent to sphere, roughly "north"
    HMM_Vec3 world_up = HMM_V3(0, 1, 0);
    if (fabsf(HMM_DotV3(a->local_up, world_up)) > 0.99f)
        world_up = HMM_V3(1, 0, 0);
    a->forward = HMM_NormV3(HMM_SubV3(world_up, HMM_MulV3F(a->local_up, HMM_DotV3(world_up, a->local_up))));
    a->yaw = 0;
    a->on_ground = true;

    // Generate mesh
    generate_agent_mesh(a);

    a->active = true;
    a->state = AGENT_STATE_IDLE;

    printf("[AGENT] Spawned '%s' at (%.1f, %.1f, %.1f)\n",
           a->name, a->position.X, a->position.Y, a->position.Z);
    fflush(stdout);
}

void ai_agent_system_update(AiAgentSystem* sys, HexTerrain* ht,
                            const double player_pos[3], float dt) {
    for (int i = 0; i < sys->count; i++) {
        AiAgent* a = &sys->agents[i];
        if (!a->active) continue;

        // Distance to player
        double dx = a->pos_d[0] - player_pos[0];
        double dy = a->pos_d[1] - player_pos[1];
        double dz = a->pos_d[2] - player_pos[2];
        float dist = (float)sqrt(dx*dx + dy*dy + dz*dz);

        // Sleep/wake based on distance
        if (dist > a->sleep_range && a->state != AGENT_STATE_SLEEPING) {
            a->state = AGENT_STATE_SLEEPING;
            printf("[AGENT] '%s' sleeping (player too far: %.0fm)\n", a->name, dist);
            fflush(stdout);
            continue;
        }
        if (dist <= a->sleep_range && a->state == AGENT_STATE_SLEEPING) {
            a->state = AGENT_STATE_IDLE;
            printf("[AGENT] '%s' woke up (player nearby: %.0fm)\n", a->name, dist);
            fflush(stdout);
        }
        if (a->state == AGENT_STATE_SLEEPING) continue;

        a->state_timer += dt;
        a->anim_time += dt;
        ai_emotions_update(&a->emotions, dt, a->state);

        // Gravity: pull toward planet center
        a->position = HMM_V3((float)a->pos_d[0], (float)a->pos_d[1], (float)a->pos_d[2]);
        // Use double for radius to avoid float precision bobbing at large distances
        double r_d = sqrt(a->pos_d[0]*a->pos_d[0] + a->pos_d[1]*a->pos_d[1] + a->pos_d[2]*a->pos_d[2]);
        float r = (float)r_d;
        if (r > 0.001f) {
            a->local_up = HMM_MulV3F(a->position, 1.0f / r);
        }

        // Ground detection via hex terrain
        int gcol, grow, layer;
        hex_terrain_world_to_hex(ht, a->position, &gcol, &grow, &layer);
        // Compute ground_r in double to match r_d precision
        double ground_r_d = (double)hex_terrain_ground_height(ht, gcol, grow);
        float ground_r = (float)ground_r_d;

        // Store current hex position for other systems
        a->hex_q = gcol;
        a->hex_r = grow;
        a->hex_layer = layer;
        if (ground_r > 0.0f) {
            float base = hex_terrain_col_base_r(ht, gcol, grow);
            a->hex_ground_layer = (base > 0.0f) ? (int)((ground_r - base) / HEX_HEIGHT) : 0;
        }

        if (ground_r_d > 0.0) {
            // Ground comparison in double precision to avoid float bobbing
            double height_above_d = r_d - ground_r_d;
            float height_above = (float)height_above_d;
            bool jumping_up = (a->state == AGENT_STATE_JUMPING &&
                              HMM_DotV3(a->velocity, a->local_up) > 0);
            // Allow falling from any height up to 20m (comfortable cliff jump)
            bool falling = (height_above > 0.5f && a->state != AGENT_STATE_SLEEPING);

            if (jumping_up || falling) {
                // In a jump arc — apply gravity, don't snap
                a->on_ground = false;
                float grav = 10.0f * dt;
                a->velocity = HMM_SubV3(a->velocity, HMM_MulV3F(a->local_up, grav));
                a->pos_d[0] += a->velocity.X * dt;
                a->pos_d[1] += a->velocity.Y * dt;
                a->pos_d[2] += a->velocity.Z * dt;
            } else {
                // Not jumping — snap to ground (in double precision)
                a->on_ground = true;
                double feet_r = ground_r_d + 0.05;
                double norm = 1.0 / r_d;
                a->pos_d[0] = a->pos_d[0] * norm * feet_r;
                a->pos_d[1] = a->pos_d[1] * norm * feet_r;
                a->pos_d[2] = a->pos_d[2] * norm * feet_r;
                a->velocity = HMM_V3(0, 0, 0);
            }
        }

        // ---- Pathfinding-based walking (+ jumping + stuck detection) ----
        // Timeout or invalid path: give up
        if ((a->state == AGENT_STATE_WALKING || a->state == AGENT_STATE_JUMPING)
            && (!a->path.valid || a->state_timer > 30.0f)) {
            a->state = AGENT_STATE_IDLE;
            a->has_target = false;
            a->path.valid = false;
            a->stuck_timer = 0.0f;
        }
        if ((a->state == AGENT_STATE_WALKING || a->state == AGENT_STATE_JUMPING)
            && a->path.valid) {
            int next_q, next_r;
            if (ai_path_next(&a->path, &next_q, &next_r)) {
                // Get world position & ground height of next step
                float wx, wy, wz, ux, uy, uz;
                float next_ground = hex_terrain_ground_height(ht, next_q, next_r);
                hex_terrain_hex_to_world(ht, next_q, next_r,
                    (int)((next_ground
                           - hex_terrain_col_base_r(ht, next_q, next_r)) / HEX_HEIGHT),
                    &wx, &wy, &wz, &ux, &uy, &uz);

                // Distance to next step (tangent plane)
                HMM_Vec3 target = HMM_V3(wx, wy, wz);
                HMM_Vec3 to_target = HMM_SubV3(target, a->position);
                float dot_up = HMM_DotV3(to_target, a->local_up);
                to_target = HMM_SubV3(to_target, HMM_MulV3F(a->local_up, dot_up));
                float dist_to_step = HMM_LenV3(to_target);

                float height_diff = next_ground - ground_r;

                // ---- Stuck detection ----
                // If not making progress toward next step, we're stuck
                if (a->on_ground && a->state == AGENT_STATE_WALKING) {
                    if (a->last_dist_to_step > 0.0f &&
                        dist_to_step >= a->last_dist_to_step - 0.01f) {
                        a->stuck_timer += dt;
                    } else {
                        a->stuck_timer = 0.0f;
                    }
                    a->last_dist_to_step = dist_to_step;

                    if (a->stuck_timer > 0.5f) {
                        // Stuck! Decide what to do based on height difference
                        if (height_diff > 0.3f && height_diff <= HEX_HEIGHT) {
                            // 1-block wall: jump over it
                            float jump_v = sqrtf(2.0f * 10.0f * (height_diff + 0.5f));
                            a->velocity = HMM_MulV3F(a->local_up, jump_v);
                            a->on_ground = false;
                            a->state = AGENT_STATE_JUMPING;
                            a->state_timer = 0.0f;
                            a->anim_time = 0.0f;
                            a->stuck_timer = 0.0f;
                        } else if (height_diff > HEX_HEIGHT) {
                            // Too high to jump — place a step block and climb
                            HexHitResult hit;
                            memset(&hit, 0, sizeof(hit));
                            hit.valid = true;
                            hit.place_gcol = gcol;
                            hit.place_grow = grow;
                            // Place at current ground + 1 layer (create a step)
                            float base_r = hex_terrain_col_base_r(ht, gcol, grow);
                            hit.place_layer = (int)((ground_r - base_r) / HEX_HEIGHT) + 1;
                            if (hex_terrain_place(ht, &hit, VOXEL_STONE)) {
                                printf("[AGENT] '%s' placed step block at q=%d r=%d layer=%d\n",
                                       a->name, gcol, grow, hit.place_layer);
                                fflush(stdout);
                            }
                            a->stuck_timer = 0.0f;
                        } else {
                            // Same height but stuck (collision with something)
                            // Skip this path step or abandon path
                            a->stuck_timer = 0.0f;
                            ai_path_advance(&a->path);
                            if (ai_path_done(&a->path)) {
                                a->state = AGENT_STATE_IDLE;
                                a->has_target = false;
                            }
                        }
                    }
                }

                // ---- Proactive jump (before getting stuck) ----
                if (height_diff > 0.3f && height_diff <= HEX_HEIGHT
                    && a->on_ground && a->state == AGENT_STATE_WALKING
                    && dist_to_step < 2.0f) {
                    float jump_v = sqrtf(2.0f * 10.0f * (height_diff + 0.5f));
                    a->velocity = HMM_MulV3F(a->local_up, jump_v);
                    a->on_ground = false;
                    a->state = AGENT_STATE_JUMPING;
                    a->state_timer = 0.0f;
                    a->anim_time = 0.0f;
                    a->stuck_timer = 0.0f;
                }

                // Return to walking after landing from jump
                if (a->state == AGENT_STATE_JUMPING && a->on_ground
                    && a->state_timer > 0.1f) {
                    a->state = AGENT_STATE_WALKING;
                    a->stuck_timer = 0.0f;
                    a->last_dist_to_step = 0.0f;
                }

                if (dist_to_step < 0.5f) {
                    // Reached this step, advance
                    ai_path_advance(&a->path);
                    a->stuck_timer = 0.0f;
                    a->last_dist_to_step = 0.0f;
                    if (ai_path_done(&a->path)) {
                        a->state = AGENT_STATE_IDLE;
                        a->has_target = false;
                    }
                } else {
                    // Walk/move toward step
                    HMM_Vec3 move_dir = HMM_MulV3F(to_target, 1.0f / dist_to_step);
                    // Speed tiers: normal 2m/s, sprint 4m/s at >50m, run 8m/s if following at >30m
                    float speed = a->move_speed;
                    bool following = (a->has_target && a->target_q == -9999);
                    if (following && dist > 30.0f) speed = 8.0f;
                    else if (dist > 50.0f) speed = 4.0f;
                    float step = speed * dt;
                    if (step > dist_to_step) step = dist_to_step;
                    a->pos_d[0] += move_dir.X * step;
                    a->pos_d[1] += move_dir.Y * step;
                    a->pos_d[2] += move_dir.Z * step;

                    // Tiny hop every ~8 steps to break out of ground-snap oscillation
                    if (a->on_ground && fmodf(a->anim_time, 2.5f) < dt) {
                        a->velocity = HMM_MulV3F(a->local_up, 1.5f);
                        a->on_ground = false;
                    }

                    // Smooth turn toward movement direction
                    a->forward = HMM_NormV3(HMM_LerpV3(a->forward, 8.0f * dt, move_dir));
                }
            } else {
                a->state = AGENT_STATE_IDLE;
            }
        }

        // ---- AI action execution ----
        ai_npc_poll(&a->ai);
        if (a->ai.state == AI_SUCCESS) {
            printf("[AGENT] '%s' got AI_SUCCESS (dialogue=%s, %d actions)\n",
                   a->name, a->ai.dialogue[0] ? "yes" : "EMPTY", a->ai.action_count);
            fflush(stdout);
            // Enqueue actions from LLM response
            ai_actions_enqueue(&a->executor, a->ai.actions, a->ai.action_count);
            // Play talking animation when agent responds with dialogue
            if (a->ai.dialogue[0]) {
                a->state = AGENT_STATE_TALKING;
                a->state_timer = 0.0f;
                a->anim_time = 0.0f;
            }
            // Don't reset ai.state here — let main.c read it for chat log first
        } else if (a->ai.state == AI_ERROR) {
            printf("[AGENT] '%s' got AI_ERROR: %s\n", a->name, a->ai.error);
            fflush(stdout);
            // Don't reset — let main.c read the error first
        }

        // Process action queue — walk to target before executing place/break
        if (ai_actions_busy(&a->executor) && a->state != AGENT_STATE_WALKING) {
            AiAction* next_act = &a->executor.queue[a->executor.queue_head];

            // For place/break actions, walk to the location first
            if ((next_act->type == AI_ACTION_PLACE_BLOCK ||
                 next_act->type == AI_ACTION_BREAK_BLOCK) &&
                a->state != AGENT_STATE_WORKING) {

                // Check if we're close enough to the target
                float wx, wy, wz, ux, uy, uz;
                hex_terrain_hex_to_world(ht, next_act->q, next_act->r, next_act->layer,
                    &wx, &wy, &wz, &ux, &uy, &uz);
                HMM_Vec3 target = HMM_V3(wx, wy, wz);
                float dist_to_block = HMM_LenV3(HMM_SubV3(target, a->position));

                if (dist_to_block > 3.0f) {
                    // Too far — pathfind there first
                    if (ai_pathfind(&a->path, ht, gcol, grow, next_act->q, next_act->r)) {
                        a->state = AGENT_STATE_WALKING;
                    } else {
                        // Can't reach — skip this action
                        printf("[AGENT] '%s' can't reach q=%d r=%d, skipping\n",
                               a->name, next_act->q, next_act->r);
                        a->executor.queue_head++;
                        a->executor.actions_failed++;
                    }
                } else {
                    // Close enough — execute the action
                    a->state = AGENT_STATE_WORKING;
                    ai_actions_update(&a->executor, ht, dt);
                }
            } else if (next_act->type == AI_ACTION_MOVE_TO) {
                // Move_to: pathfind to hex coords
                if (ai_pathfind(&a->path, ht, gcol, grow, next_act->q, next_act->r)) {
                    a->state = AGENT_STATE_WALKING;
                }
                a->executor.queue_head++;
                a->executor.actions_executed++;
            } else {
                // Say, look_at, etc — execute immediately
                ai_actions_update(&a->executor, ht, dt);
                a->state = AGENT_STATE_WORKING;
            }
        } else if (a->state == AGENT_STATE_WORKING && !ai_actions_busy(&a->executor)) {
            a->state = AGENT_STATE_IDLE;
        } else if (a->state == AGENT_STATE_TALKING && a->state_timer > 4.0f) {
            // Return to idle after talking animation plays
            a->state = AGENT_STATE_IDLE;
        } else if (a->state == AGENT_STATE_WAVING && a->state_timer > 3.0f) {
            a->state = AGENT_STATE_IDLE;
        }

        // ---- Conversation timeout ----
        // Stay put while in conversation; expire after 15s idle or 60s typing
        if (a->in_conversation) {
            a->convo_timer += dt;
            // If AI is pending (waiting for response), don't time out
            if (a->ai.state == AI_PENDING) {
                a->convo_timer = 0.0f;
            }
            // 60s if chat window is open (player typing), 15s otherwise
            float timeout = 15.0f; // will be overridden by main.c if chat is open
            if (a->convo_timer > timeout) {
                a->in_conversation = false;
                a->convo_timer = 0.0f;
            }

            // Smoothly face the player while in conversation
            if (a->state == AGENT_STATE_IDLE || a->state == AGENT_STATE_TALKING ||
                a->state == AGENT_STATE_WAVING) {
                HMM_Vec3 ppos = HMM_V3((float)player_pos[0], (float)player_pos[1], (float)player_pos[2]);
                HMM_Vec3 to_p = HMM_SubV3(ppos, a->position);
                float dp = HMM_DotV3(to_p, a->local_up);
                to_p = HMM_SubV3(to_p, HMM_MulV3F(a->local_up, dp));
                float pl = HMM_LenV3(to_p);
                if (pl > 0.001f) {
                    HMM_Vec3 target_fwd = HMM_MulV3F(to_p, 1.0f / pl);
                    a->forward = HMM_NormV3(HMM_LerpV3(a->forward, 5.0f * dt, target_fwd));
                }
            }
        }

        // ---- Script execution ----
        if (!ai_script_idle(&a->script_runner)) {
            // Auto-set anchor to agent's current hex pos if not set yet
            if (!a->script_runner.anchor_set) {
                a->script_runner.anchor_q = gcol;
                a->script_runner.anchor_r = grow;
                a->script_runner.anchor_layer = a->hex_ground_layer;
                a->script_runner.anchor_set = true;
            }

            ScriptActionType sact = ai_script_update(&a->script_runner, dt);
            if (sact != SCRIPT_ACT_NONE) {
                const AiScriptStep* step = ai_script_current_step(&a->script_runner);
                if (step) {
                    // Resolve relative coords to absolute
                    int abs_q = step->q, abs_r = step->r, abs_layer = step->layer;
                    if (sact == SCRIPT_ACT_MOVE_REL || sact == SCRIPT_ACT_PLACE_REL ||
                        sact == SCRIPT_ACT_BREAK_REL) {
                        abs_q = a->script_runner.anchor_q + step->dq;
                        abs_r = a->script_runner.anchor_r + step->dr;
                        abs_layer = a->script_runner.anchor_layer + step->dlayer;
                    }

                    switch (sact) {
                        case SCRIPT_ACT_MOVE_REL:
                        case SCRIPT_ACT_MOVE_TO:
                            if (ai_pathfind(&a->path, ht, gcol, grow, abs_q, abs_r)) {
                                a->state = AGENT_STATE_WALKING;
                                a->state_timer = 0.0f;
                                a->script_runner.state = SCRIPT_MOVING;
                            } else {
                                ai_script_advance(&a->script_runner);
                            }
                            break;
                        case SCRIPT_ACT_PLACE_REL:
                        case SCRIPT_ACT_PLACE:
                        case SCRIPT_ACT_BREAK_REL:
                        case SCRIPT_ACT_BREAK: {
                            // Phase 1: face the target block, enter BUILDING state
                            float bx, by, bz, bux, buy, buz;
                            hex_terrain_hex_to_world(ht, abs_q, abs_r, abs_layer,
                                &bx, &by, &bz, &bux, &buy, &buz);
                            HMM_Vec3 block_pos = HMM_V3(bx, by, bz);
                            HMM_Vec3 to_block = HMM_SubV3(block_pos, a->position);
                            float db = HMM_DotV3(to_block, a->local_up);
                            to_block = HMM_SubV3(to_block, HMM_MulV3F(a->local_up, db));
                            float bl = HMM_LenV3(to_block);
                            if (bl > 0.001f) {
                                HMM_Vec3 target_fwd = HMM_MulV3F(to_block, 1.0f / bl);
                                a->forward = HMM_NormV3(HMM_LerpV3(a->forward, 0.5f, target_fwd));
                            }
                            a->script_runner.state = SCRIPT_BUILDING;
                            a->script_runner.wait_timer = 0.5f; // cooldown
                            if (sact == SCRIPT_ACT_PLACE || sact == SCRIPT_ACT_PLACE_REL) {
                                a->state = AGENT_STATE_PLACING;
                            } else {
                                a->state = AGENT_STATE_WORKING;
                            }
                            a->state_timer = 0.0f;
                        } break;
                        case SCRIPT_ACT_SCAN:
                            ai_sensors_scan(ht, gcol, grow,
                                step->radius > 0 ? step->radius : 5, &a->last_scan);
                            ai_script_advance(&a->script_runner);
                            break;
                        case SCRIPT_ACT_SCAN_REPORT:
                            ai_sensors_scan(ht, gcol, grow, 5, &a->last_scan);
                            ai_sensors_report(&a->last_scan, a->sensor_report,
                                             sizeof(a->sensor_report));
                            ai_script_advance(&a->script_runner);
                            break;
                        case SCRIPT_ACT_JUMP: {
                            float jump_v = sqrtf(2.0f * 10.0f * 1.5f);
                            a->velocity = HMM_MulV3F(a->local_up, jump_v);
                            a->on_ground = false;
                            a->state = AGENT_STATE_JUMPING;
                            a->state_timer = 0.0f;
                            ai_script_advance(&a->script_runner);
                        } break;
                        case SCRIPT_ACT_SAY:
                            // Say is handled by main.c reading the step text
                            ai_script_advance(&a->script_runner);
                            break;
                        default:
                            ai_script_advance(&a->script_runner);
                            break;
                    }
                }
            }
            // If script was MOVING and we've arrived (path done), advance
            if (a->script_runner.state == SCRIPT_MOVING &&
                a->state == AGENT_STATE_IDLE) {
                a->script_runner.state = SCRIPT_RUNNING;
                ai_script_advance(&a->script_runner);
            }

            // BUILDING state: cooldown timer, then execute place/break
            if (a->script_runner.state == SCRIPT_BUILDING) {
                a->script_runner.wait_timer -= dt;

                // Keep facing the target during cooldown
                const AiScriptStep* bstep = ai_script_current_step(&a->script_runner);
                if (bstep) {
                    int bq = bstep->q, br = bstep->r, bl = bstep->layer;
                    if (bstep->action == SCRIPT_ACT_PLACE_REL || bstep->action == SCRIPT_ACT_BREAK_REL) {
                        bq = a->script_runner.anchor_q + bstep->dq;
                        br = a->script_runner.anchor_r + bstep->dr;
                        bl = a->script_runner.anchor_layer + bstep->dlayer;
                    }
                    float bwx, bwy, bwz, bux, buy, buz;
                    hex_terrain_hex_to_world(ht, bq, br, bl, &bwx, &bwy, &bwz, &bux, &buy, &buz);
                    HMM_Vec3 to_b = HMM_SubV3(HMM_V3(bwx, bwy, bwz), a->position);
                    float dp = HMM_DotV3(to_b, a->local_up);
                    to_b = HMM_SubV3(to_b, HMM_MulV3F(a->local_up, dp));
                    float tl = HMM_LenV3(to_b);
                    if (tl > 0.001f)
                        a->forward = HMM_NormV3(HMM_LerpV3(a->forward, 8.0f * dt, HMM_MulV3F(to_b, 1.0f / tl)));
                }

                if (a->script_runner.wait_timer <= 0.0f && bstep) {
                    int bq = bstep->q, br = bstep->r, bl = bstep->layer;
                    if (bstep->action == SCRIPT_ACT_PLACE_REL || bstep->action == SCRIPT_ACT_BREAK_REL) {
                        bq = a->script_runner.anchor_q + bstep->dq;
                        br = a->script_runner.anchor_r + bstep->dr;
                        bl = a->script_runner.anchor_layer + bstep->dlayer;
                    }

                    if (bstep->action == SCRIPT_ACT_PLACE || bstep->action == SCRIPT_ACT_PLACE_REL) {
                        // Remap layer relative to actual ground at anchor
                        // LLM gives absolute layers based on what it was told;
                        // correct by the delta between told vs actual ground
                        int layer_offset = bl - a->script_runner.anchor_layer;
                        int target_layer = a->hex_ground_layer + layer_offset;

                        HexHitResult hit = {0};
                        hit.valid = true;
                        hit.place_gcol = bq;
                        hit.place_grow = br;
                        hit.place_layer = target_layer;
                        uint8_t vtype = VOXEL_STONE;
                        if (strcmp(bstep->block, "dirt") == 0)  vtype = VOXEL_DIRT;
                        else if (strcmp(bstep->block, "grass") == 0) vtype = VOXEL_GRASS;
                        else if (strcmp(bstep->block, "sand") == 0)  vtype = VOXEL_SAND;
                        else if (strcmp(bstep->block, "ice") == 0)   vtype = VOXEL_ICE;
                        else if (strcmp(bstep->block, "torch") == 0) vtype = VOXEL_TORCH;
                        bool ok = hex_terrain_place(ht, &hit, vtype);
                        printf("[BUILD] %s place %s q=%d r=%d layer=%d (offset=%d): %s\n",
                               a->name, bstep->block[0] ? bstep->block : "stone",
                               bq, br, target_layer, layer_offset,
                               ok ? "OK" : "FAILED");
                        fflush(stdout);
                    } else {
                        int layer_offset = bl - a->script_runner.anchor_layer;
                        int target_layer = a->hex_ground_layer + layer_offset;

                        HexHitResult hit = {0};
                        hit.valid = true;
                        hit.gcol = bq;
                        hit.grow = br;
                        hit.layer = target_layer;
                        bool ok = hex_terrain_break(ht, &hit);
                        printf("[BUILD] %s break q=%d r=%d layer=%d (offset=%d): %s\n",
                               a->name, bq, br, target_layer, layer_offset,
                               ok ? "OK" : "FAILED");
                        fflush(stdout);
                    }

                    a->script_runner.state = SCRIPT_RUNNING;
                    ai_script_advance(&a->script_runner);
                }
            }
        }

        // ---- Follow player mode (target_q == -9999) ----
        if (a->has_target && a->target_q == -9999) {
            // Keep conversation alive while following (don't timeout)
            a->in_conversation = true;
            a->convo_timer = 0.0f;

            HMM_Vec3 ppos = HMM_V3((float)player_pos[0], (float)player_pos[1], (float)player_pos[2]);
            float pdist = HMM_LenV3(HMM_SubV3(ppos, a->position));

            // Re-pathfind when idle or when current path is outdated
            if (pdist > 5.0f && (a->state == AGENT_STATE_IDLE ||
                (a->state == AGENT_STATE_WALKING && a->state_timer > 2.0f))) {
                int pq, pr, pl;
                hex_terrain_world_to_hex(ht, ppos, &pq, &pr, &pl);
                if (ai_pathfind(&a->path, ht, gcol, grow, pq, pr)) {
                    a->state = AGENT_STATE_WALKING;
                    a->state_timer = 0.0f;
                    a->anim_time = 0.0f;
                }
            }
            // If close enough, idle and face player
            if (pdist <= 5.0f && a->state == AGENT_STATE_WALKING) {
                a->state = AGENT_STATE_IDLE;
                a->path.valid = false;
                HMM_Vec3 to_player = HMM_SubV3(ppos, a->position);
                float dot_p = HMM_DotV3(to_player, a->local_up);
                to_player = HMM_SubV3(to_player, HMM_MulV3F(a->local_up, dot_p));
                float plen = HMM_LenV3(to_player);
                if (plen > 0.001f) {
                    HMM_Vec3 target_fwd = HMM_MulV3F(to_player, 1.0f / plen);
                    a->forward = HMM_NormV3(HMM_LerpV3(a->forward, 5.0f * dt, target_fwd));
                }
            }
        }

        // ---- Boredom-driven building ----
        // When bored enough, self-prompt to build something creative
        if (ai_emotions_wants_to_build(&a->emotions)
            && a->state == AGENT_STATE_IDLE
            && a->ai.state == AI_IDLE
            && !a->in_conversation
            && ai_script_idle(&a->script_runner)) {

            // Pick a random creative prompt
            static const char* build_ideas[] = {
                "I'm bored. Build something cool nearby - a small sculpture or landmark.",
                "Nothing to do. Build a tiny watchtower or decorative pillar.",
                "I should practice my craft. Build a small stone arch.",
                "This spot needs something. Build a little stone bench or marker.",
                "Time to create! Build a small decorative hex column with torches.",
                "I'll build a cairn - a stack of stones as a trail marker.",
            };
            int pick = rand() % 6;

            // Scan surroundings for context
            ai_sensors_scan(ht, gcol, grow, 4, &a->last_scan);
            ai_sensors_report(&a->last_scan, a->sensor_report, sizeof(a->sensor_report));

            // Set up as a build request
            a->ai.script_mode = true;
            a->state = AGENT_STATE_WORKING;
            a->state_timer = 0.0f;

            char mood[256];
            ai_emotions_describe(&a->emotions, mood, sizeof(mood));
            snprintf(a->ai.context_prefix, sizeof(a->ai.context_prefix),
                "%s\n[SURROUNDINGS]\n%s\n"
                "[BUILD INSTRUCTIONS]\n"
                "You are at hex grid position q=%d, r=%d, ground_layer=%d.\n"
                "Output a JSON build script. Start with a SHORT 1-sentence comment, then the script.\n"
                "Format between [SCRIPT] and [/SCRIPT] tags.\n"
                "Keep it small - 5 to 15 steps max. You're building for fun, not a mega project.\n"
                "Available actions: place, break, move_to, say, wait, jump, scan\n"
                "Available blocks: stone, dirt, grass, sand, ice, torch\n"
                "Your position: q=%d r=%d. Ground layer is %d.\n"
                "Start blocks AT layer %d.\n",
                mood, a->sensor_report,
                gcol, grow, a->hex_ground_layer,
                gcol, grow, a->hex_ground_layer,
                a->hex_ground_layer);

            ai_npc_send(&a->ai, build_ideas[pick]);
            a->emotions.boredom = 0.0f;
            a->emotions.idle_time = 0.0f;

            printf("[AGENT] '%s' got bored, self-prompting: \"%s\"\n",
                   a->name, build_ideas[pick]);
            fflush(stdout);
        }

        // ---- Idle wandering ----
        // When idle with no pending actions or AI requests, wander randomly
        // Don't wander during conversation, or while running a script
        if (a->state == AGENT_STATE_IDLE && !ai_actions_busy(&a->executor)
            && a->ai.state == AI_IDLE && !a->in_conversation
            && ai_script_idle(&a->script_runner)
            && a->state_timer > 3.0f + (float)(i % 5)) {
            // Pick a random nearby walkable hex (offset by 3-8 cells)
            int wander_q = gcol, wander_r = grow;
            int steps = 3 + (rand() % 6);
            for (int s = 0; s < steps; s++) {
                int dir = rand() % 6;
                int nq, nr;
                // Offset-coordinate hex neighbors
                int parity = wander_r & 1;
                static const int even_dq[6] = {+1, 0,-1,-1, 0,+1};
                static const int even_dr[6] = { 0,-1,-1, 0,+1,+1};
                static const int odd_dq[6]  = {+1,+1, 0,-1,+1, 0};
                static const int odd_dr[6]  = { 0,-1,-1, 0,+1,+1};
                // Fix: use proper offset coords
                if (parity == 0) {
                    nq = wander_q + even_dq[dir];
                    nr = wander_r + even_dr[dir];
                } else {
                    nq = wander_q + odd_dq[dir];
                    nr = wander_r + odd_dr[dir];
                }
                // Check walkability (has ground)
                float gr = hex_terrain_ground_height(ht, nq, nr);
                if (gr > 0.0f) {
                    wander_q = nq;
                    wander_r = nr;
                }
            }
            // Only pathfind if we actually moved from current cell
            if (wander_q != gcol || wander_r != grow) {
                if (ai_pathfind(&a->path, ht, gcol, grow, wander_q, wander_r)) {
                    a->state = AGENT_STATE_WALKING;
                    a->state_timer = 0.0f;
                    a->anim_time = 0.0f;
                }
            }
        }
    }
}

void ai_agent_system_render(AiAgentSystem* sys, HMM_Mat4 vp,
                            const double world_origin[3],
                            HMM_Vec3 sun_dir) {
    // Agents use the same planet.glsl pipeline (LodVertex: pos, normal, color)
    // The pipeline and shader must be set by the caller before this call.
    for (int i = 0; i < sys->count; i++) {
        AiAgent* a = &sys->agents[i];
        if (!a->active || !a->mesh_valid || a->state == AGENT_STATE_SLEEPING) continue;

        // Camera-relative position
        float rx = (float)(a->pos_d[0] - world_origin[0]);
        float ry = (float)(a->pos_d[1] - world_origin[1]);
        float rz = (float)(a->pos_d[2] - world_origin[2]);

        // Build model matrix: translate to position, orient to local up + facing
        // The mesh is generated with Y-up at origin, feet at Y=0.
        // We need to rotate it so Y aligns with local_up and Z aligns with forward.
        HMM_Vec3 up = a->local_up;
        HMM_Vec3 fwd = a->forward;
        HMM_Vec3 right = HMM_NormV3(HMM_Cross(fwd, up));
        fwd = HMM_NormV3(HMM_Cross(up, right)); // re-orthogonalize

        HMM_Mat4 model = {
            .Elements = {
                {right.X, right.Y, right.Z, 0},
                {up.X,    up.Y,    up.Z,    0},
                {fwd.X,   fwd.Y,   fwd.Z,   0},
                {rx,      ry,      rz,      1},
            }
        };

        HMM_Mat4 mvp = HMM_MulM4(vp, model);
        (void)mvp;
        (void)sun_dir;

        // TODO: Draw with planet pipeline once we have access to the shader uniforms.
        // For now, agents are update-only (no rendering until wired into render.c).
    }
}

void ai_agent_system_shutdown(AiAgentSystem* sys) {
    for (int i = 0; i < sys->count; i++) {
        AiAgent* a = &sys->agents[i];
        if (a->mesh_valid && a->mesh_verts) {
            free(a->mesh_verts);
            a->mesh_verts = NULL;
        }
        ai_npc_shutdown(&a->ai);
    }
    sys->count = 0;
}

// ---- Raycast against agent bounding boxes ----

int ai_agent_raycast(const AiAgentSystem* sys, HMM_Vec3 ray_origin,
                     HMM_Vec3 ray_dir, float max_dist) {
    int best = -1;
    float best_t = max_dist;

    for (int i = 0; i < sys->count; i++) {
        const AiAgent* a = &sys->agents[i];
        if (!a->active || a->state == AGENT_STATE_SLEEPING) continue;

        // Agent bounding box in local space: [-0.3, 0.3] x [0, 1.5] x [-0.3, 0.3]
        // Transform ray into agent's local space
        HMM_Vec3 aup = a->local_up;
        HMM_Vec3 afwd = a->forward;
        HMM_Vec3 art = HMM_NormV3(HMM_Cross(afwd, aup));
        afwd = HMM_NormV3(HMM_Cross(aup, art));

        // Agent position (with 0.5m lift to match rendering)
        HMM_Vec3 apos = HMM_V3(
            (float)a->pos_d[0] + aup.X * 0.5f,
            (float)a->pos_d[1] + aup.Y * 0.5f,
            (float)a->pos_d[2] + aup.Z * 0.5f);

        // Ray origin relative to agent
        HMM_Vec3 rel = HMM_SubV3(ray_origin, apos);

        // Project into local space
        float ox = HMM_DotV3(rel, art);
        float oy = HMM_DotV3(rel, aup);
        float oz = HMM_DotV3(rel, afwd);
        float dx = HMM_DotV3(ray_dir, art);
        float dy = HMM_DotV3(ray_dir, aup);
        float dz = HMM_DotV3(ray_dir, afwd);

        // AABB: x [-0.3, 0.3], y [0, 1.5], z [-0.3, 0.3]
        float bmin[3] = {-0.3f, 0.0f, -0.3f};
        float bmax[3] = { 0.3f, 1.5f,  0.3f};
        float o[3] = {ox, oy, oz};
        float d[3] = {dx, dy, dz};

        float tmin = 0.0f, tmax = best_t;
        for (int ax = 0; ax < 3; ax++) {
            if (fabsf(d[ax]) < 1e-8f) {
                if (o[ax] < bmin[ax] || o[ax] > bmax[ax]) { tmin = tmax + 1; break; }
            } else {
                float t1 = (bmin[ax] - o[ax]) / d[ax];
                float t2 = (bmax[ax] - o[ax]) / d[ax];
                if (t1 > t2) { float tmp = t1; t1 = t2; t2 = tmp; }
                if (t1 > tmin) tmin = t1;
                if (t2 < tmax) tmax = t2;
            }
        }

        if (tmin <= tmax && tmin < best_t && tmin >= 0.0f) {
            best_t = tmin;
            best = i;
        }
    }
    return best;
}

// ---- Position persistence ----

typedef struct {
    char id[64];
    double pos[3];
    float yaw;
    int state;
} AgentSaveEntry;

void ai_agent_save_positions(const AiAgentSystem* sys, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) return;

    int count = sys->count;
    fwrite(&count, sizeof(int), 1, f);

    for (int i = 0; i < sys->count; i++) {
        const AiAgent* a = &sys->agents[i];
        AgentSaveEntry entry;
        memset(&entry, 0, sizeof(entry));
        snprintf(entry.id, sizeof(entry.id), "%s", a->id);
        entry.pos[0] = a->pos_d[0];
        entry.pos[1] = a->pos_d[1];
        entry.pos[2] = a->pos_d[2];
        entry.yaw = a->yaw;
        entry.state = (int)a->state;
        fwrite(&entry, sizeof(AgentSaveEntry), 1, f);
    }

    fclose(f);
    // printf("[AGENT] Saved %d agent positions to %s\n", count, path);
    fflush(stdout);
}

void ai_agent_load_positions(AiAgentSystem* sys, const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return;

    int count = 0;
    fread(&count, sizeof(int), 1, f);

    for (int i = 0; i < count; i++) {
        AgentSaveEntry entry;
        if (fread(&entry, sizeof(AgentSaveEntry), 1, f) != 1) break;

        // Find matching agent by ID
        for (int j = 0; j < sys->count; j++) {
            AiAgent* a = &sys->agents[j];
            if (strcmp(a->id, entry.id) == 0) {
                a->pos_d[0] = entry.pos[0];
                a->pos_d[1] = entry.pos[1];
                a->pos_d[2] = entry.pos[2];
                a->position = HMM_V3((float)entry.pos[0], (float)entry.pos[1], (float)entry.pos[2]);
                a->yaw = entry.yaw;
                float len = HMM_LenV3(a->position);
                if (len > 0.001f)
                    a->local_up = HMM_MulV3F(a->position, 1.0f / len);
                a->active = true;
                printf("[AGENT] Restored '%s' position from save\n", a->name);
                break;
            }
        }
    }

    fclose(f);
    fflush(stdout);
}

int ai_agent_nearest(const AiAgentSystem* sys, const double pos[3], float max_dist) {
    int best = -1;
    float best_dist = max_dist;
    for (int i = 0; i < sys->count; i++) {
        const AiAgent* a = &sys->agents[i];
        if (!a->active) continue;
        double dx = a->pos_d[0] - pos[0];
        double dy = a->pos_d[1] - pos[1];
        double dz = a->pos_d[2] - pos[2];
        float d = (float)sqrt(dx*dx + dy*dy + dz*dz);
        if (d < best_dist) { best_dist = d; best = i; }
    }
    return best;
}

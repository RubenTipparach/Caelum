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

    // Generate all parts into g_verts (single pass)
    g_vert_count = 0;

    // Primary: torso, shoulders, pelvis, legs
    hex_prism(0, torsoY, 0, torsoR, torsoH, stocky, 0.8f, pc[0], pc[1], pc[2]);
    hex_prism(-padSpacing, padY, 0, padR, padH, 0.8f, 0.8f, pc[0], pc[1], pc[2]);
    hex_prism(padSpacing, padY, 0, padR, padH, 0.8f, 0.8f, pc[0], pc[1], pc[2]);
    hex_prism(0, pelvisY, 0, pelvisR, pelvisH, stocky, 0.7f, pc[0], pc[1], pc[2]);
    hex_prism(-legSpacing, legY, 0, legR, legH, stocky*0.65f, 0.65f, pc[0], pc[1], pc[2]);
    hex_prism(legSpacing, legY, 0, legR, legH, stocky*0.65f, 0.65f, pc[0], pc[1], pc[2]);

    // Accent: neck, head, arms
    hex_prism(0, neckY, 0, neckR, neckH, 1.0f, 1.0f, ac[0], ac[1], ac[2]);
    hex_prism(0, headY, 0, headR, headH, 1.0f, 1.0f, ac[0], ac[1], ac[2]);
    hex_prism(-armSpacing, armY, 0, armR, armH, 0.55f, 0.55f, ac[0], ac[1], ac[2]);
    hex_prism(armSpacing, armY, 0, armR, armH, 0.55f, 0.55f, ac[0], ac[1], ac[2]);

    // Eyes
    float eyeR = headR * 0.22f;
    float eyeDepth = headR * 0.3f;
    float eyeSpacing = headR * 0.45f;
    float eyeY = headY + headH * 0.08f;
    float eyeZ = headR * 0.7f;
    hex_prism_z(-eyeSpacing, eyeY, eyeZ, eyeR, eyeDepth, 1,1,1);
    hex_prism_z(eyeSpacing, eyeY, eyeZ, eyeR, eyeDepth, 1,1,1);

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

    // Set up agent's AI with server connection (port 8080, no need to launch — already running)
    memset(&a->ai, 0, sizeof(AiNpc));
    a->ai.server_port = 8080;
    a->ai.server_running = true; // assume server launched by main ai_npc_init
    a->ai.provider = AI_PROVIDER_LOCAL;

    // Load config to check if we're using Claude
    {
        FILE* cf = fopen("config/ai_config.yaml", "r");
        if (cf) {
            char line[256];
            while (fgets(line, sizeof(line), cf)) {
                if (strstr(line, "provider:") && strstr(line, "claude")) {
                    a->ai.provider = AI_PROVIDER_CLAUDE;
                }
                char* api = strstr(line, "api_key:");
                if (api) {
                    api += 8;
                    while (*api == ' ' || *api == '"') api++;
                    int k = 0;
                    while (*api && *api != '"' && *api != '\n' && k < 127)
                        a->ai.claude_api_key[k++] = *api++;
                    a->ai.claude_api_key[k] = '\0';
                }
            }
            fclose(cf);
        }
    }

    // Build system prompt from directives + personality
    {
        char directives[512] = "";
        FILE* df = fopen("config/ai_directives.txt", "r");
        if (df) {
            int dlen = (int)fread(directives, 1, sizeof(directives) - 1, df);
            directives[dlen] = '\0';
            fclose(df);
        }
        snprintf(a->ai.system_prompt, sizeof(a->ai.system_prompt),
            "%s\n\nYour name is %s. %s\n"
            "Respond with dialogue and a list of actions in JSON format.",
            directives, a->name, a->directive);
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

        // Gravity: pull toward planet center
        a->position = HMM_V3((float)a->pos_d[0], (float)a->pos_d[1], (float)a->pos_d[2]);
        float r = HMM_LenV3(a->position);
        if (r > 0.001f) {
            a->local_up = HMM_MulV3F(a->position, 1.0f / r);
        }

        // Ground detection via hex terrain
        int gcol, grow, layer;
        hex_terrain_world_to_hex(ht, a->position, &gcol, &grow, &layer);
        float ground_r = hex_terrain_ground_height(ht, gcol, grow);

        if (ground_r > 0.0f && r <= ground_r + 0.1f) {
            // On ground — snap feet to surface (agent pos = feet + height offset)
            a->on_ground = true;
            float feet_r = ground_r + 0.05f; // tiny offset above ground
            a->pos_d[0] = a->local_up.X * feet_r;
            a->pos_d[1] = a->local_up.Y * feet_r;
            a->pos_d[2] = a->local_up.Z * feet_r;
            a->velocity = HMM_V3(0, 0, 0);
        } else if (ground_r > 0.0f) {
            // In air — apply gravity
            a->on_ground = false;
            float grav = 10.0f * dt;
            a->velocity = HMM_SubV3(a->velocity, HMM_MulV3F(a->local_up, grav));
            a->pos_d[0] += a->velocity.X * dt;
            a->pos_d[1] += a->velocity.Y * dt;
            a->pos_d[2] += a->velocity.Z * dt;
        }

        // ---- Pathfinding-based walking ----
        if (a->state == AGENT_STATE_WALKING && a->path.valid) {
            int next_q, next_r;
            if (ai_path_next(&a->path, &next_q, &next_r)) {
                // Get world position of next step
                float wx, wy, wz, ux, uy, uz;
                hex_terrain_hex_to_world(ht, next_q, next_r,
                    (int)((hex_terrain_ground_height(ht, next_q, next_r)
                           - hex_terrain_col_base_r(ht, next_q, next_r)) / HEX_HEIGHT),
                    &wx, &wy, &wz, &ux, &uy, &uz);

                // Move toward target
                HMM_Vec3 target = HMM_V3(wx, wy, wz);
                HMM_Vec3 to_target = HMM_SubV3(target, a->position);
                // Project onto tangent plane
                float dot = HMM_DotV3(to_target, a->local_up);
                to_target = HMM_SubV3(to_target, HMM_MulV3F(a->local_up, dot));
                float dist_to_step = HMM_LenV3(to_target);

                if (dist_to_step < 0.5f) {
                    // Reached this step, advance
                    ai_path_advance(&a->path);
                    if (ai_path_done(&a->path)) {
                        a->state = AGENT_STATE_IDLE;
                        a->has_target = false;
                    }
                } else {
                    // Walk toward step
                    HMM_Vec3 move_dir = HMM_MulV3F(to_target, 1.0f / dist_to_step);
                    float step = a->move_speed * dt;
                    if (step > dist_to_step) step = dist_to_step;
                    a->pos_d[0] += move_dir.X * step;
                    a->pos_d[1] += move_dir.Y * step;
                    a->pos_d[2] += move_dir.Z * step;
                    // Update facing
                    a->forward = move_dir;
                }
            } else {
                a->state = AGENT_STATE_IDLE;
            }
        }

        // ---- AI action execution ----
        ai_npc_poll(&a->ai);
        if (a->ai.state == AI_SUCCESS) {
            // Enqueue actions from LLM response
            ai_actions_enqueue(&a->executor, a->ai.actions, a->ai.action_count);
            a->ai.state = AI_IDLE;
        } else if (a->ai.state == AI_ERROR) {
            a->ai.state = AI_IDLE;
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
    printf("[AGENT] Saved %d agent positions to %s\n", count, path);
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

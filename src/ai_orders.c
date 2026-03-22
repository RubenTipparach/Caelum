#include "ai_orders.h"
#include "planet.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <direct.h>
#endif
#include <ctype.h>

// Direction names to hex direction index
// Hex dirs: E=0, NE=1, NW=2, W=3, SW=4, SE=5
static int parse_direction(const char* s) {
    if (strcmp(s, "e") == 0 || strcmp(s, "east") == 0)       return 0;
    if (strcmp(s, "ne") == 0 || strcmp(s, "northeast") == 0) return 1;
    if (strcmp(s, "nw") == 0 || strcmp(s, "northwest") == 0) return 2;
    if (strcmp(s, "n") == 0 || strcmp(s, "north") == 0)      return 1; // NE as proxy for N
    if (strcmp(s, "w") == 0 || strcmp(s, "west") == 0)       return 3;
    if (strcmp(s, "sw") == 0 || strcmp(s, "southwest") == 0) return 4;
    if (strcmp(s, "se") == 0 || strcmp(s, "southeast") == 0) return 5;
    if (strcmp(s, "s") == 0 || strcmp(s, "south") == 0)      return 4; // SW as proxy for S
    if (strcmp(s, "forward") == 0 || strcmp(s, "ahead") == 0) return -1;
    return -2; // unknown
}

// Hex neighbor step
static void hex_step(int q, int r, int dir, int* nq, int* nr) {
    static const int even_dq[6] = {+1, +1,  0, -1, -1,  0};
    static const int even_dr[6] = { 0, -1, -1, -1,  0, +1};
    static const int odd_dq[6]  = {+1, +1,  0, -1, -1,  0};
    static const int odd_dr[6]  = {+1,  0, -1,  0, +1, +1};
    if (q % 2 == 0) {
        *nq = q + even_dq[dir % 6];
        *nr = r + even_dr[dir % 6];
    } else {
        *nq = q + odd_dq[dir % 6];
        *nr = r + odd_dr[dir % 6];
    }
}

static void to_lower(const char* src, char* dst, int max) {
    int i = 0;
    for (; src[i] && i < max - 1; i++)
        dst[i] = (src[i] >= 'A' && src[i] <= 'Z') ? src[i] + 32 : src[i];
    dst[i] = '\0';
}

bool ai_order_parse(const char* text, AiOrder* out) {
    memset(out, 0, sizeof(AiOrder));
    snprintf(out->block, sizeof(out->block), "stone");

    char lower[256];
    to_lower(text, lower, sizeof(lower));

    // "follow me" / "follow" / "come with me" / "come here"
    if (strstr(lower, "follow") || strstr(lower, "come with") ||
        strstr(lower, "come here") || strstr(lower, "come along")) {
        out->type = ORDER_FOLLOW;
        return true;
    }

    // Stop commands — "stay", "stop", "wait here", "stop building", "stop following", etc.
    if (strstr(lower, "stay") || strstr(lower, "wait here") || strstr(lower, "hold") ||
        strstr(lower, "stop building") || strstr(lower, "stop following") ||
        strstr(lower, "stop construct") || strstr(lower, "stop work") ||
        strstr(lower, "stop moving") || strstr(lower, "stop walk") ||
        strstr(lower, "cancel") || strstr(lower, "abort") || strstr(lower, "nevermind") ||
        strstr(lower, "knock it off") || strstr(lower, "that's enough") ||
        strcmp(lower, "stop") == 0) {
        out->type = ORDER_STAY;
        return true;
    }

    // "go <dir> <N>" / "walk <dir> <N>" / "move <dir> <N>"
    if (strstr(lower, "go ") == lower || strstr(lower, "walk ") == lower ||
        strstr(lower, "move ") == lower) {
        out->type = ORDER_GO;
        // Parse: command dir count
        char cmd[16], dir_str[16];
        int count = 5;
        int parsed = sscanf(lower, "%15s %15s %d", cmd, dir_str, &count);
        if (parsed >= 2) {
            out->direction = parse_direction(dir_str);
            if (out->direction == -2) return false;
            out->count = count > 0 ? count : 5;
            return true;
        }
        return false;
    }

    // "build wall <dir> <N>" / "build stairs <dir> <N>"
    if (strstr(lower, "build ") == lower) {
        char cmd[16], what[16], dir_str[16];
        int count = 3;
        int parsed = sscanf(lower, "%15s %15s %15s %d", cmd, what, dir_str, &count);
        if (parsed >= 3) {
            if (strcmp(what, "wall") == 0) out->type = ORDER_BUILD_WALL;
            else if (strcmp(what, "stairs") == 0 || strcmp(what, "staircase") == 0)
                out->type = ORDER_BUILD_STAIRS;
            else return false;

            out->direction = parse_direction(dir_str);
            if (out->direction == -2) return false;
            out->count = count > 0 ? count : 3;

            // Check for block type: "build wall north 5 dirt"
            char block_str[16] = "";
            if (sscanf(lower, "%*s %*s %*s %*d %15s", block_str) == 1 && block_str[0]) {
                snprintf(out->block, sizeof(out->block), "%s", block_str);
            }
            return true;
        }
        return false;
    }

    // "clear area <N>" / "clear ahead <N>" / "clear <dir> <N>"
    if (strstr(lower, "clear ") == lower) {
        char cmd[16], what[16];
        int count = 3;
        int parsed = sscanf(lower, "%15s %15s %d", cmd, what, &count);
        if (parsed >= 2) {
            if (strcmp(what, "area") == 0 || strcmp(what, "around") == 0) {
                out->type = ORDER_CLEAR_AREA;
                out->count = count > 0 ? count : 3;
                return true;
            }
            if (strcmp(what, "ahead") == 0 || strcmp(what, "forward") == 0) {
                out->type = ORDER_CLEAR_AHEAD;
                out->direction = -1;
                out->count = count > 0 ? count : 5;
                return true;
            }
            // "clear <dir> <N>"
            int dir = parse_direction(what);
            if (dir >= 0) {
                out->type = ORDER_CLEAR_AHEAD;
                out->direction = dir;
                out->count = count > 0 ? count : 5;
                return true;
            }
        }
        return false;
    }

    return false;
}

static const char* dir_name(int dir) {
    static const char* names[] = {"E", "NE", "NW", "W", "SW", "SE"};
    if (dir >= 0 && dir < 6) return names[dir];
    return "forward";
}

bool ai_order_to_script(const AiOrder* order, const HexTerrain* ht,
                        int agent_q, int agent_r, int agent_ground_layer,
                        AiScript* out) {
    memset(out, 0, sizeof(AiScript));

    int dir = order->direction;
    if (dir < 0) dir = 0; // default to E if forward isn't resolved

    switch (order->type) {
        case ORDER_GO: {
            snprintf(out->name, sizeof(out->name), "go_%s_%d",
                     dir_name(dir), order->count);
            snprintf(out->description, sizeof(out->description),
                     "Walk %d cells %s", order->count, dir_name(dir));

            // Generate move steps along direction
            int q = agent_q, r = agent_r;
            for (int i = 0; i < order->count && out->step_count < AI_SCRIPT_MAX_STEPS; i++) {
                int nq, nr;
                hex_step(q, r, dir, &nq, &nr);
                AiScriptStep* s = &out->steps[out->step_count++];
                s->action = SCRIPT_ACT_MOVE_TO;
                s->q = nq; s->r = nr;
                q = nq; r = nr;
            }
            // Final move_to to the last cell
            if (out->step_count > 0) {
                // Simplify: just one move_to to the final destination
                AiScriptStep final_step = out->steps[out->step_count - 1];
                out->step_count = 1;
                out->steps[0] = final_step;
            }
            return out->step_count > 0;
        }

        case ORDER_BUILD_WALL: {
            snprintf(out->name, sizeof(out->name), "wall_%s_%d",
                     dir_name(dir), order->count);
            snprintf(out->description, sizeof(out->description),
                     "Build %d-block %s wall going %s", order->count, order->block, dir_name(dir));

            int q = agent_q, r = agent_r;
            float base_r = hex_terrain_col_base_r(ht, q, r);
            int base_layer = (base_r > 0) ? agent_ground_layer : 0;

            for (int i = 0; i < order->count && out->step_count < AI_SCRIPT_MAX_STEPS - 2; i++) {
                int nq, nr;
                hex_step(q, r, dir, &nq, &nr);

                // Move to the cell
                AiScriptStep* move = &out->steps[out->step_count++];
                move->action = SCRIPT_ACT_MOVE_TO;
                move->q = nq; move->r = nr;

                // Place wall block at ground+1 (one block high wall)
                AiScriptStep* place = &out->steps[out->step_count++];
                place->action = SCRIPT_ACT_PLACE;
                place->q = nq; place->r = nr;
                place->layer = base_layer + 1;
                snprintf(place->block, sizeof(place->block), "%s", order->block);

                q = nq; r = nr;
            }

            // Say done
            AiScriptStep* say = &out->steps[out->step_count++];
            say->action = SCRIPT_ACT_SAY;
            snprintf(say->text, sizeof(say->text), "Wall's done! %d blocks placed.", order->count);

            return out->step_count > 0;
        }

        case ORDER_BUILD_STAIRS: {
            snprintf(out->name, sizeof(out->name), "stairs_%s_%d",
                     dir_name(dir), order->count);
            snprintf(out->description, sizeof(out->description),
                     "Build %d-step staircase going %s", order->count, dir_name(dir));

            int q = agent_q, r = agent_r;
            float base_r = hex_terrain_col_base_r(ht, q, r);
            int base_layer = (base_r > 0) ? agent_ground_layer : 0;

            for (int i = 0; i < order->count && out->step_count < AI_SCRIPT_MAX_STEPS - 2; i++) {
                int nq, nr;
                hex_step(q, r, dir, &nq, &nr);

                // Move to the cell
                AiScriptStep* move = &out->steps[out->step_count++];
                move->action = SCRIPT_ACT_MOVE_TO;
                move->q = nq; move->r = nr;

                // Place step block (each one higher than the last)
                AiScriptStep* place = &out->steps[out->step_count++];
                place->action = SCRIPT_ACT_PLACE;
                place->q = nq; place->r = nr;
                place->layer = base_layer + i + 1;
                snprintf(place->block, sizeof(place->block), "%s", order->block);

                q = nq; r = nr;
            }

            AiScriptStep* say = &out->steps[out->step_count++];
            say->action = SCRIPT_ACT_SAY;
            snprintf(say->text, sizeof(say->text), "Staircase done! %d steps up.", order->count);

            return out->step_count > 0;
        }

        case ORDER_CLEAR_AREA: {
            snprintf(out->name, sizeof(out->name), "clear_area_%d", order->count);
            snprintf(out->description, sizeof(out->description),
                     "Clear blocks in %d-cell radius", order->count);

            int radius = order->count;
            for (int dq = -radius; dq <= radius && out->step_count < AI_SCRIPT_MAX_STEPS - 1; dq++) {
                for (int dr = -radius; dr <= radius && out->step_count < AI_SCRIPT_MAX_STEPS - 1; dr++) {
                    int q = agent_q + dq, r = agent_r + dr;
                    float gr = hex_terrain_ground_height(ht, q, r);
                    if (gr <= 0.0f) continue;

                    float base = hex_terrain_col_base_r(ht, q, r);
                    int top = (base > 0) ? (int)((gr - base) / HEX_HEIGHT) : 0;

                    // Only clear blocks above the agent's ground level
                    if (top > agent_ground_layer) {
                        AiScriptStep* brk = &out->steps[out->step_count++];
                        brk->action = SCRIPT_ACT_BREAK;
                        brk->q = q; brk->r = r;
                        brk->layer = top;
                    }
                }
            }

            AiScriptStep* say = &out->steps[out->step_count++];
            say->action = SCRIPT_ACT_SAY;
            snprintf(say->text, sizeof(say->text), "Area cleared! Broke %d blocks.", out->step_count - 1);

            return out->step_count > 0;
        }

        case ORDER_CLEAR_AHEAD: {
            snprintf(out->name, sizeof(out->name), "clear_%s_%d",
                     dir_name(dir), order->count);
            snprintf(out->description, sizeof(out->description),
                     "Clear blocks %s for %d cells", dir_name(dir), order->count);

            int q = agent_q, r = agent_r;
            for (int i = 0; i < order->count && out->step_count < AI_SCRIPT_MAX_STEPS - 1; i++) {
                int nq, nr;
                hex_step(q, r, dir, &nq, &nr);

                float gr = hex_terrain_ground_height(ht, nq, nr);
                if (gr <= 0.0f) break;

                float base = hex_terrain_col_base_r(ht, nq, nr);
                int top = (base > 0) ? (int)((gr - base) / HEX_HEIGHT) : 0;

                // Break blocks above agent level to clear a path
                for (int l = agent_ground_layer + 1; l <= top + 2 &&
                     out->step_count < AI_SCRIPT_MAX_STEPS - 1; l++) {
                    uint8_t voxel = hex_terrain_get_voxel(ht, nq, nr, l);
                    if (voxel != VOXEL_AIR && voxel != VOXEL_BEDROCK) {
                        AiScriptStep* brk = &out->steps[out->step_count++];
                        brk->action = SCRIPT_ACT_BREAK;
                        brk->q = nq; brk->r = nr;
                        brk->layer = l;
                    }
                }
                q = nq; r = nr;
            }

            AiScriptStep* say = &out->steps[out->step_count++];
            say->action = SCRIPT_ACT_SAY;
            snprintf(say->text, sizeof(say->text), "Path cleared!");

            return out->step_count > 0;
        }

        default:
            return false;
    }
}

void ai_order_save(const AiOrder* order, const AiScript* script,
                   AiMemory* mem) {
    // Save script JSON to agent's learning folder
    char path[512];
    snprintf(path, sizeof(path), "%s/orders", mem->agent_dir);

#ifdef _WIN32
    _mkdir(path);
#else
    mkdir(path, 0755);
#endif

    snprintf(path, sizeof(path), "%s/orders/%s.json", mem->agent_dir, script->name);
    ai_script_save(script, path);

    // Journal the order
    char entry[256];
    snprintf(entry, sizeof(entry), "Received order: %s (%d steps)",
             script->description, script->step_count);
    ai_memory_journal(mem, entry);
}

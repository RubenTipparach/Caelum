#include "ai_pathfind.h"
#include "planet.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Hex neighbor offsets (flat-top, offset coordinates matching hex_terrain)
// Even columns: neighbors at (q±1,r), (q,r±1), (q±1,r-1)
// Odd columns:  neighbors at (q±1,r), (q,r±1), (q±1,r+1)
static void hex_neighbors(int q, int r, int out_q[6], int out_r[6]) {
    if (q % 2 == 0) {
        // Even column
        out_q[0] = q+1; out_r[0] = r;
        out_q[1] = q+1; out_r[1] = r-1;
        out_q[2] = q;   out_r[2] = r-1;
        out_q[3] = q-1; out_r[3] = r-1;
        out_q[4] = q-1; out_r[4] = r;
        out_q[5] = q;   out_r[5] = r+1;
    } else {
        // Odd column
        out_q[0] = q+1; out_r[0] = r+1;
        out_q[1] = q+1; out_r[1] = r;
        out_q[2] = q;   out_r[2] = r-1;
        out_q[3] = q-1; out_r[3] = r;
        out_q[4] = q-1; out_r[4] = r+1;
        out_q[5] = q;   out_r[5] = r+1;
    }
}

// Simple hex distance heuristic
static int hex_dist(int q1, int r1, int q2, int r2) {
    int dq = abs(q2 - q1);
    int dr = abs(r2 - r1);
    return dq + (dr > dq / 2 ? dr - dq / 2 : 0);
}

// Check if a hex cell is walkable: must have ground and 2 layers of air above
static bool is_walkable(const HexTerrain* ht, int q, int r) {
    float ground = hex_terrain_ground_height(ht, q, r);
    if (ground <= 0.0f) return false; // no chunk loaded or no ground

    // Check headroom (need 2 layers of air for the agent to fit)
    float base_r = hex_terrain_col_base_r(ht, q, r);
    if (base_r <= 0.0f) return false;
    int layer = (int)((ground - base_r) / HEX_HEIGHT);
    return hex_terrain_has_headroom(ht, q, r, layer + 1, 2);
}

// ---- A* implementation ----

#define ASTAR_MAX_NODES 4096

typedef struct {
    int q, r;
    int parent;     // index into nodes array
    float g, f;     // cost from start, estimated total cost
    bool closed;
} AstarNode;

// Simple open set using linear scan (fine for <4K nodes)
static AstarNode g_nodes[ASTAR_MAX_NODES];
static int g_node_count;

static int find_node(int q, int r) {
    for (int i = 0; i < g_node_count; i++) {
        if (g_nodes[i].q == q && g_nodes[i].r == r) return i;
    }
    return -1;
}

static int lowest_f_open(void) {
    int best = -1;
    float best_f = 1e9f;
    for (int i = 0; i < g_node_count; i++) {
        if (!g_nodes[i].closed && g_nodes[i].f < best_f) {
            best_f = g_nodes[i].f;
            best = i;
        }
    }
    return best;
}

bool ai_pathfind(AiPath* path, const HexTerrain* ht,
                 int start_q, int start_r,
                 int end_q, int end_r) {
    memset(path, 0, sizeof(AiPath));

    // Quick check: is destination reachable?
    if (!is_walkable(ht, end_q, end_r)) return false;

    g_node_count = 0;

    // Add start node
    g_nodes[g_node_count++] = (AstarNode){
        .q = start_q, .r = start_r,
        .parent = -1, .g = 0,
        .f = (float)hex_dist(start_q, start_r, end_q, end_r),
    };

    int iterations = 0;
    while (iterations++ < ASTAR_MAX_NODES) {
        int current = lowest_f_open();
        if (current < 0) break; // no path

        AstarNode* cur = &g_nodes[current];

        // Reached goal?
        if (cur->q == end_q && cur->r == end_r) {
            // Reconstruct path
            int trace[AI_PATH_MAX_STEPS];
            int trace_len = 0;
            int idx = current;
            while (idx >= 0 && trace_len < AI_PATH_MAX_STEPS) {
                trace[trace_len++] = idx;
                idx = g_nodes[idx].parent;
            }
            // Reverse into path
            path->step_count = trace_len;
            for (int i = 0; i < trace_len; i++) {
                int j = trace_len - 1 - i;
                path->steps[i].q = g_nodes[trace[j]].q;
                path->steps[i].r = g_nodes[trace[j]].r;
            }
            path->current_step = 0;
            path->valid = true;
            return true;
        }

        cur->closed = true;

        // Expand neighbors
        int nq[6], nr[6];
        hex_neighbors(cur->q, cur->r, nq, nr);

        for (int i = 0; i < 6; i++) {
            if (!is_walkable(ht, nq[i], nr[i])) continue;

            float step_cost = 1.0f;
            // Penalize height changes (climbing is harder)
            float h1 = hex_terrain_ground_height(ht, cur->q, cur->r);
            float h2 = hex_terrain_ground_height(ht, nq[i], nr[i]);
            float dh = fabsf(h2 - h1);
            if (dh > 2.0f * HEX_HEIGHT) continue; // too steep, can't walk
            step_cost += dh * 0.5f;

            float new_g = cur->g + step_cost;

            int existing = find_node(nq[i], nr[i]);
            if (existing >= 0) {
                if (g_nodes[existing].closed) continue;
                if (new_g < g_nodes[existing].g) {
                    g_nodes[existing].g = new_g;
                    g_nodes[existing].f = new_g + (float)hex_dist(nq[i], nr[i], end_q, end_r);
                    g_nodes[existing].parent = current;
                }
            } else if (g_node_count < ASTAR_MAX_NODES) {
                g_nodes[g_node_count++] = (AstarNode){
                    .q = nq[i], .r = nr[i],
                    .parent = current,
                    .g = new_g,
                    .f = new_g + (float)hex_dist(nq[i], nr[i], end_q, end_r),
                };
            }
        }
    }

    return false; // no path found
}

bool ai_path_next(const AiPath* path, int* out_q, int* out_r) {
    if (!path->valid || path->current_step >= path->step_count) return false;
    *out_q = path->steps[path->current_step].q;
    *out_r = path->steps[path->current_step].r;
    return true;
}

void ai_path_advance(AiPath* path) {
    if (path->valid && path->current_step < path->step_count)
        path->current_step++;
}

bool ai_path_done(const AiPath* path) {
    return !path->valid || path->current_step >= path->step_count;
}

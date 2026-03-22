#include "ai_sensors.h"
#include "planet.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// Hex neighbor offsets (offset coordinates)
static void hex_neighbor(int q, int r, int dir, int* nq, int* nr) {
    static const int even_dq[6] = {+1, +1,  0, -1, -1,  0};
    static const int even_dr[6] = { 0, -1, -1, -1,  0, +1};
    static const int odd_dq[6]  = {+1, +1,  0, -1, -1,  0};
    static const int odd_dr[6]  = {+1,  0, -1,  0, +1, +1};
    if (q % 2 == 0) {
        *nq = q + even_dq[dir];
        *nr = r + even_dr[dir];
    } else {
        *nq = q + odd_dq[dir];
        *nr = r + odd_dr[dir];
    }
}

static const char* voxel_name(uint8_t type) {
    switch (type) {
        case VOXEL_AIR:       return "air";
        case VOXEL_WATER:     return "water";
        case VOXEL_SAND:      return "sand";
        case VOXEL_DIRT:      return "dirt";
        case VOXEL_GRASS:     return "grass";
        case VOXEL_STONE:     return "stone";
        case VOXEL_ICE:       return "ice";
        case VOXEL_BEDROCK:   return "bedrock";
        case VOXEL_TORCH:     return "torch";
        case VOXEL_MOON_ROCK: return "moon_rock";
        default:              return "unknown";
    }
}

static int hex_distance(int q1, int r1, int q2, int r2) {
    int dq = abs(q2 - q1);
    int dr = abs(r2 - r1);
    return dq + (dr > dq / 2 ? dr - dq / 2 : 0);
}

void ai_sensors_scan(const HexTerrain* ht, int center_q, int center_r,
                     int radius, AiScanResult* out) {
    memset(out, 0, sizeof(AiScanResult));
    out->agent_q = center_q;
    out->agent_r = center_r;
    out->agent_ground_r = hex_terrain_ground_height(ht, center_q, center_r);

    float base = hex_terrain_col_base_r(ht, center_q, center_r);
    if (base > 0.0f)
        out->agent_ground_layer = (int)((out->agent_ground_r - base) / HEX_HEIGHT);

    // Spiral scan outward from center
    for (int dq = -radius; dq <= radius && out->count < AI_SCAN_MAX_CELLS; dq++) {
        for (int dr = -radius; dr <= radius && out->count < AI_SCAN_MAX_CELLS; dr++) {
            int q = center_q + dq;
            int r = center_r + dr;
            int dist = hex_distance(center_q, center_r, q, r);
            if (dist > radius) continue;

            float gr = hex_terrain_ground_height(ht, q, r);
            if (gr <= 0.0f) continue; // chunk not loaded

            float col_base = hex_terrain_col_base_r(ht, q, r);
            int layer = (col_base > 0.0f) ? (int)((gr - col_base) / HEX_HEIGHT) : 0;

            // Get surface block type
            uint8_t surface = hex_terrain_get_voxel(ht, q, r, layer);
            if (surface == VOXEL_AIR && layer > 0)
                surface = hex_terrain_get_voxel(ht, q, r, layer - 1);

            // Check walkability
            bool walkable = (gr > 0.0f) && hex_terrain_has_headroom(ht, q, r, layer + 1, 2);

            AiHexObs* obs = &out->cells[out->count++];
            obs->q = q;
            obs->r = r;
            obs->ground_layer = layer;
            obs->ground_r = gr;
            obs->surface_type = surface;
            obs->height_diff = gr - out->agent_ground_r;
            obs->distance = (float)dist;
            obs->walkable = walkable;
        }
    }
}

int ai_sensors_report(const AiScanResult* scan, char* buf, int buf_size) {
    int pos = 0;

    pos += snprintf(buf + pos, buf_size - pos,
        "=== SURROUNDINGS (you are at q=%d, r=%d, ground_layer=%d) ===\n",
        scan->agent_q, scan->agent_r, scan->agent_ground_layer);

    // Summarize terrain composition
    int counts[VOXEL_TYPE_COUNT] = {0};
    int walkable_count = 0, blocked_count = 0;
    float min_h = 1e9f, max_h = -1e9f;
    for (int i = 0; i < scan->count; i++) {
        const AiHexObs* c = &scan->cells[i];
        if (c->surface_type < VOXEL_TYPE_COUNT) counts[c->surface_type]++;
        if (c->walkable) walkable_count++; else blocked_count++;
        if (c->height_diff < min_h) min_h = c->height_diff;
        if (c->height_diff > max_h) max_h = c->height_diff;
    }

    pos += snprintf(buf + pos, buf_size - pos,
        "Terrain: %d cells scanned, %d walkable, %d blocked\n",
        scan->count, walkable_count, blocked_count);
    pos += snprintf(buf + pos, buf_size - pos,
        "Elevation range: %.1fm below to %.1fm above you\n", -min_h, max_h);

    // Surface types present
    pos += snprintf(buf + pos, buf_size - pos, "Surface types: ");
    bool first = true;
    for (int t = 1; t < VOXEL_TYPE_COUNT; t++) {
        if (counts[t] > 0) {
            pos += snprintf(buf + pos, buf_size - pos, "%s%s(%d)",
                           first ? "" : ", ", voxel_name((uint8_t)t), counts[t]);
            first = false;
        }
    }
    pos += snprintf(buf + pos, buf_size - pos, "\n");

    // Notable features: cliffs, water, obstacles
    pos += snprintf(buf + pos, buf_size - pos, "\nNearby cells (q,r | type | height_diff | walkable):\n");
    for (int i = 0; i < scan->count && pos < buf_size - 80; i++) {
        const AiHexObs* c = &scan->cells[i];
        if (c->distance <= 3.0f || !c->walkable ||
            fabsf(c->height_diff) > 1.0f) {
            pos += snprintf(buf + pos, buf_size - pos,
                "  (%d,%d) %s %+.1fm %s\n",
                c->q, c->r, voxel_name(c->surface_type),
                c->height_diff, c->walkable ? "ok" : "BLOCKED");
        }
    }

    return pos;
}

int ai_sensors_column_report(const HexTerrain* ht, int q, int r,
                             char* buf, int buf_size) {
    int pos = 0;
    float base_r = hex_terrain_col_base_r(ht, q, r);
    if (base_r <= 0.0f) {
        return snprintf(buf, buf_size, "Column (%d,%d): not loaded\n", q, r);
    }

    float ground = hex_terrain_ground_height(ht, q, r);
    int ground_layer = (int)((ground - base_r) / HEX_HEIGHT);

    pos += snprintf(buf + pos, buf_size - pos,
        "Column (%d,%d): ground_layer=%d\n", q, r, ground_layer);

    // Show voxels around ground level
    int start = ground_layer > 3 ? ground_layer - 3 : 0;
    int end = ground_layer + 5;
    if (end > MAX_VOXEL_HEIGHT) end = MAX_VOXEL_HEIGHT;

    for (int l = end - 1; l >= start && pos < buf_size - 40; l--) {
        uint8_t v = hex_terrain_get_voxel(ht, q, r, l);
        const char* marker = (l == ground_layer) ? " <-- ground" : "";
        pos += snprintf(buf + pos, buf_size - pos,
            "  layer %d: %s%s\n", l, voxel_name(v), marker);
    }

    return pos;
}

int ai_sensors_directional_report(const HexTerrain* ht,
                                  int center_q, int center_r,
                                  int depth, char* buf, int buf_size) {
    int pos = 0;
    static const char* dir_names[6] = {"E", "NE", "NW", "W", "SW", "SE"};

    float center_gr = hex_terrain_ground_height(ht, center_q, center_r);

    pos += snprintf(buf + pos, buf_size - pos,
        "=== DIRECTIONS from (%d,%d) ===\n", center_q, center_r);

    for (int d = 0; d < 6 && pos < buf_size - 100; d++) {
        pos += snprintf(buf + pos, buf_size - pos, "%s: ", dir_names[d]);
        int q = center_q, r = center_r;
        for (int step = 0; step < depth && pos < buf_size - 40; step++) {
            int nq, nr;
            hex_neighbor(q, r, d, &nq, &nr);
            float gr = hex_terrain_ground_height(ht, nq, nr);
            if (gr <= 0.0f) {
                pos += snprintf(buf + pos, buf_size - pos, "unloaded ");
                break;
            }
            float dh = gr - center_gr;
            uint8_t surf = hex_terrain_get_voxel(ht, nq, nr,
                (int)((gr - hex_terrain_col_base_r(ht, nq, nr)) / HEX_HEIGHT));
            if (surf == VOXEL_AIR) surf = hex_terrain_get_voxel(ht, nq, nr,
                (int)((gr - hex_terrain_col_base_r(ht, nq, nr)) / HEX_HEIGHT) - 1);

            bool walk = hex_terrain_has_headroom(ht, nq, nr,
                (int)((gr - hex_terrain_col_base_r(ht, nq, nr)) / HEX_HEIGHT) + 1, 2);

            pos += snprintf(buf + pos, buf_size - pos, "%s(%+.0f%s) ",
                           voxel_name(surf), dh, walk ? "" : "!");
            q = nq; r = nr;
        }
        pos += snprintf(buf + pos, buf_size - pos, "\n");
    }

    return pos;
}

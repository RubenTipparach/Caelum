#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>

typedef enum { PARAM_FLOAT, PARAM_VEC3 } ParamType;

typedef struct {
    const char* name;
    size_t offset;
    ParamType type;
} ParamEntry;

static const ParamEntry param_table[] = {
    { "atmosphere_height", offsetof(VisualConfig, atmosphere_height), PARAM_FLOAT },
    { "rayleigh_scale",    offsetof(VisualConfig, rayleigh_scale),    PARAM_FLOAT },
    { "mie_scale",         offsetof(VisualConfig, mie_scale),         PARAM_FLOAT },
    { "mie_g",             offsetof(VisualConfig, mie_g),             PARAM_FLOAT },
    { "sun_intensity",     offsetof(VisualConfig, sun_intensity),     PARAM_FLOAT },
    { "scale_height",      offsetof(VisualConfig, scale_height),      PARAM_FLOAT },
    { "fog_scale_height",  offsetof(VisualConfig, fog_scale_height),  PARAM_FLOAT },
    { "dusk_sun_color",    offsetof(VisualConfig, dusk_sun_color),    PARAM_VEC3  },
    { "day_sun_color",     offsetof(VisualConfig, day_sun_color),     PARAM_VEC3  },
    { "hex_fade_start",    offsetof(VisualConfig, hex_fade_start),    PARAM_FLOAT },
    { "hex_fade_end",      offsetof(VisualConfig, hex_fade_end),      PARAM_FLOAT },
    { "lod_split_factor",  offsetof(VisualConfig, lod_split_factor),  PARAM_FLOAT },
};

#define PARAM_COUNT (sizeof(param_table) / sizeof(param_table[0]))

VisualConfig config_defaults(void) {
    return (VisualConfig){
        .atmosphere_height = 50000.0f,
        .rayleigh_scale    = 0.015f,
        .mie_scale         = 0.01f,
        .mie_g             = 0.85f,
        .sun_intensity     = 5.0f,
        .scale_height      = 0.25f,
        .fog_scale_height  = 0.25f,
        .dusk_sun_color    = {{ 1.3f, 0.45f, 0.12f }},
        .day_sun_color     = {{ 1.0f, 0.98f, 0.95f }},
        .hex_fade_start    = 600.0f,
        .hex_fade_end      = 1200.0f,
        .lod_split_factor  = 8.0f,
    };
}

// Trim leading whitespace in-place, return pointer to first non-space
static char* trim_left(char* s) {
    while (*s && isspace((unsigned char)*s)) s++;
    return s;
}

// Trim trailing whitespace in-place
static void trim_right(char* s) {
    char* end = s + strlen(s) - 1;
    while (end >= s && isspace((unsigned char)*end)) *end-- = '\0';
}

bool config_load(VisualConfig* cfg, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) {
        printf("[CONFIG] Could not open '%s', using defaults\n", path);
        fflush(stdout);
        return false;
    }

    char line[256];
    int line_num = 0;
    int params_loaded = 0;

    while (fgets(line, sizeof(line), f)) {
        line_num++;

        // Strip comment
        char* hash = strchr(line, '#');
        if (hash) *hash = '\0';

        // Trim
        char* trimmed = trim_left(line);
        trim_right(trimmed);
        if (*trimmed == '\0') continue;

        // Split on first ':'
        char* colon = strchr(trimmed, ':');
        if (!colon) {
            printf("[CONFIG] line %d: no ':' found, skipping\n", line_num);
            continue;
        }

        *colon = '\0';
        char* key = trimmed;
        char* val = trim_left(colon + 1);
        trim_right(key);

        // Find matching parameter
        bool found = false;
        for (size_t i = 0; i < PARAM_COUNT; i++) {
            if (strcmp(key, param_table[i].name) != 0) continue;
            found = true;

            char* base = (char*)cfg + param_table[i].offset;

            if (param_table[i].type == PARAM_FLOAT) {
                float v = strtof(val, NULL);
                *(float*)base = v;
            } else { // PARAM_VEC3
                float x = 0, y = 0, z = 0;
                if (sscanf(val, "%f , %f , %f", &x, &y, &z) == 3) {
                    HMM_Vec3* vec = (HMM_Vec3*)base;
                    vec->X = x; vec->Y = y; vec->Z = z;
                } else {
                    printf("[CONFIG] line %d: bad vec3 format for '%s' (expected: x, y, z)\n",
                           line_num, key);
                }
            }
            params_loaded++;
            break;
        }

        if (!found) {
            printf("[CONFIG] line %d: unknown key '%s'\n", line_num, key);
        }
    }

    fclose(f);
    printf("[CONFIG] Loaded %d parameters from '%s'\n", params_loaded, path);
    fflush(stdout);
    return true;
}

bool config_save(const VisualConfig* cfg, const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) {
        printf("[CONFIG] Could not write '%s'\n", path);
        return false;
    }

    fprintf(f, "# config.yaml - Hex Planets visual config\n");
    fprintf(f, "# Press R in-game to hot-reload\n\n");

    fprintf(f, "# === Atmosphere (fullscreen ray march) ===\n");
    fprintf(f, "atmosphere_height: %.0f\n", cfg->atmosphere_height);
    fprintf(f, "rayleigh_scale: %.4f\n", cfg->rayleigh_scale);
    fprintf(f, "mie_scale: %.4f\n", cfg->mie_scale);
    fprintf(f, "mie_g: %.2f\n", cfg->mie_g);
    fprintf(f, "sun_intensity: %.1f\n", cfg->sun_intensity);
    fprintf(f, "scale_height: %.2f\n", cfg->scale_height);

    fprintf(f, "\n# === Terrain fog (aerial perspective) ===\n");
    fprintf(f, "fog_scale_height: %.2f\n", cfg->fog_scale_height);

    fprintf(f, "\n# === Dusk tinting ===\n");
    fprintf(f, "dusk_sun_color: %.2f, %.2f, %.2f\n",
            cfg->dusk_sun_color.X, cfg->dusk_sun_color.Y, cfg->dusk_sun_color.Z);
    fprintf(f, "day_sun_color: %.2f, %.2f, %.2f\n",
            cfg->day_sun_color.X, cfg->day_sun_color.Y, cfg->day_sun_color.Z);

    fprintf(f, "\n# === Hex terrain distance fade ===\n");
    fprintf(f, "hex_fade_start: %.0f\n", cfg->hex_fade_start);
    fprintf(f, "hex_fade_end: %.0f\n", cfg->hex_fade_end);

    fprintf(f, "\n# === LOD ===\n");
    fprintf(f, "lod_split_factor: %.1f\n", cfg->lod_split_factor);

    fclose(f);
    printf("[CONFIG] Saved config to '%s'\n", path);
    fflush(stdout);
    return true;
}

#ifndef CONFIG_H
#define CONFIG_H

#include "HandmadeMath.h"
#include <stdbool.h>

typedef struct VisualConfig {
    // Atmosphere (fullscreen ray march)
    float atmosphere_height;    // offset above surface (meters)
    float rayleigh_scale;       // Rayleigh scattering coefficient
    float mie_scale;            // Mie scattering coefficient
    float mie_g;                // Mie phase asymmetry (0=isotropic, 1=forward)
    float sun_intensity;        // sun brightness multiplier
    float scale_height;         // density falloff (0-1)

    // Terrain fog (aerial perspective)
    float fog_scale_height;     // should match scale_height for consistency

    // Dusk tinting
    HMM_Vec3 dusk_sun_color;    // warm sunset color
    HMM_Vec3 day_sun_color;     // neutral daylight color

    // Hex terrain distance fade
    float hex_fade_start;       // distance where texture begins fading (meters)
    float hex_fade_end;         // distance where fully flat color (meters)

    // LOD
    float lod_split_factor;     // split threshold multiplier (higher = more detail)
} VisualConfig;

// Returns hardcoded defaults (same values currently in the code)
VisualConfig config_defaults(void);

// Load config from file. Returns true on success.
// On failure, logs a warning and leaves cfg unchanged.
bool config_load(VisualConfig* cfg, const char* path);

// Write current config to file (for generating a starter config.yaml)
bool config_save(const VisualConfig* cfg, const char* path);

#endif

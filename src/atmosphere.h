#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "HandmadeMath.h"
#include "sokol_gfx.h"

typedef struct AtmosphereConfig {
    float planet_radius;        // Effective surface radius (planet.radius + sea_level)
    float atmosphere_radius;    // Outer edge of atmosphere
    float rayleigh_scale;       // Rayleigh scattering coefficient
    float mie_scale;            // Mie scattering coefficient
    float mie_g;                // Mie phase asymmetry (0.85 = strong forward scatter)
    float sun_intensity;        // Sun brightness multiplier
} AtmosphereConfig;

typedef struct Atmosphere {
    sg_pipeline pipeline;
    AtmosphereConfig config;
} Atmosphere;

// Create default config tuned for a planet of the given effective surface radius
AtmosphereConfig atmosphere_default_config(float effective_surface_radius);

// Initialize atmosphere pipeline (shader + additive blend)
void atmosphere_init(Atmosphere* atmos, AtmosphereConfig config);

// Render atmosphere as fullscreen pass (call between sky and terrain draws)
// Uses the shared fullscreen triangle vertex buffer from the sky pass
void atmosphere_render(const Atmosphere* atmos, const sg_bindings* fullscreen_bind,
                       HMM_Vec3 camera_pos, HMM_Vec3 sun_dir, HMM_Mat4 inv_vp);

void atmosphere_destroy(Atmosphere* atmos);

#endif

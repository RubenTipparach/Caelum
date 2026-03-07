#ifndef SOLAR_CONFIG_H
#define SOLAR_CONFIG_H

#include "HandmadeMath.h"

#define SOLAR_MAX_MOONS 10

typedef enum {
    MOON_SHAPE_ELLIPSOID = 0,
    MOON_SHAPE_CAPSULE   = 1,
} MoonShapeType;

typedef struct TenebrisConfig {
    float base_radius;         /* 796,000 m */
    float surface_gravity;     /* 10.0 m/s^2 */
    int   noise_seed;          /* 42 */
    int   sea_level;           /* 24 (layers) */
    float layer_thickness;     /* 1.0 m */
} TenebrisConfig;

typedef struct MoonConfig {
    const char* name;
    /* Shape */
    MoonShapeType shape_type;
    float base_radius;
    float ellipsoid_scale[3];
    float surface_gravity;
    /* Noise */
    int   noise_seed;
    float noise_frequency;
    float noise_amplitude;
    int   noise_octaves;
    /* Orbit */
    float semi_major_axis;     /* distance from Tenebris center (m) */
    float eccentricity;
    float inclination_deg;
    float lon_ascending_deg;
    float period;              /* orbital period (seconds) */
    /* Color */
    HMM_Vec3 base_color;
    HMM_Vec3 highlight_color;
    HMM_Vec3 shadow_color;
} MoonConfig;

typedef struct SolarSystemConfig {
    TenebrisConfig tenebris;
    MoonConfig moons[SOLAR_MAX_MOONS];
    int moon_count;
} SolarSystemConfig;

/* Returns the default config (all current hardcoded values) */
SolarSystemConfig solar_system_default_config(void);

#endif

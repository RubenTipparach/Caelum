#include "solar_config.h"

SolarSystemConfig solar_system_default_config(void) {
    SolarSystemConfig cfg = {0};

    /* ---- Tenebris (main planet) ---- */
    cfg.tenebris.base_radius     = 796000.0f;
    cfg.tenebris.surface_gravity = 20.0f;
    cfg.tenebris.noise_seed      = 42;
    cfg.tenebris.sea_level       = 24;
    cfg.tenebris.layer_thickness = 1.0f;
    cfg.tenebris.lod_fade_start  = 250.0f;
    cfg.tenebris.lod_fade_end    = 400.0f;

    /* ---- 10 Moons ---- */
    /*
     * All moons are grey rocky bodies — cool blue-gray tones to match moon.png texture.
     * Orbital periods: Kepler's 3rd law scaling (T ~ a^1.5).
     * noise_amp = 0.01 for all moons (1% of radius).
     * Capsule shape for axis ratio > 1.5: Dravok (#3), Cryx (#5), Vexis (#8).
     */
    cfg.moon_count = SOLAR_MAX_MOONS;

    /* --- 2 giant moons (50-80km radius) --- */
    cfg.moons[0] = (MoonConfig){
        .name = "Gorrath", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 80000.0f, .ellipsoid_scale = {1.0f, 0.9f, 1.05f},
        .surface_gravity = 12.0f,
        .noise_seed = 100, .noise_frequency = 0.3f, .noise_amplitude = 0.01f, .noise_octaves = 4,
        .semi_major_axis = 5000000.0f, .eccentricity = 0.02f,
        .inclination_deg = 5.0f, .lon_ascending_deg = 0.0f, .period = 14230.0f,
        .base_color = {{0.38f, 0.38f, 0.42f}},
        .highlight_color = {{0.52f, 0.52f, 0.57f}},
        .shadow_color = {{0.18f, 0.18f, 0.22f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[1] = (MoonConfig){
        .name = "Atheron", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 55000.0f, .ellipsoid_scale = {1.05f, 0.95f, 1.0f},
        .surface_gravity = 9.0f,
        .noise_seed = 150, .noise_frequency = 0.4f, .noise_amplitude = 0.01f, .noise_octaves = 4,
        .semi_major_axis = 10000000.0f, .eccentricity = 0.03f,
        .inclination_deg = 10.0f, .lon_ascending_deg = 90.0f, .period = 40250.0f,
        .base_color = {{0.40f, 0.40f, 0.44f}},
        .highlight_color = {{0.55f, 0.55f, 0.60f}},
        .shadow_color = {{0.20f, 0.20f, 0.24f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    /* --- 3 large moons (5-7km radius) --- */
    cfg.moons[2] = (MoonConfig){
        .name = "Kelthos", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 7000.0f, .ellipsoid_scale = {0.95f, 1.10f, 1.0f},
        .surface_gravity = 4.0f,
        .noise_seed = 300, .noise_frequency = 0.5f, .noise_amplitude = 0.01f, .noise_octaves = 4,
        .semi_major_axis = 3000000.0f, .eccentricity = 0.04f,
        .inclination_deg = 12.0f, .lon_ascending_deg = 45.0f, .period = 6600.0f,
        .base_color = {{0.36f, 0.37f, 0.40f}},
        .highlight_color = {{0.50f, 0.51f, 0.55f}},
        .shadow_color = {{0.17f, 0.17f, 0.20f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[3] = (MoonConfig){
        .name = "Dravok", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 6000.0f, .ellipsoid_scale = {1.20f, 0.80f, 0.95f},
        .surface_gravity = 3.6f,
        .noise_seed = 400, .noise_frequency = 0.18f, .noise_amplitude = 0.01f, .noise_octaves = 3,
        .semi_major_axis = 6500000.0f, .eccentricity = 0.06f,
        .inclination_deg = 18.0f, .lon_ascending_deg = 200.0f, .period = 21100.0f,
        .base_color = {{0.34f, 0.34f, 0.38f}},
        .highlight_color = {{0.48f, 0.48f, 0.53f}},
        .shadow_color = {{0.16f, 0.16f, 0.19f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[4] = (MoonConfig){
        .name = "Serath", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 5000.0f, .ellipsoid_scale = {1.0f, 1.0f, 1.15f},
        .surface_gravity = 3.0f,
        .noise_seed = 500, .noise_frequency = 0.4f, .noise_amplitude = 0.01f, .noise_octaves = 3,
        .semi_major_axis = 8500000.0f, .eccentricity = 0.03f,
        .inclination_deg = 8.0f, .lon_ascending_deg = 300.0f, .period = 31560.0f,
        .base_color = {{0.42f, 0.42f, 0.46f}},
        .highlight_color = {{0.57f, 0.57f, 0.62f}},
        .shadow_color = {{0.22f, 0.22f, 0.25f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    /* --- 5 small moons (1-4km radius) --- */
    cfg.moons[5] = (MoonConfig){
        .name = "Cryx", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 3000.0f, .ellipsoid_scale = {1.3f, 0.7f, 1.0f},
        .surface_gravity = 2.0f,
        .noise_seed = 600, .noise_frequency = 0.15f, .noise_amplitude = 0.01f, .noise_octaves = 2,
        .semi_major_axis = 2000000.0f, .eccentricity = 0.08f,
        .inclination_deg = 25.0f, .lon_ascending_deg = 60.0f, .period = 3600.0f,
        .base_color = {{0.32f, 0.32f, 0.36f}},
        .highlight_color = {{0.46f, 0.46f, 0.51f}},
        .shadow_color = {{0.15f, 0.15f, 0.18f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[6] = (MoonConfig){
        .name = "Nyctra", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 2000.0f, .ellipsoid_scale = {0.8f, 1.2f, 1.1f},
        .surface_gravity = 1.6f,
        .noise_seed = 700, .noise_frequency = 0.18f, .noise_amplitude = 0.01f, .noise_octaves = 2,
        .semi_major_axis = 4500000.0f, .eccentricity = 0.05f,
        .inclination_deg = 30.0f, .lon_ascending_deg = 150.0f, .period = 12150.0f,
        .base_color = {{0.37f, 0.38f, 0.42f}},
        .highlight_color = {{0.52f, 0.53f, 0.58f}},
        .shadow_color = {{0.18f, 0.18f, 0.22f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[7] = (MoonConfig){
        .name = "Thalwen", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 4000.0f, .ellipsoid_scale = {1.1f, 0.9f, 1.2f},
        .surface_gravity = 2.4f,
        .noise_seed = 800, .noise_frequency = 0.15f, .noise_amplitude = 0.01f, .noise_octaves = 3,
        .semi_major_axis = 7500000.0f, .eccentricity = 0.07f,
        .inclination_deg = 15.0f, .lon_ascending_deg = 270.0f, .period = 26200.0f,
        .base_color = {{0.35f, 0.35f, 0.40f}},
        .highlight_color = {{0.50f, 0.50f, 0.56f}},
        .shadow_color = {{0.17f, 0.17f, 0.21f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[8] = (MoonConfig){
        .name = "Vexis", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 1500.0f, .ellipsoid_scale = {1.4f, 0.6f, 1.0f},
        .surface_gravity = 1.2f,
        .noise_seed = 900, .noise_frequency = 0.21f, .noise_amplitude = 0.01f, .noise_octaves = 2,
        .semi_major_axis = 11000000.0f, .eccentricity = 0.09f,
        .inclination_deg = 40.0f, .lon_ascending_deg = 90.0f, .period = 46440.0f,
        .base_color = {{0.33f, 0.33f, 0.38f}},
        .highlight_color = {{0.47f, 0.47f, 0.53f}},
        .shadow_color = {{0.15f, 0.15f, 0.19f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    cfg.moons[9] = (MoonConfig){
        .name = "Zephyros", .shape_type = MOON_SHAPE_ELLIPSOID,
        .base_radius = 1000.0f, .ellipsoid_scale = {1.0f, 1.3f, 0.8f},
        .surface_gravity = 1.0f,
        .noise_seed = 1000, .noise_frequency = 0.21f, .noise_amplitude = 0.01f, .noise_octaves = 2,
        .semi_major_axis = 13000000.0f, .eccentricity = 0.02f,
        .inclination_deg = 10.0f, .lon_ascending_deg = 330.0f, .period = 59720.0f,
        .base_color = {{0.39f, 0.40f, 0.44f}},
        .highlight_color = {{0.54f, 0.55f, 0.60f}},
        .shadow_color = {{0.19f, 0.20f, 0.23f}},
        .lod_fade_start = 250.0f, .lod_fade_end = 400.0f,
    };

    return cfg;
}

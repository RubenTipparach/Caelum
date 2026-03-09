#ifndef TERRAIN_NOISE_H
#define TERRAIN_NOISE_H

// Shared terrain noise functions used by both hex_terrain.c and erosion.c.
// Header-only (static inline) to guarantee identical noise sampling.

#include "HandmadeMath.h"
#include "FastNoiseLite.h"
#include <math.h>

// ---- Terrain height constants ----
#define TERRAIN_SEA_LEVEL_M   4000.0f
#define TERRAIN_AMPLITUDE_M   8000.0f
#define TERRAIN_MIN_M         500.0f

static inline float ht_smoothstepf(float edge0, float edge1, float x) {
    float t = (x - edge0) / (edge1 - edge0);
    t = fminf(1.0f, fmaxf(0.0f, t));
    return t * t * (3.0f - 2.0f * t);
}

static inline fnl_state ht_create_continental_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 3;
    noise.frequency = 0.6f;
    noise.seed = seed;
    return noise;
}

static inline fnl_state ht_create_mountain_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_RIDGED;
    noise.octaves = 5;
    noise.frequency = 1.5f;
    noise.seed = seed + 4000;
    return noise;
}

static inline fnl_state ht_create_warp_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 3;
    noise.frequency = 4.0f;
    noise.seed = seed + 1000;
    return noise;
}

static inline fnl_state ht_create_detail_noise(int seed) {
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_RIDGED;
    noise.octaves = 3;
    noise.frequency = 16.0f;
    noise.seed = seed + 2000;
    return noise;
}

static inline float ht_sample_terrain_noise(fnl_state* continental, fnl_state* mountain,
                                             fnl_state* warp_noise, fnl_state* detail_noise,
                                             HMM_Vec3 unit_pos) {
    float scale = 3.0f;
    float px = unit_pos.X * scale;
    float py = unit_pos.Y * scale;
    float pz = unit_pos.Z * scale;

    float warp_strength = 0.5f;
    float wx = fnlGetNoise3D(warp_noise, px + 5.2f, py + 1.3f, pz + 3.7f);
    float wy = fnlGetNoise3D(warp_noise, px + 9.1f, py + 4.8f, pz + 7.2f);
    float wz = fnlGetNoise3D(warp_noise, px + 2.6f, py + 8.4f, pz + 0.9f);
    float wpx = px + wx * warp_strength;
    float wpy = py + wy * warp_strength;
    float wpz = pz + wz * warp_strength;

    float continent = fnlGetNoise3D(continental, wpx, wpy, wpz);

    float mountain_raw = fnlGetNoise3D(mountain, wpx, wpy, wpz);
    float mountain_val = (mountain_raw + 1.0f) * 0.5f;
    mountain_val *= mountain_val;
    float land_factor = ht_smoothstepf(-0.05f, 0.35f, continent);
    float mountain_height = mountain_val * land_factor;

    float detail = fnlGetNoise3D(detail_noise, px, py, pz);
    float detail_weight = 0.05f + land_factor * 0.10f;

    float height = continent * 0.55f - 0.22f;
    height += mountain_height * 0.45f;
    height += detail * detail_weight;

    // Power redistribution: steepen mountains, flatten plains/ocean floors.
    if (height > 0.0f) {
        height = powf(height, 1.35f);
    } else {
        height = -powf(-height, 0.8f);
    }

    if (height > 1.0f) height = 1.0f;
    if (height < -1.0f) height = -1.0f;
    return height;
}

static inline float ht_sample_height_m(fnl_state* continental, fnl_state* mountain,
                                        fnl_state* warp_noise, fnl_state* detail_noise,
                                        HMM_Vec3 unit_pos) {
    float n = ht_sample_terrain_noise(continental, mountain, warp_noise, detail_noise, unit_pos);
    float height = TERRAIN_MIN_M + (n + 1.0f) * 0.5f * TERRAIN_AMPLITUDE_M;
    if (height < 0.0f) height = 0.0f;
    return height;
}

#endif

#ifndef EROSION_H
#define EROSION_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "job_system.h"

// ---- Configuration ----
#define EROSION_REGION_SAMPLES    256       // Grid resolution per region (interior)
#define EROSION_SAMPLE_SPACING    3.0f      // Meters between samples
#define EROSION_REGION_SIZE_M     (EROSION_REGION_SAMPLES * EROSION_SAMPLE_SPACING)  // 768m
#define EROSION_OVERLAP           32        // Border overlap samples per side (96m)
#define EROSION_TOTAL_SAMPLES     (EROSION_REGION_SAMPLES + 2 * EROSION_OVERLAP)    // 320
#define EROSION_MAX_CACHED        32        // LRU cache slots
#define EROSION_DEFAULT_ITERS     5000      // Default droplets per region (tunable)

// ---- Erosion parameters (tunable) ----
typedef struct ErosionParams {
    int   iterations;        // Droplets per region (50000)
    int   max_lifetime;      // Max steps per droplet (80)
    int   brush_radius;      // Deposition/erosion brush (3)
    float inertia;           // 0.3  (0=instant turn, 1=no turn)
    float sediment_capacity; // 8.0  (max sediment factor)
    float erosion_rate;      // 0.7  (fraction of capacity deficit eroded)
    float deposition_rate;   // 0.02
    float evaporation_rate;  // 0.02
    float min_slope;         // 0.01 (minimum slope for capacity calc)
    float gravity;           // 10.0
} ErosionParams;

// ---- Region key (stable on unit sphere, tangent-frame-independent) ----
typedef struct ErosionRegionKey {
    int32_t lon_idx;
    int32_t lat_idx;
} ErosionRegionKey;

// ---- Region state ----
typedef enum {
    EROSION_STATE_EMPTY,       // Slot unused
    EROSION_STATE_COMPUTING,   // Job submitted, not yet done
    EROSION_STATE_READY,       // Delta map available for sampling
} ErosionRegionState;

typedef struct ErosionRegion {
    ErosionRegionKey key;
    ErosionRegionState state;

    // Region geometry (deterministic from key)
    HMM_Vec3 center_unit;       // Unit sphere center direction
    HMM_Vec3 east, north, up;   // Tangent frame at region center

    // Height delta map: EROSION_REGION_SAMPLES x EROSION_REGION_SAMPLES
    // delta = eroded_height - raw_noise_height (meters)
    float* delta_heights;

    // LRU tracking
    uint64_t last_access_frame;

    // Job completion (written by worker, read by main thread)
    volatile int job_completed;
    void* pending_job;          // ErosionJob* (opaque, freed after transfer)
} ErosionRegion;

// ---- Erosion cache (embedded in HexTerrain) ----
typedef struct ErosionCache {
    ErosionRegion regions[EROSION_MAX_CACHED];
    int seed;
    float planet_radius;
    ErosionParams params;
    uint32_t params_hash;       // FNV-1a of params for cache invalidation
    uint64_t current_frame;
    bool new_region_ready;      // Set when a region transitions to READY this frame
    char cache_dir[256];
} ErosionCache;

// ---- Public API ----

void erosion_init(ErosionCache* cache, int seed, float planet_radius, const char* cache_dir);
void erosion_destroy(ErosionCache* cache);

ErosionParams erosion_default_params(void);
void erosion_set_params(ErosionCache* cache, const ErosionParams* params);

// Called once per frame from hex_terrain_update. Prefetches nearby regions,
// polls computing jobs, saves completed results to disk.
void erosion_update(ErosionCache* cache, HMM_Vec3 camera_unit_pos, JobSystem* jobs);

// Thread-safe read-only lookup. Returns height delta in meters at unit_pos.
// Returns 0.0 if the region hasn't been computed yet.
float erosion_sample_delta(const ErosionCache* cache, HMM_Vec3 unit_pos);

// Precompute erosion for a grid of regions around a spawn point.
// Checks disk cache first, only computes missing regions.
// Blocks until all regions are ready. Call at startup before first frame.
void erosion_precompute(ErosionCache* cache, HMM_Vec3 center_unit,
                        JobSystem* jobs, int radius);

#endif

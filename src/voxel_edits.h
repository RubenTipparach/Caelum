#ifndef VOXEL_EDITS_H
#define VOXEL_EDITS_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "job_system.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

// ---- Configuration ----
#define EDIT_SECTOR_SIZE_M      256.0f      // Surface meters per sector side
#define EDIT_MAX_CACHED         32          // LRU cache slots
#define EDIT_MAX_PER_SECTOR     4096        // Max edits per sector file
#define EDIT_PRELOAD_RADIUS     2           // Preload sectors in NxN neighborhood

// ---- Sector key (cube-map face + grid index, no polar singularity) ----
typedef struct EditSectorKey {
    int8_t face;        // 0-5: +X, -X, +Y, -Y, +Z, -Z
    int8_t pad;
    int16_t u_idx;      // grid index on face U axis
    int16_t v_idx;      // grid index on face V axis
} EditSectorKey;

// ---- Single voxel edit ----
// Position stored as 64-bit unit-sphere direction (no polar singularity).
// Sector file provides the "chunk address"; this stores the precise hex position.
typedef struct VoxelEdit {
    double ux, uy, uz;       // Unit-sphere direction (64-bit, stable everywhere)
    int32_t layer;           // World layer index (absolute)
    uint8_t voxel_type;      // VoxelType to set (VOXEL_AIR = break)
    uint8_t pad[3];          // Alignment (32 bytes total)
} VoxelEdit;

// ---- Sector state ----
typedef enum {
    EDIT_STATE_EMPTY,
    EDIT_STATE_LOADING,
    EDIT_STATE_READY,
} EditSectorState;

typedef struct EditSector {
    EditSectorKey key;
    EditSectorState state;

    VoxelEdit* edits;
    int edit_count;
    int edit_capacity;

    bool dirty;
    uint64_t last_access_frame;

    volatile int load_completed;
    VoxelEdit* loaded_edits;
    int loaded_count;
} EditSector;

// ---- Edit cache ----
typedef struct EditCache {
    EditSector sectors[EDIT_MAX_CACHED];
    int seed;
    float planet_radius;
    uint64_t current_frame;

    char cache_dir[256];
    float flush_timer;

#ifdef _WIN32
    CRITICAL_SECTION mutex;
#else
    pthread_mutex_t mutex;
#endif
} EditCache;

// ---- Public API ----

void edit_cache_init(EditCache* cache, int seed, float planet_radius, const char* cache_dir);
void edit_cache_destroy(EditCache* cache);

// Called once per frame. Preloads nearby sectors, flushes dirty periodically.
void edit_cache_update(EditCache* cache, HMM_Vec3 camera_unit_pos);

// Record a voxel edit. Called from hex_terrain_break/place on main thread.
void edit_cache_record(EditCache* cache, double ux, double uy, double uz,
                       int32_t layer, uint8_t voxel_type);

// Query edits overlapping a cube-map sector key range.
// Loads sectors from disk as needed. Returns number of edits written.
int edit_cache_query_area(const EditCache* cache,
                          int8_t face,
                          int16_t u_idx_min, int16_t u_idx_max,
                          int16_t v_idx_min, int16_t v_idx_max,
                          VoxelEdit* out_edits, int max_edits);

// Force-flush all dirty sectors to disk synchronously (called on shutdown).
void edit_cache_flush_all(EditCache* cache);

// ---- Cube-map helpers (shared with hex_terrain.c) ----

// Project a unit vector to cube-map (face, u, v). u,v in [-1, 1].
void edit_unit_to_cubemap(HMM_Vec3 unit, int8_t* face, float* u, float* v);

// Compute sector key from a unit vector.
EditSectorKey edit_sector_key(HMM_Vec3 unit, float planet_radius);

#endif

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
#define EDIT_FLUSH_INTERVAL_S   5.0f        // Auto-flush dirty sectors every N seconds
#define EDIT_PRELOAD_RADIUS     2           // Preload sectors in NxN neighborhood

// ---- Sector key (stable on unit sphere) ----
typedef struct EditSectorKey {
    int32_t lon_idx;
    int32_t lat_idx;
} EditSectorKey;

// ---- Single voxel edit ----
typedef struct VoxelEdit {
    float lon;              // Spherical longitude (radians)
    float lat;              // Spherical latitude (radians)
    int32_t layer;          // World layer index (absolute)
    uint8_t voxel_type;     // VoxelType to set (VOXEL_AIR = break)
    uint8_t pad[3];         // Alignment padding (20 bytes total)
} VoxelEdit;

// ---- Sector state ----
typedef enum {
    EDIT_STATE_EMPTY,       // Slot unused
    EDIT_STATE_LOADING,     // Being loaded from disk on worker thread
    EDIT_STATE_READY,       // Available for read/write
} EditSectorState;

typedef struct EditSector {
    EditSectorKey key;
    EditSectorState state;

    VoxelEdit* edits;       // Dynamic array of edits
    int edit_count;
    int edit_capacity;

    bool dirty;             // Has unsaved changes
    uint64_t last_access_frame;

    // Loading completion (written by worker, read by main thread)
    volatile int load_completed;
    VoxelEdit* loaded_edits;    // Loaded from disk (transferred on completion)
    int loaded_count;
} EditSector;

// ---- Edit cache ----
typedef struct EditCache {
    EditSector sectors[EDIT_MAX_CACHED];
    int seed;
    float planet_radius;
    uint64_t current_frame;

    char cache_dir[256];

    // Background flush state
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
void edit_cache_record(EditCache* cache, float lon, float lat,
                       int32_t layer, uint8_t voxel_type);

// Thread-safe query for edits in a bounding box. Called from worker threads.
// Returns number of edits written to out_edits (up to max_edits).
int edit_cache_query_area(const EditCache* cache,
                          float lon_min, float lon_max,
                          float lat_min, float lat_max,
                          VoxelEdit* out_edits, int max_edits);

// Force-flush all dirty sectors to disk synchronously (called on shutdown).
void edit_cache_flush_all(EditCache* cache);

#endif

#include "voxel_edits.h"
#include "math_utils.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#define mkdir_p(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <unistd.h>
#define mkdir_p(path) mkdir(path, 0755)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---- File format ----

#define EDIT_FILE_MAGIC   0x54494445  // "EDIT"
#define EDIT_FILE_VERSION 1

typedef struct EditFileHeader {
    uint32_t magic;
    uint32_t version;
    int32_t  seed;
    int32_t  lon_idx;
    int32_t  lat_idx;
    int32_t  edit_count;
    int32_t  reserved[2];   // 32 bytes total
} EditFileHeader;

// ---- Sector key computation ----

static float edit_sector_angular_size(float planet_radius) {
    return EDIT_SECTOR_SIZE_M / planet_radius;
}

static EditSectorKey edit_sector_key_from_lonlat(float lon, float lat, float planet_radius) {
    float sector_rad = edit_sector_angular_size(planet_radius);
    EditSectorKey key;
    key.lon_idx = (int32_t)floorf(lon / sector_rad);
    key.lat_idx = (int32_t)floorf(lat / sector_rad);
    return key;
}

static EditSectorKey edit_sector_key_from_unit(HMM_Vec3 unit_pos, float planet_radius) {
    float lat = asinf(unit_pos.Y);
    float lon = atan2f(unit_pos.Z, unit_pos.X);
    return edit_sector_key_from_lonlat(lon, lat, planet_radius);
}

// ---- Directory helpers ----

static void edit_ensure_dir(const char* path) {
    char buf[512];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char* p = buf + 1; *p; p++) {
        if (*p == '/' || *p == '\\') {
            *p = '\0';
            mkdir_p(buf);
            *p = '/';
        }
    }
    mkdir_p(buf);
}

static void edit_sector_path(const EditCache* cache, const EditSectorKey* key, char* out, int out_size) {
    snprintf(out, out_size, "%s/%d/%d_%d.bin",
             cache->cache_dir, cache->seed, key->lon_idx, key->lat_idx);
}

// ---- Disk I/O ----

static bool edit_save_sector(const EditCache* cache, const EditSector* sector) {
    char path[512];
    edit_sector_path(cache, &sector->key, path, sizeof(path));

    // Ensure directory exists
    char dir[512];
    snprintf(dir, sizeof(dir), "%s/%d", cache->cache_dir, cache->seed);
    edit_ensure_dir(dir);

    FILE* f = fopen(path, "wb");
    if (!f) return false;

    EditFileHeader hdr = {0};
    hdr.magic = EDIT_FILE_MAGIC;
    hdr.version = EDIT_FILE_VERSION;
    hdr.seed = cache->seed;
    hdr.lon_idx = sector->key.lon_idx;
    hdr.lat_idx = sector->key.lat_idx;
    hdr.edit_count = sector->edit_count;

    fwrite(&hdr, sizeof(hdr), 1, f);
    if (sector->edit_count > 0 && sector->edits) {
        fwrite(sector->edits, sizeof(VoxelEdit), sector->edit_count, f);
    }
    fclose(f);
    return true;
}

static bool edit_load_sector_from_disk(const EditCache* cache, const EditSectorKey* key,
                                        VoxelEdit** out_edits, int* out_count) {
    char path[512];
    edit_sector_path(cache, key, path, sizeof(path));

    FILE* f = fopen(path, "rb");
    if (!f) {
        *out_edits = NULL;
        *out_count = 0;
        return false;
    }

    EditFileHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1 ||
        hdr.magic != EDIT_FILE_MAGIC ||
        hdr.version != EDIT_FILE_VERSION ||
        hdr.seed != cache->seed) {
        fclose(f);
        *out_edits = NULL;
        *out_count = 0;
        return false;
    }

    int count = hdr.edit_count;
    if (count < 0 || count > EDIT_MAX_PER_SECTOR) {
        fclose(f);
        *out_edits = NULL;
        *out_count = 0;
        return false;
    }

    VoxelEdit* edits = NULL;
    if (count > 0) {
        edits = (VoxelEdit*)malloc((size_t)count * sizeof(VoxelEdit));
        if (!edits || fread(edits, sizeof(VoxelEdit), count, f) != (size_t)count) {
            free(edits);
            fclose(f);
            *out_edits = NULL;
            *out_count = 0;
            return false;
        }
    }

    fclose(f);
    *out_edits = edits;
    *out_count = count;
    return true;
}

// ---- Forward declarations ----
static void edit_load_into_slot(EditCache* cache, int slot, const EditSectorKey* key);

// ---- LRU cache management ----

static int edit_find_slot(const EditCache* cache, const EditSectorKey* key) {
    for (int i = 0; i < EDIT_MAX_CACHED; i++) {
        const EditSector* s = &cache->sectors[i];
        if (s->state != EDIT_STATE_EMPTY &&
            s->key.lon_idx == key->lon_idx &&
            s->key.lat_idx == key->lat_idx) {
            return i;
        }
    }
    return -1;
}

static int edit_find_or_evict_slot(EditCache* cache) {
    // Find empty slot
    for (int i = 0; i < EDIT_MAX_CACHED; i++) {
        if (cache->sectors[i].state == EDIT_STATE_EMPTY)
            return i;
    }

    // Evict least-recently-used READY slot
    int best = -1;
    uint64_t oldest = UINT64_MAX;
    for (int i = 0; i < EDIT_MAX_CACHED; i++) {
        if (cache->sectors[i].state == EDIT_STATE_READY &&
            cache->sectors[i].last_access_frame < oldest) {
            oldest = cache->sectors[i].last_access_frame;
            best = i;
        }
    }

    if (best >= 0) {
        EditSector* s = &cache->sectors[best];
        // Save dirty sector before evicting
        if (s->dirty) {
            edit_save_sector(cache, s);
            s->dirty = false;
        }
        free(s->edits);
        s->edits = NULL;
        s->edit_count = 0;
        s->edit_capacity = 0;
        s->state = EDIT_STATE_EMPTY;
    }

    return best;
}

static void edit_sector_add(EditSector* sector, const VoxelEdit* edit) {
    if (sector->edit_count >= EDIT_MAX_PER_SECTOR) return;

    if (sector->edit_count >= sector->edit_capacity) {
        int new_cap = sector->edit_capacity < 16 ? 16 : sector->edit_capacity * 2;
        if (new_cap > EDIT_MAX_PER_SECTOR) new_cap = EDIT_MAX_PER_SECTOR;
        VoxelEdit* new_edits = (VoxelEdit*)realloc(sector->edits, (size_t)new_cap * sizeof(VoxelEdit));
        if (!new_edits) return;
        sector->edits = new_edits;
        sector->edit_capacity = new_cap;
    }

    sector->edits[sector->edit_count++] = *edit;
    sector->dirty = true;
}

// ---- Public API ----

void edit_cache_init(EditCache* cache, int seed, float planet_radius, const char* cache_dir) {
    memset(cache, 0, sizeof(EditCache));
    cache->seed = seed;
    cache->planet_radius = planet_radius;
    snprintf(cache->cache_dir, sizeof(cache->cache_dir), "%s", cache_dir);

    edit_ensure_dir(cache_dir);

#ifdef _WIN32
    InitializeCriticalSection(&cache->mutex);
#else
    pthread_mutex_init(&cache->mutex, NULL);
#endif

    printf("[EDITS] Initialized: seed=%d, sector=%.0fm\n", seed, EDIT_SECTOR_SIZE_M);
    fflush(stdout);
}

void edit_cache_destroy(EditCache* cache) {
    for (int i = 0; i < EDIT_MAX_CACHED; i++) {
        EditSector* s = &cache->sectors[i];
        free(s->edits);
        s->edits = NULL;
        s->state = EDIT_STATE_EMPTY;
    }

#ifdef _WIN32
    DeleteCriticalSection(&cache->mutex);
#else
    pthread_mutex_destroy(&cache->mutex);
#endif
}

void edit_cache_record(EditCache* cache, float lon, float lat,
                       int32_t layer, uint8_t voxel_type) {
    EditSectorKey key = edit_sector_key_from_lonlat(lon, lat, cache->planet_radius);

#ifdef _WIN32
    EnterCriticalSection(&cache->mutex);
#else
    pthread_mutex_lock(&cache->mutex);
#endif

    int slot = edit_find_slot(cache, &key);
    if (slot < 0) {
        // Sector not cached — find or evict a slot, load synchronously
        slot = edit_find_or_evict_slot(cache);
        if (slot < 0) {
#ifdef _WIN32
            LeaveCriticalSection(&cache->mutex);
#else
            pthread_mutex_unlock(&cache->mutex);
#endif
            return;
        }

        EditSector* s = &cache->sectors[slot];
        s->key = key;
        s->state = EDIT_STATE_READY;
        s->last_access_frame = cache->current_frame;
        s->dirty = false;
        s->edit_count = 0;
        s->edit_capacity = 0;
        s->edits = NULL;

        // Try loading from disk (synchronous — edits are small files)
        VoxelEdit* loaded = NULL;
        int loaded_count = 0;
        edit_load_sector_from_disk(cache, &key, &loaded, &loaded_count);
        if (loaded && loaded_count > 0) {
            s->edits = loaded;
            s->edit_count = loaded_count;
            s->edit_capacity = loaded_count;
        } else {
            free(loaded);
        }
    }

    EditSector* s = &cache->sectors[slot];
    s->last_access_frame = cache->current_frame;

    VoxelEdit edit;
    edit.lon = lon;
    edit.lat = lat;
    edit.layer = layer;
    edit.voxel_type = voxel_type;
    memset(edit.pad, 0, sizeof(edit.pad));

    edit_sector_add(s, &edit);

#ifdef _WIN32
    LeaveCriticalSection(&cache->mutex);
#else
    pthread_mutex_unlock(&cache->mutex);
#endif
}

int edit_cache_query_area(const EditCache* cache,
                          float lon_min, float lon_max,
                          float lat_min, float lat_max,
                          VoxelEdit* out_edits, int max_edits) {
    int count = 0;

    // Cast away const for mutex — query is logically const
    EditCache* mutable_cache = (EditCache*)cache;

#ifdef _WIN32
    EnterCriticalSection(&mutable_cache->mutex);
#else
    pthread_mutex_lock(&mutable_cache->mutex);
#endif

    float sector_rad = edit_sector_angular_size(cache->planet_radius);

    // Compute which sector keys overlap the query bounding box
    int32_t key_lon_min = (int32_t)floorf(lon_min / sector_rad);
    int32_t key_lon_max = (int32_t)floorf(lon_max / sector_rad);
    int32_t key_lat_min = (int32_t)floorf(lat_min / sector_rad);
    int32_t key_lat_max = (int32_t)floorf(lat_max / sector_rad);

    // For each sector key in the query box, ensure it's loaded
    for (int32_t ky = key_lat_min; ky <= key_lat_max && count < max_edits; ky++) {
        for (int32_t kx = key_lon_min; kx <= key_lon_max && count < max_edits; kx++) {
            EditSectorKey key;
            key.lon_idx = kx;
            key.lat_idx = ky;

            // Check if already in cache
            int slot = edit_find_slot(mutable_cache, &key);
            if (slot < 0) {
                // Not cached — try loading from disk on the fly
                // (edit files are tiny, <10KB, safe to do synchronously)
                int free_slot = edit_find_or_evict_slot(mutable_cache);
                if (free_slot >= 0) {
                    edit_load_into_slot(mutable_cache, free_slot, &key);
                    slot = free_slot;
                }
            }

            if (slot < 0) continue;
            const EditSector* s = &cache->sectors[slot];
            if (s->state != EDIT_STATE_READY || s->edit_count == 0) continue;

            // Check individual edits in this sector
            for (int e = 0; e < s->edit_count && count < max_edits; e++) {
                const VoxelEdit* ed = &s->edits[e];
                if (ed->lon >= lon_min && ed->lon <= lon_max &&
                    ed->lat >= lat_min && ed->lat <= lat_max) {
                    out_edits[count++] = *ed;
                }
            }
        }
    }

#ifdef _WIN32
    LeaveCriticalSection(&mutable_cache->mutex);
#else
    pthread_mutex_unlock(&mutable_cache->mutex);
#endif

    return count;
}

// Load a sector synchronously into a cache slot
static void edit_load_into_slot(EditCache* cache, int slot, const EditSectorKey* key) {
    EditSector* s = &cache->sectors[slot];
    memset(s, 0, sizeof(EditSector));
    s->key = *key;
    s->state = EDIT_STATE_READY;
    s->last_access_frame = cache->current_frame;

    VoxelEdit* loaded = NULL;
    int loaded_count = 0;
    if (edit_load_sector_from_disk(cache, key, &loaded, &loaded_count)) {
        s->edits = loaded;
        s->edit_count = loaded_count;
        s->edit_capacity = loaded_count;
    } else {
        free(loaded);
    }
}

void edit_cache_update(EditCache* cache, HMM_Vec3 camera_unit_pos) {
    cache->current_frame++;

#ifdef _WIN32
    EnterCriticalSection(&cache->mutex);
#else
    pthread_mutex_lock(&cache->mutex);
#endif

    // Preload nearby sectors (synchronous — files are < 10KB)
    {
        EditSectorKey center = edit_sector_key_from_unit(camera_unit_pos, cache->planet_radius);
        int loads = 0;

        for (int dy = -EDIT_PRELOAD_RADIUS; dy <= EDIT_PRELOAD_RADIUS && loads < 4; dy++) {
            for (int dx = -EDIT_PRELOAD_RADIUS; dx <= EDIT_PRELOAD_RADIUS && loads < 4; dx++) {
                EditSectorKey key;
                key.lon_idx = center.lon_idx + dx;
                key.lat_idx = center.lat_idx + dy;

                int slot = edit_find_slot(cache, &key);
                if (slot >= 0) {
                    cache->sectors[slot].last_access_frame = cache->current_frame;
                    continue;
                }

                int free_slot = edit_find_or_evict_slot(cache);
                if (free_slot < 0) continue;

                edit_load_into_slot(cache, free_slot, &key);
                loads++;
            }
        }
    }

    // Flush dirty sectors every ~300 frames (~5 seconds at 60fps)
    cache->flush_timer += 1.0f;
    if (cache->flush_timer >= 300.0f) {
        cache->flush_timer = 0.0f;

        for (int i = 0; i < EDIT_MAX_CACHED; i++) {
            EditSector* s = &cache->sectors[i];
            if (s->state == EDIT_STATE_READY && s->dirty) {
                edit_save_sector(cache, s);
                s->dirty = false;
            }
        }
    }

#ifdef _WIN32
    LeaveCriticalSection(&cache->mutex);
#else
    pthread_mutex_unlock(&cache->mutex);
#endif
}

void edit_cache_flush_all(EditCache* cache) {
#ifdef _WIN32
    EnterCriticalSection(&cache->mutex);
#else
    pthread_mutex_lock(&cache->mutex);
#endif

    for (int i = 0; i < EDIT_MAX_CACHED; i++) {
        EditSector* s = &cache->sectors[i];
        if (s->state == EDIT_STATE_READY && s->dirty) {
            edit_save_sector(cache, s);
            s->dirty = false;
        }
    }

#ifdef _WIN32
    LeaveCriticalSection(&cache->mutex);
#else
    pthread_mutex_unlock(&cache->mutex);
#endif

    printf("[EDITS] Flushed all dirty sectors to disk\n");
    fflush(stdout);
}

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

// ---- File format (v2: cube-map sectors, unit-vector edits) ----

#define EDIT_FILE_MAGIC   0x54494445  // "EDIT"
#define EDIT_FILE_VERSION 2

typedef struct EditFileHeader {
    uint32_t magic;         // 4
    uint32_t version;       // 4
    int32_t  seed;          // 4
    int32_t  edit_count;    // 4
    int8_t   face;          // 1
    int8_t   pad_h[3];      // 3
    int16_t  u_idx;         // 2
    int16_t  v_idx;         // 2
    int32_t  reserved[2];   // 8
} EditFileHeader;           // 32 bytes total

// ---- Cube-map projection ----

void edit_unit_to_cubemap(HMM_Vec3 unit, int8_t* face, float* u, float* v) {
    float ax = fabsf(unit.X), ay = fabsf(unit.Y), az = fabsf(unit.Z);
    if (ay >= ax && ay >= az) {
        *face = (unit.Y > 0.0f) ? 2 : 3;
        float inv = 1.0f / ay;
        *u = unit.X * inv;
        *v = unit.Z * inv;
    } else if (ax >= az) {
        *face = (unit.X > 0.0f) ? 0 : 1;
        float inv = 1.0f / ax;
        *u = unit.Y * inv;
        *v = unit.Z * inv;
    } else {
        *face = (unit.Z > 0.0f) ? 4 : 5;
        float inv = 1.0f / az;
        *u = unit.X * inv;
        *v = unit.Y * inv;
    }
}

static float edit_sector_uv_size(float planet_radius) {
    return EDIT_SECTOR_SIZE_M / planet_radius;
}

EditSectorKey edit_sector_key(HMM_Vec3 unit, float planet_radius) {
    int8_t face;
    float u, v;
    edit_unit_to_cubemap(unit, &face, &u, &v);

    float sector_uv = edit_sector_uv_size(planet_radius);
    EditSectorKey key;
    key.face = face;
    key.pad = 0;
    key.u_idx = (int16_t)floorf(u / sector_uv);
    key.v_idx = (int16_t)floorf(v / sector_uv);
    return key;
}

static bool edit_key_eq(const EditSectorKey* a, const EditSectorKey* b) {
    return a->face == b->face && a->u_idx == b->u_idx && a->v_idx == b->v_idx;
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

static void edit_sector_path(const EditCache* cache, const EditSectorKey* key,
                              char* out, int out_size) {
    snprintf(out, out_size, "%s/%d/f%d_%d_%d.bin",
             cache->cache_dir, cache->seed,
             (int)key->face, (int)key->u_idx, (int)key->v_idx);
}

// ---- Disk I/O ----

static bool edit_save_sector(const EditCache* cache, const EditSector* sector) {
    char path[512];
    edit_sector_path(cache, &sector->key, path, sizeof(path));

    char dir[512];
    snprintf(dir, sizeof(dir), "%s/%d", cache->cache_dir, cache->seed);
    edit_ensure_dir(dir);

    FILE* f = fopen(path, "wb");
    if (!f) return false;

    EditFileHeader hdr = {0};
    hdr.magic = EDIT_FILE_MAGIC;
    hdr.version = EDIT_FILE_VERSION;
    hdr.seed = cache->seed;
    hdr.edit_count = sector->edit_count;
    hdr.face = sector->key.face;
    hdr.u_idx = sector->key.u_idx;
    hdr.v_idx = sector->key.v_idx;

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
        if (s->state != EDIT_STATE_EMPTY && edit_key_eq(&s->key, key)) {
            return i;
        }
    }
    return -1;
}

static int edit_find_or_evict_slot(EditCache* cache) {
    for (int i = 0; i < EDIT_MAX_CACHED; i++) {
        if (cache->sectors[i].state == EDIT_STATE_EMPTY)
            return i;
    }

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

    printf("[EDITS] Initialized: seed=%d, sector=%.0fm, cube-map sectors\n",
           seed, EDIT_SECTOR_SIZE_M);
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

void edit_cache_record(EditCache* cache, double ux, double uy, double uz,
                       int32_t layer, uint8_t voxel_type) {
    HMM_Vec3 unit = {{(float)ux, (float)uy, (float)uz}};
    EditSectorKey key = edit_sector_key(unit, cache->planet_radius);

#ifdef _WIN32
    EnterCriticalSection(&cache->mutex);
#else
    pthread_mutex_lock(&cache->mutex);
#endif

    int slot = edit_find_slot(cache, &key);
    if (slot < 0) {
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
    edit.ux = ux;
    edit.uy = uy;
    edit.uz = uz;
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
                          int8_t face,
                          int16_t u_idx_min, int16_t u_idx_max,
                          int16_t v_idx_min, int16_t v_idx_max,
                          VoxelEdit* out_edits, int max_edits) {
    int count = 0;

    EditCache* mutable_cache = (EditCache*)cache;

#ifdef _WIN32
    EnterCriticalSection(&mutable_cache->mutex);
#else
    pthread_mutex_lock(&mutable_cache->mutex);
#endif

    for (int16_t vy = v_idx_min; vy <= v_idx_max && count < max_edits; vy++) {
        for (int16_t ux = u_idx_min; ux <= u_idx_max && count < max_edits; ux++) {
            EditSectorKey key;
            key.face = face;
            key.pad = 0;
            key.u_idx = ux;
            key.v_idx = vy;

            int slot = edit_find_slot(mutable_cache, &key);
            if (slot < 0) {
                int free_slot = edit_find_or_evict_slot(mutable_cache);
                if (free_slot >= 0) {
                    edit_load_into_slot(mutable_cache, free_slot, &key);
                    slot = free_slot;
                }
            }

            if (slot < 0) continue;
            const EditSector* s = &cache->sectors[slot];
            if (s->state != EDIT_STATE_READY || s->edit_count == 0) continue;

            // Copy all edits from this sector (caller does fine filtering)
            for (int e = 0; e < s->edit_count && count < max_edits; e++) {
                out_edits[count++] = s->edits[e];
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

void edit_cache_update(EditCache* cache, HMM_Vec3 camera_unit_pos) {
    cache->current_frame++;

#ifdef _WIN32
    EnterCriticalSection(&cache->mutex);
#else
    pthread_mutex_lock(&cache->mutex);
#endif

    // Preload nearby sectors
    {
        EditSectorKey center = edit_sector_key(camera_unit_pos, cache->planet_radius);
        int loads = 0;

        for (int dy = -EDIT_PRELOAD_RADIUS; dy <= EDIT_PRELOAD_RADIUS && loads < 4; dy++) {
            for (int dx = -EDIT_PRELOAD_RADIUS; dx <= EDIT_PRELOAD_RADIUS && loads < 4; dx++) {
                EditSectorKey key;
                key.face = center.face;
                key.pad = 0;
                key.u_idx = center.u_idx + (int16_t)dx;
                key.v_idx = center.v_idx + (int16_t)dy;

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

    // Flush dirty sectors periodically (~5 seconds at 60fps)
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

#include "erosion.h"
#include "terrain_noise.h"
#include "math_utils.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
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

// ---- Deterministic PRNG (splitmix64) ----

typedef struct {
    uint64_t state;
} ErosionRng;

static uint64_t erosion_rng_next(ErosionRng* rng) {
    rng->state += 0x9E3779B97F4A7C15ULL;
    uint64_t z = rng->state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static float erosion_rng_float(ErosionRng* rng) {
    return (float)(erosion_rng_next(rng) >> 40) / (float)(1ULL << 24);
}

static ErosionRng erosion_rng_create(int seed, int32_t lon_idx, int32_t lat_idx) {
    ErosionRng rng;
    rng.state = (uint64_t)(uint32_t)seed * 6364136223846793005ULL
              + (uint64_t)(uint32_t)lon_idx * 1442695040888963407ULL
              + (uint64_t)(uint32_t)lat_idx * 2862933555777941757ULL;
    erosion_rng_next(&rng);
    erosion_rng_next(&rng);
    return rng;
}

// ---- FNV-1a hash for params ----

static uint32_t fnv1a_hash(const void* data, size_t len) {
    const uint8_t* bytes = (const uint8_t*)data;
    uint32_t hash = 2166136261u;
    for (size_t i = 0; i < len; i++) {
        hash ^= bytes[i];
        hash *= 16777619u;
    }
    return hash;
}

// ---- Region key computation ----

static float erosion_region_angular_size(float planet_radius) {
    return EROSION_REGION_SIZE_M / planet_radius;
}

static ErosionRegionKey erosion_region_key(HMM_Vec3 unit_pos, float planet_radius) {
    float region_rad = erosion_region_angular_size(planet_radius);
    float lat = asinf(unit_pos.Y);
    float lon = atan2f(unit_pos.Z, unit_pos.X);
    ErosionRegionKey key;
    key.lon_idx = (int32_t)floorf(lon / region_rad);
    key.lat_idx = (int32_t)floorf(lat / region_rad);
    return key;
}

// ---- Region tangent frame (deterministic from key) ----

static HMM_Vec3 erosion_key_to_center(const ErosionRegionKey* key, float planet_radius) {
    float region_rad = erosion_region_angular_size(planet_radius);
    float lat = (key->lat_idx + 0.5f) * region_rad;
    float lon = (key->lon_idx + 0.5f) * region_rad;
    float cos_lat = cosf(lat);
    HMM_Vec3 center;
    center.X = cos_lat * cosf(lon);
    center.Y = sinf(lat);
    center.Z = cos_lat * sinf(lon);
    return center;
}

static void erosion_compute_tangent_frame(ErosionRegion* r) {
    r->up = r->center_unit;
    HMM_Vec3 ref = (HMM_Vec3){{0.0f, 1.0f, 0.0f}};
    if (fabsf(r->center_unit.Y) > 0.999f)
        ref = (HMM_Vec3){{1.0f, 0.0f, 0.0f}};
    r->east = vec3_normalize(vec3_cross(ref, r->center_unit));
    r->north = vec3_cross(r->center_unit, r->east);
}

// ---- Disk cache I/O ----

#define EROSION_FILE_MAGIC 0x584F5245  // "EROX"
#define EROSION_FILE_VERSION 1

typedef struct {
    uint32_t magic;
    uint32_t version;
    int32_t  seed;
    uint32_t params_hash;
    int32_t  lon_idx;
    int32_t  lat_idx;
    int32_t  grid_size;
    int32_t  reserved;
} ErosionFileHeader;

static void erosion_ensure_dir(const char* dir) {
    char buf[512];
    size_t len = strlen(dir);
    if (len >= sizeof(buf)) return;
    memcpy(buf, dir, len + 1);
    for (size_t i = 1; i < len; i++) {
        if (buf[i] == '/' || buf[i] == '\\') {
            buf[i] = '\0';
            mkdir_p(buf);
            buf[i] = '/';
        }
    }
    mkdir_p(buf);
}

static void erosion_file_path(const ErosionCache* cache, const ErosionRegionKey* key,
                               char* out_path, int path_len) {
    snprintf(out_path, path_len, "%s/%d_%08x/%d_%d.bin",
             cache->cache_dir, cache->seed, cache->params_hash,
             key->lon_idx, key->lat_idx);
}

static bool erosion_save_region(const ErosionCache* cache, const ErosionRegion* r) {
    char path[512];
    erosion_file_path(cache, &r->key, path, sizeof(path));

    // Ensure parent directory exists
    char dir[512];
    snprintf(dir, sizeof(dir), "%s/%d_%08x", cache->cache_dir, cache->seed, cache->params_hash);
    erosion_ensure_dir(dir);

    FILE* f = fopen(path, "wb");
    if (!f) return false;

    ErosionFileHeader hdr = {0};
    hdr.magic = EROSION_FILE_MAGIC;
    hdr.version = EROSION_FILE_VERSION;
    hdr.seed = cache->seed;
    hdr.params_hash = cache->params_hash;
    hdr.lon_idx = r->key.lon_idx;
    hdr.lat_idx = r->key.lat_idx;
    hdr.grid_size = EROSION_REGION_SAMPLES;

    size_t written = fwrite(&hdr, sizeof(hdr), 1, f);
    if (written != 1) { fclose(f); return false; }

    size_t data_size = (size_t)EROSION_REGION_SAMPLES * EROSION_REGION_SAMPLES * sizeof(float);
    written = fwrite(r->delta_heights, 1, data_size, f);
    fclose(f);
    return written == data_size;
}

static bool erosion_load_region(const ErosionCache* cache, ErosionRegion* r) {
    char path[512];
    erosion_file_path(cache, &r->key, path, sizeof(path));

    FILE* f = fopen(path, "rb");
    if (!f) return false;

    ErosionFileHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) { fclose(f); return false; }

    if (hdr.magic != EROSION_FILE_MAGIC ||
        hdr.version != EROSION_FILE_VERSION ||
        hdr.seed != cache->seed ||
        hdr.params_hash != cache->params_hash ||
        hdr.grid_size != EROSION_REGION_SAMPLES) {
        fclose(f);
        return false;
    }

    size_t data_size = (size_t)EROSION_REGION_SAMPLES * EROSION_REGION_SAMPLES * sizeof(float);
    if (!r->delta_heights) {
        r->delta_heights = (float*)malloc(data_size);
        if (!r->delta_heights) { fclose(f); return false; }
    }

    size_t read = fread(r->delta_heights, 1, data_size, f);
    fclose(f);
    return read == data_size;
}

// ---- Erosion job (runs on worker thread) ----

typedef struct ErosionJob {
    ErosionRegionKey key;
    HMM_Vec3 center_unit;
    HMM_Vec3 east, north, up;
    float planet_radius;
    int seed;
    ErosionParams params;

    // Output (transferred to ErosionRegion on completion)
    float* delta_out;

    // Completion flag
    volatile int completed;
    int region_slot;
} ErosionJob;

// ---- Brush application ----

static void apply_brush(float* heightmap, int grid_size, float px, float py,
                         float amount, int brush_radius) {
    int cx = (int)px, cy = (int)py;
    float total_weight = 0.0f;

    for (int dy = -brush_radius; dy <= brush_radius; dy++) {
        for (int dx = -brush_radius; dx <= brush_radius; dx++) {
            int gx = cx + dx, gy = cy + dy;
            if (gx < 0 || gx >= grid_size || gy < 0 || gy >= grid_size) continue;
            float dist = sqrtf((float)(dx * dx + dy * dy));
            if (dist > (float)brush_radius) continue;
            total_weight += 1.0f - dist / (float)brush_radius;
        }
    }
    if (total_weight < 0.0001f) return;

    for (int dy = -brush_radius; dy <= brush_radius; dy++) {
        for (int dx = -brush_radius; dx <= brush_radius; dx++) {
            int gx = cx + dx, gy = cy + dy;
            if (gx < 0 || gx >= grid_size || gy < 0 || gy >= grid_size) continue;
            float dist = sqrtf((float)(dx * dx + dy * dy));
            if (dist > (float)brush_radius) continue;
            float w = (1.0f - dist / (float)brush_radius) / total_weight;
            heightmap[gy * grid_size + gx] += amount * w;
        }
    }
}

// ---- Bilinear height sampling ----

static float bilinear_sample(const float* heightmap, int grid_size, float px, float py) {
    int ix = (int)px, iy = (int)py;
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (ix >= grid_size - 1) ix = grid_size - 2;
    if (iy >= grid_size - 1) iy = grid_size - 2;
    float fx = px - (float)ix;
    float fy = py - (float)iy;
    if (fx < 0.0f) fx = 0.0f;
    if (fy < 0.0f) fy = 0.0f;
    float h00 = heightmap[iy * grid_size + ix];
    float h10 = heightmap[iy * grid_size + ix + 1];
    float h01 = heightmap[(iy + 1) * grid_size + ix];
    float h11 = heightmap[(iy + 1) * grid_size + ix + 1];
    return h00 * (1.0f - fx) * (1.0f - fy)
         + h10 * fx * (1.0f - fy)
         + h01 * (1.0f - fx) * fy
         + h11 * fx * fy;
}

// ---- Core droplet simulation ----

static void run_erosion_droplets(float* heightmap, int grid_size,
                                  const ErosionParams* p, ErosionRng* rng) {
    int interior_start = EROSION_OVERLAP;
    int interior_end = EROSION_OVERLAP + EROSION_REGION_SAMPLES - 1;

    for (int i = 0; i < p->iterations; i++) {
        // Random start within interior region
        float px = erosion_rng_float(rng) * (float)(interior_end - interior_start) + (float)interior_start;
        float py = erosion_rng_float(rng) * (float)(interior_end - interior_start) + (float)interior_start;

        float dir_x = 0.0f, dir_y = 0.0f;
        float speed = 1.0f;
        float water = 1.0f;
        float sediment = 0.0f;

        for (int step = 0; step < p->max_lifetime; step++) {
            int ix = (int)px, iy = (int)py;
            if (ix < 1 || ix >= grid_size - 2 || iy < 1 || iy >= grid_size - 2)
                break;

            float fx = px - (float)ix;
            float fy = py - (float)iy;

            // Bilinear height at current position
            float h00 = heightmap[iy * grid_size + ix];
            float h10 = heightmap[iy * grid_size + ix + 1];
            float h01 = heightmap[(iy + 1) * grid_size + ix];
            float h11 = heightmap[(iy + 1) * grid_size + ix + 1];
            float old_h = h00 * (1.0f - fx) * (1.0f - fy)
                        + h10 * fx * (1.0f - fy)
                        + h01 * (1.0f - fx) * fy
                        + h11 * fx * fy;

            // Compute gradient
            float grad_x = (h10 - h00) * (1.0f - fy) + (h11 - h01) * fy;
            float grad_y = (h01 - h00) * (1.0f - fx) + (h11 - h10) * fx;

            // Update direction with inertia
            dir_x = dir_x * p->inertia - grad_x * (1.0f - p->inertia);
            dir_y = dir_y * p->inertia - grad_y * (1.0f - p->inertia);

            float len = sqrtf(dir_x * dir_x + dir_y * dir_y);
            if (len > 0.0001f) {
                dir_x /= len;
                dir_y /= len;
            } else {
                // Random direction to prevent stagnation
                float angle = erosion_rng_float(rng) * 2.0f * (float)M_PI;
                dir_x = cosf(angle);
                dir_y = sinf(angle);
            }

            // Move droplet
            float new_px = px + dir_x;
            float new_py = py + dir_y;

            // Check bounds
            if (new_px < 1 || new_px >= grid_size - 2 || new_py < 1 || new_py >= grid_size - 2)
                break;

            float new_h = bilinear_sample(heightmap, grid_size, new_px, new_py);
            float h_diff = new_h - old_h;

            // Sediment capacity
            float capacity = fmaxf(-h_diff, p->min_slope) * speed * water * p->sediment_capacity;

            if (sediment > capacity || h_diff > 0.0f) {
                // Deposit sediment
                float deposit;
                if (h_diff > 0.0f) {
                    deposit = fminf(h_diff, sediment);
                } else {
                    deposit = (sediment - capacity) * p->deposition_rate;
                }
                sediment -= deposit;
                apply_brush(heightmap, grid_size, px, py, deposit, p->brush_radius);
            } else {
                // Erode terrain
                float erode_amount = fminf((capacity - sediment) * p->erosion_rate, -h_diff);
                sediment += erode_amount;
                apply_brush(heightmap, grid_size, px, py, -erode_amount, p->brush_radius);
            }

            // Update speed and water
            speed = sqrtf(fmaxf(speed * speed + h_diff * p->gravity, 0.0f));
            water *= (1.0f - p->evaporation_rate);

            px = new_px;
            py = new_py;
        }
    }
}

// ---- Erosion job worker (runs on job system thread) ----

static void erosion_job_worker(void* data) {
    ErosionJob* job = (ErosionJob*)data;
    int N = EROSION_TOTAL_SAMPLES;

    // Allocate working heightmap
    float* heightmap = (float*)malloc((size_t)N * N * sizeof(float));
    float* raw_copy = (float*)malloc((size_t)N * N * sizeof(float));
    if (!heightmap || !raw_copy) {
        free(heightmap);
        free(raw_copy);
        job->completed = 1;
        return;
    }

    // Phase 1: Sample noise heights onto the grid
    fnl_state cont = ht_create_continental_noise(job->seed);
    fnl_state mtn  = ht_create_mountain_noise(job->seed);
    fnl_state wrp  = ht_create_warp_noise(job->seed);
    fnl_state dtl  = ht_create_detail_noise(job->seed);

    float half = (float)N * 0.5f;
    float pr = job->planet_radius;

    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            float dx = ((float)x - half) * EROSION_SAMPLE_SPACING;
            float dz = ((float)y - half) * EROSION_SAMPLE_SPACING;
            HMM_Vec3 world_pt = vec3_add(
                vec3_scale(job->center_unit, pr),
                vec3_add(vec3_scale(job->east, dx),
                         vec3_scale(job->north, dz)));
            HMM_Vec3 unit = vec3_normalize(world_pt);
            heightmap[y * N + x] = ht_sample_height_m(&cont, &mtn, &wrp, &dtl, unit);
        }
    }

    memcpy(raw_copy, heightmap, (size_t)N * N * sizeof(float));

    // Phase 2: Run droplet simulation
    ErosionRng rng = erosion_rng_create(job->seed, job->key.lon_idx, job->key.lat_idx);
    run_erosion_droplets(heightmap, N, &job->params, &rng);

    // Phase 3: Extract interior deltas
    int R = EROSION_REGION_SAMPLES;
    job->delta_out = (float*)malloc((size_t)R * R * sizeof(float));
    if (job->delta_out) {
        for (int y = 0; y < R; y++) {
            for (int x = 0; x < R; x++) {
                int sy = y + EROSION_OVERLAP;
                int sx = x + EROSION_OVERLAP;
                job->delta_out[y * R + x] = heightmap[sy * N + sx] - raw_copy[sy * N + sx];
            }
        }
    }

    free(heightmap);
    free(raw_copy);
    job->completed = 1;
}

// ---- LRU cache management ----

static int erosion_find_slot(const ErosionCache* cache, const ErosionRegionKey* key) {
    for (int i = 0; i < EROSION_MAX_CACHED; i++) {
        const ErosionRegion* r = &cache->regions[i];
        if (r->state != EROSION_STATE_EMPTY &&
            r->key.lon_idx == key->lon_idx &&
            r->key.lat_idx == key->lat_idx) {
            return i;
        }
    }
    return -1;
}

static int erosion_find_or_evict_slot(ErosionCache* cache) {
    // Find an empty slot first
    for (int i = 0; i < EROSION_MAX_CACHED; i++) {
        if (cache->regions[i].state == EROSION_STATE_EMPTY)
            return i;
    }

    // Evict least-recently-used READY slot
    int best = -1;
    uint64_t oldest = UINT64_MAX;
    for (int i = 0; i < EROSION_MAX_CACHED; i++) {
        if (cache->regions[i].state == EROSION_STATE_READY &&
            cache->regions[i].last_access_frame < oldest) {
            oldest = cache->regions[i].last_access_frame;
            best = i;
        }
    }

    if (best >= 0) {
        // Free the evicted region's data
        free(cache->regions[best].delta_heights);
        cache->regions[best].delta_heights = NULL;
        cache->regions[best].state = EROSION_STATE_EMPTY;
    }

    return best;
}

static void erosion_init_region(ErosionRegion* r, const ErosionRegionKey* key, float planet_radius) {
    memset(r, 0, sizeof(ErosionRegion));
    r->key = *key;
    r->center_unit = erosion_key_to_center(key, planet_radius);
    erosion_compute_tangent_frame(r);
}

static void erosion_submit_job(ErosionCache* cache, ErosionRegion* r, JobSystem* jobs) {
    ErosionJob* job = (ErosionJob*)calloc(1, sizeof(ErosionJob));
    if (!job) return;

    job->key = r->key;
    job->center_unit = r->center_unit;
    job->east = r->east;
    job->north = r->north;
    job->up = r->up;
    job->planet_radius = cache->planet_radius;
    job->seed = cache->seed;
    job->params = cache->params;
    job->completed = 0;
    job->region_slot = (int)(r - cache->regions);

    r->state = EROSION_STATE_COMPUTING;
    r->job_completed = 0;
    r->pending_job = job;

    if (!job_system_try_submit(jobs, erosion_job_worker, job)) {
        // Queue full — clean up and leave region empty
        free(job);
        r->state = EROSION_STATE_EMPTY;
        r->pending_job = NULL;
    }
}

// ---- Public API ----

ErosionParams erosion_default_params(void) {
    ErosionParams p;
    p.iterations = 5000;
    p.max_lifetime = 80;
    p.brush_radius = 3;
    p.inertia = 0.3f;
    p.sediment_capacity = 8.0f;
    p.erosion_rate = 0.7f;
    p.deposition_rate = 0.02f;
    p.evaporation_rate = 0.02f;
    p.min_slope = 0.01f;
    p.gravity = 10.0f;
    return p;
}

void erosion_set_params(ErosionCache* cache, const ErosionParams* params) {
    cache->params = *params;
    cache->params_hash = fnv1a_hash(params, sizeof(ErosionParams));
}

void erosion_init(ErosionCache* cache, int seed, float planet_radius, const char* cache_dir) {
    memset(cache, 0, sizeof(ErosionCache));
    cache->seed = seed;
    cache->planet_radius = planet_radius;
    snprintf(cache->cache_dir, sizeof(cache->cache_dir), "%s", cache_dir);

    ErosionParams p = erosion_default_params();
    erosion_set_params(cache, &p);

    erosion_ensure_dir(cache_dir);

    printf("[EROSION] Initialized: seed=%d, planet_r=%.0f, region=%.0fm (%dx%d @ %.1fm)\n",
           seed, planet_radius, EROSION_REGION_SIZE_M,
           EROSION_REGION_SAMPLES, EROSION_REGION_SAMPLES, EROSION_SAMPLE_SPACING);
    fflush(stdout);
}

void erosion_destroy(ErosionCache* cache) {
    for (int i = 0; i < EROSION_MAX_CACHED; i++) {
        ErosionRegion* r = &cache->regions[i];
        if (r->state == EROSION_STATE_COMPUTING && r->pending_job) {
            // Spin-wait for completion
            ErosionJob* job = (ErosionJob*)r->pending_job;
            while (!job->completed) { /* spin */ }
            free(job->delta_out);
            free(job);
        }
        free(r->delta_heights);
        r->delta_heights = NULL;
        r->state = EROSION_STATE_EMPTY;
    }
}

void erosion_update(ErosionCache* cache, HMM_Vec3 camera_unit_pos, JobSystem* jobs) {
    cache->current_frame++;
    cache->new_region_ready = false;

    // Poll computing regions for completion
    for (int i = 0; i < EROSION_MAX_CACHED; i++) {
        ErosionRegion* r = &cache->regions[i];
        if (r->state != EROSION_STATE_COMPUTING || !r->pending_job) continue;

        ErosionJob* job = (ErosionJob*)r->pending_job;
        if (!job->completed) continue;

        // Transfer results
        size_t data_size = (size_t)EROSION_REGION_SAMPLES * EROSION_REGION_SAMPLES * sizeof(float);
        if (job->delta_out) {
            r->delta_heights = job->delta_out;
            job->delta_out = NULL;  // ownership transferred
            r->state = EROSION_STATE_READY;
            r->last_access_frame = cache->current_frame;
            cache->new_region_ready = true;

            // Save to disk
            erosion_save_region(cache, r);

            printf("[EROSION] Region (%d,%d) computed and cached (%.1fKB)\n",
                   r->key.lon_idx, r->key.lat_idx, (float)data_size / 1024.0f);
            fflush(stdout);
        } else {
            r->state = EROSION_STATE_EMPTY;
        }

        free(job);
        r->pending_job = NULL;
    }

    // Prefetch 3x3 neighborhood around camera.
    // If jobs == NULL, skip submitting new work (terrain still loading).
    if (jobs) {
        ErosionRegionKey center_key = erosion_region_key(camera_unit_pos, cache->planet_radius);
        int jobs_submitted = 0;

        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                ErosionRegionKey key;
                key.lon_idx = center_key.lon_idx + dx;
                key.lat_idx = center_key.lat_idx + dy;

                // Already cached?
                int slot = erosion_find_slot(cache, &key);
                if (slot >= 0) {
                    cache->regions[slot].last_access_frame = cache->current_frame;
                    continue;
                }

                // Rate limit: 1 new erosion job per frame
                if (jobs_submitted >= 1) continue;

                // Find or evict a slot
                int free_slot = erosion_find_or_evict_slot(cache);
                if (free_slot < 0) continue;  // all slots occupied by COMPUTING regions

                ErosionRegion* r = &cache->regions[free_slot];
                erosion_init_region(r, &key, cache->planet_radius);
                r->last_access_frame = cache->current_frame;

                // Try disk cache first
                if (erosion_load_region(cache, r)) {
                    r->state = EROSION_STATE_READY;
                    cache->new_region_ready = true;
                } else {
                    // Submit compute job
                    erosion_submit_job(cache, r, jobs);
                    jobs_submitted++;
                }
            }
        }
    }
}

float erosion_sample_delta(const ErosionCache* cache, HMM_Vec3 unit_pos) {
    ErosionRegionKey key = erosion_region_key(unit_pos, cache->planet_radius);

    for (int i = 0; i < EROSION_MAX_CACHED; i++) {
        const ErosionRegion* r = &cache->regions[i];
        if (r->state != EROSION_STATE_READY) continue;
        if (r->key.lon_idx != key.lon_idx || r->key.lat_idx != key.lat_idx) continue;
        if (!r->delta_heights) continue;

        // Project unit_pos onto region's tangent plane
        HMM_Vec3 offset = vec3_sub(
            vec3_scale(unit_pos, cache->planet_radius),
            vec3_scale(r->center_unit, cache->planet_radius));
        float dx = vec3_dot(offset, r->east);
        float dz = vec3_dot(offset, r->north);

        // Convert to grid coordinates (0..EROSION_REGION_SAMPLES-1)
        float half_size = EROSION_REGION_SIZE_M * 0.5f;
        float gx = (dx + half_size) / EROSION_SAMPLE_SPACING;
        float gy = (dz + half_size) / EROSION_SAMPLE_SPACING;

        int ix = (int)gx, iy = (int)gy;
        if (ix < 0 || ix >= EROSION_REGION_SAMPLES - 1 ||
            iy < 0 || iy >= EROSION_REGION_SAMPLES - 1) {
            return 0.0f;
        }

        // Bilinear interpolation of delta
        float fx = gx - (float)ix;
        float fy = gy - (float)iy;
        int R = EROSION_REGION_SAMPLES;
        float d00 = r->delta_heights[iy * R + ix];
        float d10 = r->delta_heights[iy * R + ix + 1];
        float d01 = r->delta_heights[(iy + 1) * R + ix];
        float d11 = r->delta_heights[(iy + 1) * R + ix + 1];

        return d00 * (1.0f - fx) * (1.0f - fy)
             + d10 * fx * (1.0f - fy)
             + d01 * (1.0f - fx) * fy
             + d11 * fx * fy;
    }

    return 0.0f;
}

void erosion_precompute(ErosionCache* cache, HMM_Vec3 center_unit,
                        JobSystem* jobs, int radius) {
    ErosionRegionKey center_key = erosion_region_key(center_unit, cache->planet_radius);
    int side = 2 * radius + 1;
    int total = side * side;
    int loaded = 0, submitted = 0, computed = 0;

    printf("[EROSION] Precomputing %d regions (radius=%d) around (%d,%d)...\n",
           total, radius, center_key.lon_idx, center_key.lat_idx);
    fflush(stdout);

    // Phase 1: Load from disk or submit compute jobs
    for (int dy = -radius; dy <= radius; dy++) {
        for (int dx = -radius; dx <= radius; dx++) {
            ErosionRegionKey key;
            key.lon_idx = center_key.lon_idx + dx;
            key.lat_idx = center_key.lat_idx + dy;

            // Already in cache?
            int slot = erosion_find_slot(cache, &key);
            if (slot >= 0 && cache->regions[slot].state == EROSION_STATE_READY) {
                loaded++;
                continue;
            }

            // Need a slot
            int free_slot = erosion_find_or_evict_slot(cache);
            if (free_slot < 0) {
                printf("[EROSION]   Warning: no free cache slots\n");
                continue;
            }

            ErosionRegion* r = &cache->regions[free_slot];
            erosion_init_region(r, &key, cache->planet_radius);
            r->last_access_frame = cache->current_frame;

            // Try disk cache first
            if (erosion_load_region(cache, r)) {
                r->state = EROSION_STATE_READY;
                loaded++;
            } else {
                // Submit compute job
                erosion_submit_job(cache, r, jobs);
                submitted++;
            }
        }
    }

    printf("[EROSION]   %d loaded from disk, %d computing...\n", loaded, submitted);
    fflush(stdout);

    // Phase 2: Wait for all computing jobs to complete
    while (computed < submitted) {
        for (int i = 0; i < EROSION_MAX_CACHED; i++) {
            ErosionRegion* r = &cache->regions[i];
            if (r->state != EROSION_STATE_COMPUTING || !r->pending_job) continue;

            ErosionJob* job = (ErosionJob*)r->pending_job;
            if (!job->completed) continue;

            // Transfer results
            if (job->delta_out) {
                r->delta_heights = job->delta_out;
                job->delta_out = NULL;
                r->state = EROSION_STATE_READY;
                r->last_access_frame = cache->current_frame;

                erosion_save_region(cache, r);
                computed++;

                printf("[EROSION]   Computed region (%d,%d) — %d/%d done\n",
                       r->key.lon_idx, r->key.lat_idx,
                       loaded + computed, total);
                fflush(stdout);
            } else {
                r->state = EROSION_STATE_EMPTY;
                computed++;
            }

            free(job);
            r->pending_job = NULL;
        }

        if (computed < submitted) {
#ifdef _WIN32
            Sleep(1);
#else
            usleep(1000);
#endif
        }
    }

    printf("[EROSION] Precompute complete: %d loaded, %d computed (%d total)\n",
           loaded, computed, loaded + computed);
    fflush(stdout);
}

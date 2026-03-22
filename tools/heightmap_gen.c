// Quick tool to generate an equirectangular heightmap of the planet
// Compile: cl /I ../third_party /I ../src tools/heightmap_gen.c
// Output: heightmap.ppm (color-coded biome map with sun terminator)

#include "HandmadeMath.h"
#include "terrain_noise.h"
// FNL_IMPL after terrain_noise.h so the impl is only emitted once
#define FNL_IMPL
#include "FastNoiseLite.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define WIDTH  1024
#define HEIGHT 512
#define SEED   42
#define PI     3.14159265358979323846

int main(void) {
    fnl_state continental = ht_create_continental_noise(SEED);
    fnl_state mountain    = ht_create_mountain_noise(SEED);
    fnl_state warp        = ht_create_warp_noise(SEED);
    fnl_state detail      = ht_create_detail_noise(SEED);

    // Sun direction (matches render.c)
    HMM_Vec3 sun = HMM_NormV3((HMM_Vec3){{1.0f, 0.3f, 0.5f}});

    unsigned char* img = malloc(WIDTH * HEIGHT * 3);

    float best_grass_twilight_score = -1e9f;
    float best_lat = 0, best_lon = 0;
    float best_height = 0;

    for (int y = 0; y < HEIGHT; y++) {
        float lat = PI * (0.5f - (float)y / HEIGHT);  // +90 to -90
        for (int x = 0; x < WIDTH; x++) {
            float lon = 2.0f * PI * ((float)x / WIDTH - 0.5f);  // -180 to +180

            float cos_lat = cosf(lat);
            HMM_Vec3 pos = {{cos_lat * cosf(lon), sinf(lat), cos_lat * sinf(lon)}};

            float h_m = ht_sample_height_m(&continental, &mountain, &warp, &detail, pos);
            float rel = h_m - TERRAIN_SEA_LEVEL_M;

            // Biome color
            unsigned char r, g, b;
            if (h_m < TERRAIN_SEA_LEVEL_M) {
                // Water - depth shading
                float depth = (TERRAIN_SEA_LEVEL_M - h_m) / 3500.0f;
                if (depth > 1.0f) depth = 1.0f;
                r = (unsigned char)(30 - 20 * depth);
                g = (unsigned char)(80 - 40 * depth);
                b = (unsigned char)(180 - 60 * depth);
            } else if (rel < 200.0f) {
                r = 210; g = 190; b = 130;  // Sand
            } else if (rel < 1500.0f) {
                // Grass - shade by elevation
                float t = (rel - 200.0f) / 1300.0f;
                r = (unsigned char)(60 + 30 * t);
                g = (unsigned char)(140 - 30 * t);
                b = (unsigned char)(40 + 20 * t);
            } else if (rel < 4500.0f) {
                // Stone
                float t = (rel - 1500.0f) / 3000.0f;
                r = (unsigned char)(120 + 30 * t);
                g = (unsigned char)(110 + 20 * t);
                b = (unsigned char)(100 + 20 * t);
            } else {
                r = 240; g = 245; b = 255;  // Ice
            }

            // Sun illumination overlay
            float ndotl = HMM_DotV3(pos, sun);
            float brightness = 0.3f + 0.7f * fmaxf(0.0f, ndotl);
            r = (unsigned char)(r * brightness);
            g = (unsigned char)(g * brightness);
            b = (unsigned char)(b * brightness);

            // Draw terminator line (where ndotl ~ 0)
            if (fabsf(ndotl) < 0.015f) {
                r = 255; g = 50; b = 50;  // Red terminator line
            }

            // Score grass+twilight candidates
            if (rel >= 300.0f && rel <= 1000.0f) {
                float twilight_score = 1.0f - fabsf(ndotl) * 5.0f;  // best near terminator
                if (twilight_score > best_grass_twilight_score) {
                    best_grass_twilight_score = twilight_score;
                    best_lat = lat;
                    best_lon = lon;
                    best_height = h_m;
                }
            }

            int idx = (y * WIDTH + x) * 3;
            img[idx] = r; img[idx+1] = g; img[idx+2] = b;
        }
    }

    // Mark best spawn on map (green crosshair)
    int sx = (int)((best_lon / (2.0f * PI) + 0.5f) * WIDTH);
    int sy = (int)((0.5f - best_lat / PI) * HEIGHT);
    for (int d = -8; d <= 8; d++) {
        if (sx+d >= 0 && sx+d < WIDTH)  { int i = (sy*WIDTH + sx+d)*3; img[i]=0; img[i+1]=255; img[i+2]=0; }
        if (sy+d >= 0 && sy+d < HEIGHT) { int i = ((sy+d)*WIDTH + sx)*3; img[i]=0; img[i+1]=255; img[i+2]=0; }
    }

    // Write PPM
    FILE* f = fopen("heightmap.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(img, 1, WIDTH * HEIGHT * 3, f);
    fclose(f);
    free(img);

    // Convert spawn to unit sphere direction
    float cos_lat = cosf(best_lat);
    float spawn_x = cos_lat * cosf(best_lon);
    float spawn_y = sinf(best_lat);
    float spawn_z = cos_lat * sinf(best_lon);

    printf("=== Best Spawn Point ===\n");
    printf("Lat: %.2f deg, Lon: %.2f deg\n", best_lat * 180.0f / PI, best_lon * 180.0f / PI);
    printf("Height: %.0f m (%.0f m above sea level)\n", best_height, best_height - TERRAIN_SEA_LEVEL_M);
    printf("Unit sphere dir: (%.6f, %.6f, %.6f)\n", spawn_x, spawn_y, spawn_z);
    printf("Sun dot: %.4f (0 = terminator)\n", spawn_x*sun.X + spawn_y*sun.Y + spawn_z*sun.Z);
    printf("Heightmap written to heightmap.ppm\n");

    return 0;
}

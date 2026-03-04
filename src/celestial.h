#ifndef CELESTIAL_H
#define CELESTIAL_H

#include <stdbool.h>
#include "HandmadeMath.h"
#include "sokol_gfx.h"

/* Forward declarations */
typedef struct Camera Camera;

#define MAX_MOONS 10
#define MOON_MESH_SUBDIVISIONS 4   /* icosphere subdivision depth (4 = ~2562 verts, ~5120 tris) */
#define MOON_SOI_RADIUS 50000.0f   /* 50 km sphere of influence */

typedef struct MoonShapeParams {
    float base_radius;             /* meters (1000-10000) */
    float ellipsoid_scale[3];      /* axis stretch (e.g., {1.0, 0.8, 1.2}) */
    int noise_seed;
    float noise_frequency;         /* 0.5-2.0 */
    float noise_amplitude;         /* fraction of radius (0.05-0.2) */
    int noise_octaves;             /* 2-4 */
} MoonShapeParams;

typedef struct MoonColorPalette {
    HMM_Vec3 base_color;
    HMM_Vec3 highlight_color;      /* tops / high elevations */
    HMM_Vec3 shadow_color;         /* crevices / low elevations */
} MoonColorPalette;

typedef struct OrbitParams {
    float semi_major_axis;         /* distance from Tenebris center (m) */
    float eccentricity;            /* 0.0-0.1 (near-circular) */
    float inclination;             /* radians */
    float longitude_ascending;     /* radians (Omega) */
    float argument_periapsis;      /* radians (omega) */
    float mean_anomaly_epoch;      /* radians at t=0 */
    float period;                  /* orbital period in seconds */
} OrbitParams;

typedef struct CelestialBody {
    char name[32];
    double pos_d[3];               /* double-precision world position */
    HMM_Vec3 position;             /* float copy for rendering */
    float radius;                  /* effective radius (for SOI checks) */
    MoonShapeParams shape;
    MoonColorPalette palette;
    OrbitParams orbit;

    /* GPU mesh — vertices stored in moon-local coords (centered at origin) */
    sg_buffer gpu_buffer;
    int vertex_count;
    bool mesh_ready;
} CelestialBody;

typedef struct SolarSystem {
    CelestialBody moons[MAX_MOONS];
    int moon_count;
    double elapsed_time;           /* accumulated simulation time (seconds) */
    int gravity_body;              /* -1 = main planet, 0-9 = moon index */
} SolarSystem;

/* Initialize all 10 moons with hardcoded orbital + shape params */
void solar_system_init(SolarSystem* ss);

/* Update moon positions via Kepler solver */
void solar_system_update(SolarSystem* ss, double dt);

/* Find which body's surface the player is within 50km of.
   Returns -1 for main planet, 0-9 for moon index.
   If >50km from everything, returns -1 (main planet fallback). */
int solar_system_find_gravity_body(const SolarSystem* ss,
                                   const double pos_d[3],
                                   float planet_radius);

/* Get surface radius at a given direction on a moon (ellipsoid + noise) */
float moon_surface_radius(const MoonShapeParams* shape, HMM_Vec3 unit_dir);

/* Generate icosphere mesh for a moon and upload to GPU */
void moon_generate_mesh(CelestialBody* body);

/* Render all moons using the planet pipeline (camera-relative) */
void solar_system_render(const SolarSystem* ss,
                         const Camera* cam,
                         HMM_Vec3 sun_dir,
                         HMM_Mat4 vp_rot,
                         sg_pipeline pip,
                         float Fcoef, float far_plane, float z_bias,
                         const double world_origin[3]);

/* Draw moon name labels in screen space (call between sdtx_canvas and sdtx_draw).
   Only draws when camera is in space_mode. */
void solar_system_draw_labels(const SolarSystem* ss,
                              const Camera* cam,
                              HMM_Mat4 vp);

/* Cleanup GPU resources */
void solar_system_shutdown(SolarSystem* ss);

#endif

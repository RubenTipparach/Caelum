#include "celestial.h"
#include "camera.h"
#include "math_utils.h"
#include "FastNoiseLite.h"
#include "planet.glsl.h"
#include "sokol_app.h"
#include "util/sokol_debugtext.h"
#include "util/sokol_gl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- LodVertex (must match planet.glsl layout: pos[3] + normal[3] + color[3]) ---- */
typedef struct MoonVertex {
    float pos[3];
    float normal[3];
    float color[3];
} MoonVertex;

/* ---- Icosphere generation scratch buffers ---- */
/* Max verts after 4 subdivisions: 12 + edges*splits
   For 4 subdivisions: 10*4^4 + 2 = 2562 verts, 20*4^4 = 5120 tris */
#define MOON_MAX_VERTS  2600
#define MOON_MAX_TRIS   5200

typedef struct IcoMesh {
    HMM_Vec3 verts[MOON_MAX_VERTS];
    int tris[MOON_MAX_TRIS][3];
    int vert_count;
    int tri_count;
} IcoMesh;

/* Edge midpoint cache for subdivision (key = packed pair of vertex indices) */
typedef struct EdgeEntry {
    int v0, v1;     /* sorted vertex indices */
    int mid;        /* midpoint vertex index */
} EdgeEntry;

#define EDGE_CACHE_SIZE 16384
static EdgeEntry edge_cache[EDGE_CACHE_SIZE];
static int edge_cache_count;

static void edge_cache_clear(void) {
    edge_cache_count = 0;
    memset(edge_cache, 0xFF, sizeof(edge_cache));  /* -1 sentinel */
}

static int edge_cache_get_or_create(IcoMesh* mesh, int v0, int v1) {
    /* Sort indices for consistent key */
    int a = v0 < v1 ? v0 : v1;
    int b = v0 < v1 ? v1 : v0;
    unsigned int key = ((unsigned int)a * 7919 + (unsigned int)b) % EDGE_CACHE_SIZE;

    /* Linear probe */
    for (int i = 0; i < 8; i++) {
        int idx = (int)((key + (unsigned int)i) % EDGE_CACHE_SIZE);
        if (edge_cache[idx].v0 == a && edge_cache[idx].v1 == b) {
            return edge_cache[idx].mid;
        }
        if (edge_cache[idx].v0 == -1) {
            /* Not found — create midpoint */
            HMM_Vec3 mid = HMM_NormV3(HMM_AddV3(mesh->verts[a], mesh->verts[b]));
            int mid_idx = mesh->vert_count++;
            mesh->verts[mid_idx] = mid;
            edge_cache[idx].v0 = a;
            edge_cache[idx].v1 = b;
            edge_cache[idx].mid = mid_idx;
            return mid_idx;
        }
    }
    /* Cache full (shouldn't happen with 16k entries) — just create */
    HMM_Vec3 mid = HMM_NormV3(HMM_AddV3(mesh->verts[v0], mesh->verts[v1]));
    int mid_idx = mesh->vert_count++;
    mesh->verts[mid_idx] = mid;
    return mid_idx;
}

static void ico_subdivide(IcoMesh* mesh) {
    edge_cache_clear();
    int old_tri_count = mesh->tri_count;
    int new_tris[MOON_MAX_TRIS][3];
    int new_tri_count = 0;

    for (int i = 0; i < old_tri_count; i++) {
        int v0 = mesh->tris[i][0];
        int v1 = mesh->tris[i][1];
        int v2 = mesh->tris[i][2];

        int a = edge_cache_get_or_create(mesh, v0, v1);
        int b = edge_cache_get_or_create(mesh, v1, v2);
        int c = edge_cache_get_or_create(mesh, v2, v0);

        new_tris[new_tri_count][0] = v0; new_tris[new_tri_count][1] = a; new_tris[new_tri_count][2] = c; new_tri_count++;
        new_tris[new_tri_count][0] = a;  new_tris[new_tri_count][1] = v1; new_tris[new_tri_count][2] = b; new_tri_count++;
        new_tris[new_tri_count][0] = c;  new_tris[new_tri_count][1] = b; new_tris[new_tri_count][2] = v2; new_tri_count++;
        new_tris[new_tri_count][0] = a;  new_tris[new_tri_count][1] = b; new_tris[new_tri_count][2] = c; new_tri_count++;
    }

    mesh->tri_count = new_tri_count;
    memcpy(mesh->tris, new_tris, (size_t)new_tri_count * sizeof(int) * 3);
}

/* ---- Moon surface radius at a unit direction ---- */
float moon_surface_radius(const MoonShapeParams* shape, HMM_Vec3 unit_dir) {
    /* Apply ellipsoid scaling */
    HMM_Vec3 scaled = {{
        unit_dir.X * shape->ellipsoid_scale[0],
        unit_dir.Y * shape->ellipsoid_scale[1],
        unit_dir.Z * shape->ellipsoid_scale[2]
    }};
    float scale_len = HMM_LenV3(scaled);
    if (scale_len < 1e-6f) scale_len = 1.0f;

    /* Noise displacement */
    fnl_state noise = fnlCreateState();
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = shape->noise_octaves;
    noise.frequency = shape->noise_frequency;
    noise.seed = shape->noise_seed;

    float n = fnlGetNoise3D(&noise,
        unit_dir.X * 100.0f,
        unit_dir.Y * 100.0f,
        unit_dir.Z * 100.0f);

    /* Result: base radius * ellipsoid distortion * (1 + noise) */
    float ellipsoid_r = shape->base_radius / scale_len;
    return ellipsoid_r * (1.0f + n * shape->noise_amplitude);
}

/* ---- Generate icosphere mesh for a moon ---- */
void moon_generate_mesh(CelestialBody* body) {
    static IcoMesh mesh;  /* large — use static to avoid stack overflow */
    mesh.vert_count = 0;
    mesh.tri_count = 0;

    /* Step 1: Seed with icosahedron vertices (from math_utils) */
    for (int i = 0; i < ICO_VERTEX_COUNT; i++) {
        mesh.verts[mesh.vert_count++] = HMM_NormV3(ICO_VERTICES[i]);
    }
    for (int i = 0; i < ICO_FACE_COUNT; i++) {
        mesh.tris[mesh.tri_count][0] = ICO_FACES[i][0];
        mesh.tris[mesh.tri_count][1] = ICO_FACES[i][1];
        mesh.tris[mesh.tri_count][2] = ICO_FACES[i][2];
        mesh.tri_count++;
    }

    /* Step 2: Subdivide MOON_MESH_SUBDIVISIONS times */
    for (int s = 0; s < MOON_MESH_SUBDIVISIONS; s++) {
        ico_subdivide(&mesh);
    }

    printf("[CELESTIAL] Moon '%s': %d verts, %d tris after %d subdivisions\n",
           body->name, mesh.vert_count, mesh.tri_count, MOON_MESH_SUBDIVISIONS);

    /* Step 3: Displace vertices by shape (ellipsoid + noise) */
    float* radii = (float*)malloc((size_t)mesh.vert_count * sizeof(float));
    float min_r = 1e9f, max_r = -1e9f;
    for (int i = 0; i < mesh.vert_count; i++) {
        float r = moon_surface_radius(&body->shape, mesh.verts[i]);
        radii[i] = r;
        if (r < min_r) min_r = r;
        if (r > max_r) max_r = r;
        mesh.verts[i] = HMM_MulV3F(mesh.verts[i], r);
    }
    body->radius = (min_r + max_r) * 0.5f;

    /* Step 4: Compute face normals, then vertex normals by averaging */
    HMM_Vec3* vert_normals = (HMM_Vec3*)calloc((size_t)mesh.vert_count, sizeof(HMM_Vec3));
    for (int i = 0; i < mesh.tri_count; i++) {
        HMM_Vec3 v0 = mesh.verts[mesh.tris[i][0]];
        HMM_Vec3 v1 = mesh.verts[mesh.tris[i][1]];
        HMM_Vec3 v2 = mesh.verts[mesh.tris[i][2]];
        HMM_Vec3 edge1 = HMM_SubV3(v1, v0);
        HMM_Vec3 edge2 = HMM_SubV3(v2, v0);
        HMM_Vec3 face_n = HMM_Cross(edge1, edge2);
        /* Weight by face area (not normalized) — larger faces contribute more */
        vert_normals[mesh.tris[i][0]] = HMM_AddV3(vert_normals[mesh.tris[i][0]], face_n);
        vert_normals[mesh.tris[i][1]] = HMM_AddV3(vert_normals[mesh.tris[i][1]], face_n);
        vert_normals[mesh.tris[i][2]] = HMM_AddV3(vert_normals[mesh.tris[i][2]], face_n);
    }
    for (int i = 0; i < mesh.vert_count; i++) {
        vert_normals[i] = HMM_NormV3(vert_normals[i]);
    }

    /* Step 5: Color from palette based on height + normal.Y */
    float range = max_r - min_r;
    if (range < 1.0f) range = 1.0f;

    /* Step 6: Build triangle vertex buffer (non-indexed for simplicity) */
    int total_verts = mesh.tri_count * 3;
    MoonVertex* vbuf = (MoonVertex*)malloc((size_t)total_verts * sizeof(MoonVertex));
    int vi = 0;
    for (int i = 0; i < mesh.tri_count; i++) {
        for (int j = 0; j < 3; j++) {
            int idx = mesh.tris[i][j];
            HMM_Vec3 p = mesh.verts[idx];
            HMM_Vec3 n = vert_normals[idx];
            float height_t = (radii[idx] - min_r) / range;

            /* Blend base → highlight by height, darken in crevices via normal.Y */
            HMM_Vec3 col = HMM_LerpV3(body->palette.shadow_color,
                                        height_t, body->palette.base_color);
            float top_factor = n.Y * 0.5f + 0.5f;  /* remap -1..1 to 0..1 */
            col = HMM_LerpV3(col, height_t * top_factor * 0.5f,
                              body->palette.highlight_color);

            vbuf[vi].pos[0] = p.X;
            vbuf[vi].pos[1] = p.Y;
            vbuf[vi].pos[2] = p.Z;
            vbuf[vi].normal[0] = n.X;
            vbuf[vi].normal[1] = n.Y;
            vbuf[vi].normal[2] = n.Z;
            vbuf[vi].color[0] = col.X;
            vbuf[vi].color[1] = col.Y;
            vbuf[vi].color[2] = col.Z;
            vi++;
        }
    }

    /* Upload to GPU */
    body->gpu_buffer = sg_make_buffer(&(sg_buffer_desc){
        .data = { .ptr = vbuf, .size = (size_t)total_verts * sizeof(MoonVertex) },
        .label = body->name,
    });
    body->vertex_count = total_verts;
    body->mesh_ready = true;

    free(vbuf);
    free(vert_normals);
    free(radii);
}

/* ---- Tenebris fallback mesh generation ---- */
void solar_system_generate_planet_mesh(SolarSystem* ss, float planet_radius, int seed) {
    static IcoMesh mesh;
    mesh.vert_count = 0;
    mesh.tri_count = 0;

    /* Build icosphere (3 subdivisions = ~642 verts, ~1280 tris) */
    for (int i = 0; i < ICO_VERTEX_COUNT; i++) {
        mesh.verts[mesh.vert_count++] = HMM_NormV3(ICO_VERTICES[i]);
    }
    for (int i = 0; i < ICO_FACE_COUNT; i++) {
        mesh.tris[mesh.tri_count][0] = ICO_FACES[i][0];
        mesh.tris[mesh.tri_count][1] = ICO_FACES[i][1];
        mesh.tris[mesh.tri_count][2] = ICO_FACES[i][2];
        mesh.tri_count++;
    }
    for (int s = 0; s < 3; s++) {
        ico_subdivide(&mesh);
    }

    /* Displace by terrain noise (same noise as LOD system uses) */
    fnl_state noise = fnlCreateState();
    noise.seed = seed;
    noise.noise_type = FNL_NOISE_OPENSIMPLEX2S;
    noise.fractal_type = FNL_FRACTAL_FBM;
    noise.octaves = 6;
    noise.frequency = 1.5f / planet_radius;

    float* radii = (float*)malloc((size_t)mesh.vert_count * sizeof(float));
    float min_r = 1e9f, max_r = -1e9f;
    float layer_thickness = 0.5f;  /* approximate */
    int sea_level = 50;

    for (int v = 0; v < mesh.vert_count; v++) {
        HMM_Vec3 dir = mesh.verts[v];
        float n3d = fnlGetNoise3D(&noise, dir.X * planet_radius,
                                  dir.Y * planet_radius, dir.Z * planet_radius);
        float height_layers = (n3d + 1.0f) * 0.5f * 120.0f;  /* 0-120 layers */
        float r = planet_radius + (height_layers - (float)sea_level) * layer_thickness;
        if (r < planet_radius) r = planet_radius;  /* ocean floor clamp */
        radii[v] = r;
        if (r < min_r) min_r = r;
        if (r > max_r) max_r = r;
        mesh.verts[v] = HMM_MulV3F(dir, r);
    }

    /* Compute vertex normals */
    HMM_Vec3* vert_normals = (HMM_Vec3*)calloc((size_t)mesh.vert_count, sizeof(HMM_Vec3));
    for (int i = 0; i < mesh.tri_count; i++) {
        HMM_Vec3 v0 = mesh.verts[mesh.tris[i][0]];
        HMM_Vec3 v1 = mesh.verts[mesh.tris[i][1]];
        HMM_Vec3 v2 = mesh.verts[mesh.tris[i][2]];
        HMM_Vec3 face_n = HMM_Cross(HMM_SubV3(v1, v0), HMM_SubV3(v2, v0));
        vert_normals[mesh.tris[i][0]] = HMM_AddV3(vert_normals[mesh.tris[i][0]], face_n);
        vert_normals[mesh.tris[i][1]] = HMM_AddV3(vert_normals[mesh.tris[i][1]], face_n);
        vert_normals[mesh.tris[i][2]] = HMM_AddV3(vert_normals[mesh.tris[i][2]], face_n);
    }
    for (int i = 0; i < mesh.vert_count; i++) {
        vert_normals[i] = HMM_NormV3(vert_normals[i]);
    }

    /* Color by biome: latitude + elevation */
    float range = max_r - min_r;
    if (range < 1.0f) range = 1.0f;
    int total_verts = mesh.tri_count * 3;
    MoonVertex* vbuf = (MoonVertex*)malloc((size_t)total_verts * sizeof(MoonVertex));
    int vi = 0;
    for (int i = 0; i < mesh.tri_count; i++) {
        for (int j = 0; j < 3; j++) {
            int idx = mesh.tris[i][j];
            HMM_Vec3 p = mesh.verts[idx];
            HMM_Vec3 n = vert_normals[idx];
            float height_t = (radii[idx] - min_r) / range;

            /* Determine color based on height and latitude */
            HMM_Vec3 dir = HMM_NormV3(p);
            float lat = fabsf(dir.Y);  /* 0 = equator, 1 = pole */
            HMM_Vec3 col;
            if (radii[idx] <= planet_radius + 0.1f) {
                /* Ocean: deep blue */
                col = (HMM_Vec3){{0.05f, 0.12f, 0.35f}};
            } else if (lat > 0.85f) {
                /* Polar ice */
                col = (HMM_Vec3){{0.85f, 0.88f, 0.92f}};
            } else if (height_t > 0.7f) {
                /* Mountains: grey-brown */
                col = (HMM_Vec3){{0.45f, 0.42f, 0.38f}};
            } else if (lat > 0.6f) {
                /* Tundra */
                col = (HMM_Vec3){{0.35f, 0.40f, 0.30f}};
            } else {
                /* Lowland green/brown blend */
                float green = 0.3f + height_t * 0.3f;
                col = (HMM_Vec3){{0.18f + height_t * 0.2f, green, 0.12f}};
            }

            vbuf[vi].pos[0] = p.X;
            vbuf[vi].pos[1] = p.Y;
            vbuf[vi].pos[2] = p.Z;
            vbuf[vi].normal[0] = n.X;
            vbuf[vi].normal[1] = n.Y;
            vbuf[vi].normal[2] = n.Z;
            vbuf[vi].color[0] = col.X;
            vbuf[vi].color[1] = col.Y;
            vbuf[vi].color[2] = col.Z;
            vi++;
        }
    }

    ss->planet_buffer = sg_make_buffer(&(sg_buffer_desc){
        .data = { .ptr = vbuf, .size = (size_t)total_verts * sizeof(MoonVertex) },
        .label = "tenebris-fallback",
    });
    ss->planet_vertex_count = total_verts;
    ss->planet_mesh_ready = true;
    ss->planet_radius = planet_radius;

    printf("[CELESTIAL] Tenebris fallback mesh: %d verts (3 subdivisions)\n", total_verts);

    free(vbuf);
    free(vert_normals);
    free(radii);
}

/* ---- Kepler orbit solver ---- */
static void kepler_position(const OrbitParams* orbit, double t, double out_pos[3]) {
    /* Mean anomaly */
    double M = (double)orbit->mean_anomaly_epoch +
               (2.0 * M_PI / (double)orbit->period) * t;
    /* Wrap to [0, 2*PI] */
    M = fmod(M, 2.0 * M_PI);
    if (M < 0.0) M += 2.0 * M_PI;

    /* Eccentric anomaly via Newton iteration (3 steps, low-e converges fast) */
    double e = (double)orbit->eccentricity;
    double E = M;
    for (int i = 0; i < 5; i++) {
        E = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
    }

    /* True anomaly */
    double cos_E = cos(E);
    double sin_E = sin(E);
    double nu = atan2(sqrt(1.0 - e * e) * sin_E, cos_E - e);

    /* Distance from focus */
    double r = (double)orbit->semi_major_axis * (1.0 - e * cos_E);

    /* Position in orbital plane (perifocal frame) */
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);

    /* Rotate to world frame by Omega, i, omega */
    double O = (double)orbit->longitude_ascending;
    double inc = (double)orbit->inclination;
    double w = (double)orbit->argument_periapsis;

    double cos_O = cos(O), sin_O = sin(O);
    double cos_i = cos(inc), sin_i = sin(inc);
    double cos_w = cos(w), sin_w = sin(w);

    /* Rotation matrix components (perifocal → inertial) */
    double Px = cos_O * cos_w - sin_O * sin_w * cos_i;
    double Py = sin_O * cos_w + cos_O * sin_w * cos_i;
    double Pz = sin_w * sin_i;

    double Qx = -cos_O * sin_w - sin_O * cos_w * cos_i;
    double Qy = -sin_O * sin_w + cos_O * cos_w * cos_i;
    double Qz = cos_w * sin_i;

    out_pos[0] = Px * x_orb + Qx * y_orb;
    out_pos[1] = Pz * x_orb + Qz * y_orb;  /* Y-up: swap Y/Z */
    out_pos[2] = Py * x_orb + Qy * y_orb;
}

/* ---- Solar system initialization ---- */
void solar_system_init(SolarSystem* ss) {
    memset(ss, 0, sizeof(SolarSystem));
    ss->gravity_body = -1;
    ss->pinned_body = -1;

    /* Define 10 moons: 5 large (5-10km), 5 small (1-5km) */
    typedef struct {
        const char* name;
        float semi_major;   /* orbital distance from Tenebris (m) */
        float period;       /* orbital period (seconds) */
        float incl_deg;     /* orbital inclination (degrees) */
        float lon_asc_deg;  /* longitude of ascending node (degrees) */
        float ecc;          /* eccentricity */
        float base_r;       /* base radius (m) */
        float ell[3];       /* ellipsoid scale */
        int seed;
        float noise_freq;
        float noise_amp;
        int noise_oct;
        HMM_Vec3 base_col, hi_col, shadow_col;
    } MoonDef;

    /*
     * Orbital periods: Kepler's 3rd law scaling (T ~ a^1.5).
     * Innermost moon (2Mm) = 1 hour (3600s), scale outward from there.
     * Reference: a0 = 2,000,000m, T0 = 3600s.
     * T = T0 * (a/a0)^1.5
     */
    /*
     * All moons are grey rocky bodies — cool blue-gray tones to match moon.png texture.
     * Colors: base=mid blue-gray, highlight=light blue-gray (peaks), shadow=dark blue-gray (crevices).
     */
    MoonDef defs[MAX_MOONS] = {
        /* --- 2 giant moons (50-100km radius) --- */
        /* noise_amp = 0.01 for all moons (1% of radius, matching Tenebris 8km/800km ratio) */
        { "Gorrath",  5000000.0f, 14230.0f,  5.0f,   0.0f, 0.02f,
          80000.0f, {1.0f, 0.9f, 1.05f}, 100, 0.3f, 0.01f, 4,
          {{0.38f, 0.38f, 0.42f}}, {{0.52f, 0.52f, 0.57f}}, {{0.18f, 0.18f, 0.22f}} },

        { "Atheron", 10000000.0f, 40250.0f, 10.0f,  90.0f, 0.03f,
          55000.0f, {1.05f, 0.95f, 1.0f}, 150, 0.4f, 0.01f, 4,
          {{0.40f, 0.40f, 0.44f}}, {{0.55f, 0.55f, 0.60f}}, {{0.20f, 0.20f, 0.24f}} },

        /* --- 3 large moons (5-10km radius) --- */
        { "Kelthos",  3000000.0f, 6600.0f, 12.0f,  45.0f, 0.04f,
          7000.0f, {0.95f, 1.10f, 1.0f}, 300, 0.5f, 0.01f, 4,
          {{0.36f, 0.37f, 0.40f}}, {{0.50f, 0.51f, 0.55f}}, {{0.17f, 0.17f, 0.20f}} },

        { "Dravok",   6500000.0f, 21100.0f, 18.0f, 200.0f, 0.06f,
          6000.0f, {1.20f, 0.80f, 0.95f}, 400, 0.6f, 0.01f, 3,
          {{0.34f, 0.34f, 0.38f}}, {{0.48f, 0.48f, 0.53f}}, {{0.16f, 0.16f, 0.19f}} },

        { "Serath",   8500000.0f, 31560.0f,  8.0f, 300.0f, 0.03f,
          5000.0f, {1.0f, 1.0f, 1.15f}, 500, 0.4f, 0.01f, 3,
          {{0.42f, 0.42f, 0.46f}}, {{0.57f, 0.57f, 0.62f}}, {{0.22f, 0.22f, 0.25f}} },

        /* --- 5 small moons (1-5km radius) --- */
        { "Cryx",     2000000.0f, 3600.0f, 25.0f,  60.0f, 0.08f,
          3000.0f, {1.3f, 0.7f, 1.0f}, 600, 0.5f, 0.01f, 2,
          {{0.32f, 0.32f, 0.36f}}, {{0.46f, 0.46f, 0.51f}}, {{0.15f, 0.15f, 0.18f}} },

        { "Nyctra",   4500000.0f, 12150.0f, 30.0f, 150.0f, 0.05f,
          2000.0f, {0.8f, 1.2f, 1.1f}, 700, 0.6f, 0.01f, 2,
          {{0.37f, 0.38f, 0.42f}}, {{0.52f, 0.53f, 0.58f}}, {{0.18f, 0.18f, 0.22f}} },

        { "Thalwen",  7500000.0f, 26200.0f, 15.0f, 270.0f, 0.07f,
          4000.0f, {1.1f, 0.9f, 1.2f}, 800, 0.5f, 0.01f, 3,
          {{0.35f, 0.35f, 0.40f}}, {{0.50f, 0.50f, 0.56f}}, {{0.17f, 0.17f, 0.21f}} },

        { "Vexis",   11000000.0f, 46440.0f, 40.0f,  90.0f, 0.09f,
          1500.0f, {1.4f, 0.6f, 1.0f}, 900, 0.7f, 0.01f, 2,
          {{0.33f, 0.33f, 0.38f}}, {{0.47f, 0.47f, 0.53f}}, {{0.15f, 0.15f, 0.19f}} },

        { "Zephyros", 13000000.0f, 59720.0f, 10.0f, 330.0f, 0.02f,
          1000.0f, {1.0f, 1.3f, 0.8f}, 1000, 0.7f, 0.01f, 2,
          {{0.39f, 0.40f, 0.44f}}, {{0.54f, 0.55f, 0.60f}}, {{0.19f, 0.20f, 0.23f}} },
    };

    ss->moon_count = MAX_MOONS;

    for (int i = 0; i < MAX_MOONS; i++) {
        CelestialBody* m = &ss->moons[i];
        MoonDef* d = &defs[i];
        memset(m, 0, sizeof(CelestialBody));

        /* Copy name */
        snprintf(m->name, sizeof(m->name), "%s", d->name);

        m->shape.base_radius = d->base_r;
        m->shape.ellipsoid_scale[0] = d->ell[0];
        m->shape.ellipsoid_scale[1] = d->ell[1];
        m->shape.ellipsoid_scale[2] = d->ell[2];
        m->shape.noise_seed = d->seed;
        m->shape.noise_frequency = d->noise_freq;
        m->shape.noise_amplitude = d->noise_amp;
        m->shape.noise_octaves = d->noise_oct;

        m->palette.base_color = d->base_col;
        m->palette.highlight_color = d->hi_col;
        m->palette.shadow_color = d->shadow_col;

        m->orbit.semi_major_axis = d->semi_major;
        m->orbit.eccentricity = d->ecc;
        m->orbit.inclination = d->incl_deg * (float)(M_PI / 180.0);
        m->orbit.longitude_ascending = d->lon_asc_deg * (float)(M_PI / 180.0);
        m->orbit.argument_periapsis = (float)(i * 36.0 * M_PI / 180.0);  /* Spread evenly */
        m->orbit.mean_anomaly_epoch = (float)(i * 0.628);  /* Spread initial positions */
        m->orbit.period = d->period;

        m->radius = d->base_r;  /* Will be updated after mesh gen */

        /* Initial position */
        kepler_position(&m->orbit, 0.0, m->pos_d);
        m->prev_pos_d[0] = m->pos_d[0];
        m->prev_pos_d[1] = m->pos_d[1];
        m->prev_pos_d[2] = m->pos_d[2];
        m->position = (HMM_Vec3){{(float)m->pos_d[0], (float)m->pos_d[1], (float)m->pos_d[2]}};

        /* Generate mesh */
        moon_generate_mesh(m);

        printf("[CELESTIAL] Moon '%s': orbit=%.0fkm radius=%.0fm\n",
               m->name, d->semi_major / 1000.0f, m->radius);
    }

    /* Create sgl pipeline for depth-tested orbit lines */
    ss->orbit_pip = sgl_make_pipeline(&(sg_pipeline_desc){
        .depth = {
            .compare = SG_COMPAREFUNC_LESS_EQUAL,
            .write_enabled = true,
        },
        .colors[0].blend = {
            .enabled = true,
            .src_factor_rgb = SG_BLENDFACTOR_SRC_ALPHA,
            .dst_factor_rgb = SG_BLENDFACTOR_ONE_MINUS_SRC_ALPHA,
        },
    });
}

void solar_system_update(SolarSystem* ss, double dt) {
    ss->elapsed_time += dt;
    for (int i = 0; i < ss->moon_count; i++) {
        CelestialBody* m = &ss->moons[i];
        /* Save previous position before orbital update (for camera tracking delta) */
        m->prev_pos_d[0] = m->pos_d[0];
        m->prev_pos_d[1] = m->pos_d[1];
        m->prev_pos_d[2] = m->pos_d[2];

        kepler_position(&m->orbit, ss->elapsed_time, m->pos_d);
        m->position = (HMM_Vec3){{
            (float)m->pos_d[0],
            (float)m->pos_d[1],
            (float)m->pos_d[2]
        }};
    }
}

int solar_system_find_gravity_body(const SolarSystem* ss,
                                   const double pos_d[3],
                                   float planet_radius) {
    /* Check distance to planet surface */
    double planet_dist = sqrt(pos_d[0]*pos_d[0] + pos_d[1]*pos_d[1] + pos_d[2]*pos_d[2]);
    double planet_alt = planet_dist - (double)planet_radius;

    /* Check each moon — find closest surface within SOI */
    double best_alt = planet_alt;
    int best_body = -1;  /* -1 = main planet */

    for (int i = 0; i < ss->moon_count; i++) {
        const CelestialBody* m = &ss->moons[i];
        /* For pinned body, use frozen center (Kepler pos_d drifts away from camera) */
        double mx, my, mz;
        if (i == ss->pinned_body) {
            mx = ss->pinned_center_d[0];
            my = ss->pinned_center_d[1];
            mz = ss->pinned_center_d[2];
        } else {
            mx = m->pos_d[0];
            my = m->pos_d[1];
            mz = m->pos_d[2];
        }
        double dx = pos_d[0] - mx;
        double dy = pos_d[1] - my;
        double dz = pos_d[2] - mz;
        double dist = sqrt(dx*dx + dy*dy + dz*dz);
        double alt = dist - (double)m->radius;

        if (alt < (double)MOON_SOI_RADIUS && alt < best_alt) {
            best_alt = alt;
            best_body = i;
        }
    }

    return best_body;
}

void solar_system_render(const SolarSystem* ss,
                         const Camera* cam,
                         HMM_Vec3 sun_dir,
                         HMM_Mat4 vp_rot,
                         sg_pipeline pip,
                         float Fcoef, float far_plane, float z_bias,
                         const double world_origin[3],
                         int lod_target_body,
                         sg_view planet_tex_view,
                         sg_sampler planet_tex_smp,
                         const void* planet_fallback_fs) {
    sg_apply_pipeline(pip);

    /* When on a moon (pinned_body >= 0), cam->pos_d is moon-local.
       Reconstruct world-space camera for rendering external bodies. */
    double cam_world[3];
    if (ss->pinned_body >= 0) {
        cam_world[0] = cam->pos_d[0] + ss->pinned_center_d[0];
        cam_world[1] = cam->pos_d[1] + ss->pinned_center_d[1];
        cam_world[2] = cam->pos_d[2] + ss->pinned_center_d[2];
    } else {
        cam_world[0] = cam->pos_d[0];
        cam_world[1] = cam->pos_d[1];
        cam_world[2] = cam->pos_d[2];
    }

    for (int i = 0; i < ss->moon_count; i++) {
        if (i == lod_target_body) continue;  /* LOD tree renders this one */
        const CelestialBody* m = &ss->moons[i];
        if (!m->mesh_ready) continue;

        /* Camera offset = cam_world - moon_pos (in double)
           Moon vertices are in moon-local coords (centered at 0).
           Shader does: cam_rel = vert_pos - camera_offset = vert_pos - cam + moon = world_vert - cam */
        double cam_off_d[3] = {
            cam_world[0] - m->pos_d[0],
            cam_world[1] - m->pos_d[1],
            cam_world[2] - m->pos_d[2],
        };
        float hi_x = (float)cam_off_d[0];
        float hi_y = (float)cam_off_d[1];
        float hi_z = (float)cam_off_d[2];
        float lo_x = (float)(cam_off_d[0] - (double)hi_x);
        float lo_y = (float)(cam_off_d[1] - (double)hi_y);
        float lo_z = (float)(cam_off_d[2] - (double)hi_z);

        planet_vs_params_t vs = {
            .mvp = vp_rot,
            .camera_offset = (HMM_Vec4){{hi_x, hi_y, hi_z, 0.0f}},
            .camera_offset_low = (HMM_Vec4){{lo_x, lo_y, lo_z, 0.0f}},
            .log_depth = (HMM_Vec4){{Fcoef, far_plane, z_bias, 0.0f}},
        };
        /* Pass camera_pos = cam - moon_center so shader reconstructs moon-local
           positions: world_pos_approx = cam_rel_pos + camera_pos = (vert-cam) + (cam-moon)
           = vert-moon. Then surface_dir = normalize(vert-moon) = radial from moon center.
           This gives correct terminator, AO, and rim lighting for moons. */
        planet_fs_params_t fs = {
            .sun_direction = (HMM_Vec4){{sun_dir.X, sun_dir.Y, sun_dir.Z, 0.3f}},
            .camera_pos = (HMM_Vec4){{hi_x, hi_y, hi_z, 0.0f}},  /* cam - moon (moon-local) */
            .atmos_params = (HMM_Vec4){{1.0f, 2.0f, 0.0f, 0.0f}},  /* pR=1, aR=2 (safe aThick=1), rScale=0 (no fog) */
            .lod_debug = (HMM_Vec4){{0.0f, 0.0f, 0.0f, 0.0f}},
            .dusk_sun_color = (HMM_Vec4){{1.0f, 0.85f, 0.7f, 0.0f}},
            .day_sun_color = (HMM_Vec4){{1.0f, 0.98f, 0.95f, 0.0f}},
            .planet_tex_params = (HMM_Vec4){{0.0f, 0.0f, 0.0f, 0.0f}},  /* disabled for moons */
        };

        sg_apply_uniforms(UB_planet_vs_params, &SG_RANGE(vs));
        sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(fs));

        sg_bindings bind = {0};
        bind.vertex_buffers[0] = m->gpu_buffer;
        bind.views[0] = planet_tex_view;
        bind.samplers[0] = planet_tex_smp;
        sg_apply_bindings(&bind);
        sg_draw(0, m->vertex_count, 1);
    }

    /* Draw Tenebris fallback mesh when LOD targets a moon */
    if (lod_target_body >= 0 && ss->planet_mesh_ready) {
        /* Tenebris mesh vertices are centered at origin (planet-local coords).
           Camera offset = cam_world - planet_center = cam_world - (0,0,0) = cam_world */
        double cam_off_d[3] = {
            cam_world[0],
            cam_world[1],
            cam_world[2],
        };
        float hi_x = (float)cam_off_d[0];
        float hi_y = (float)cam_off_d[1];
        float hi_z = (float)cam_off_d[2];
        float lo_x = (float)(cam_off_d[0] - (double)hi_x);
        float lo_y = (float)(cam_off_d[1] - (double)hi_y);
        float lo_z = (float)(cam_off_d[2] - (double)hi_z);

        planet_vs_params_t vs = {
            .mvp = vp_rot,
            .camera_offset = (HMM_Vec4){{hi_x, hi_y, hi_z, 0.0f}},
            .camera_offset_low = (HMM_Vec4){{lo_x, lo_y, lo_z, 0.0f}},
            .log_depth = (HMM_Vec4){{Fcoef, far_plane, z_bias, 0.0f}},
        };
        /* Use the same fs_params as the LOD path (passed from render.c) so that
           Tenebris looks identical whether viewed from space or from a moon's SOI.
           Only override camera_pos for the planet-centered reference frame. */
        planet_fs_params_t fs = *(const planet_fs_params_t*)planet_fallback_fs;
        fs.camera_pos = (HMM_Vec4){{hi_x, hi_y, hi_z, 0.0f}};

        sg_apply_uniforms(UB_planet_vs_params, &SG_RANGE(vs));
        sg_apply_uniforms(UB_planet_fs_params, &SG_RANGE(fs));

        sg_bindings bind = {0};
        bind.vertex_buffers[0] = ss->planet_buffer;
        bind.views[0] = planet_tex_view;
        bind.samplers[0] = planet_tex_smp;
        sg_apply_bindings(&bind);
        sg_draw(0, ss->planet_vertex_count, 1);
    }
}

void solar_system_draw_labels(const SolarSystem* ss,
                              const Camera* cam,
                              HMM_Mat4 vp_rot) {
    if (!cam->space_mode) return;

    float sw = sapp_widthf();
    float sh = sapp_heightf();

    /*
     * sdtx coordinate system:
     *   sdtx_canvas(sw*0.5, sh*0.5) was called, so the canvas is half-pixel sized.
     *   sdtx_origin(0.5, 0.5) was set.
     *   sdtx_pos(x, y) positions cursor in character-cell units (8 canvas pixels each).
     *   The origin offset is ADDED to pos, so pos(0,0) actually draws at cell (0.5, 0.5).
     *
     * To place text at a specific screen pixel:
     *   canvas_x = pixel_x * (canvas_w / screen_w) = pixel_x * 0.5
     *   cell_x   = canvas_x / 8.0  (each cell = 8 canvas pixels)
     *   sdtx_pos needs cell_x - origin_x, cell_y - origin_y
     */
    float origin_x = 0.5f;
    float origin_y = 0.5f;

    for (int i = 0; i < ss->moon_count; i++) {
        const CelestialBody* m = &ss->moons[i];

        /* Double-precision camera-relative position */
        double dx = m->pos_d[0] - cam->pos_d[0];
        double dy = m->pos_d[1] - cam->pos_d[1];
        double dz = m->pos_d[2] - cam->pos_d[2];
        float rel_x = (float)dx;
        float rel_y = (float)dy;
        float rel_z = (float)dz;

        /* Project through rotation-only VP (camera-relative input) */
        HMM_Vec4 clip = HMM_MulM4V4(vp_rot, (HMM_Vec4){{rel_x, rel_y, rel_z, 1.0f}});

        /* Behind camera? */
        if (clip.W <= 0.0f) continue;

        /* NDC [-1, 1] */
        float ndc_x = clip.X / clip.W;
        float ndc_y = clip.Y / clip.W;

        /* Off screen? (generous margin for labels near edge) */
        if (ndc_x < -1.2f || ndc_x > 1.2f || ndc_y < -1.2f || ndc_y > 1.2f) continue;

        /* NDC → pixel */
        float px = (ndc_x * 0.5f + 0.5f) * sw;
        float py = (-ndc_y * 0.5f + 0.5f) * sh;

        /* Pixel → sdtx cell (canvas = sw*0.5, cell = 8 canvas px) */
        float cell_x = (px * 0.5f) / 8.0f - origin_x;
        float cell_y = (py * 0.5f) / 8.0f - origin_y;

        /* Distance */
        double dist_d = sqrt(dx*dx + dy*dy + dz*dz);
        float dist_km = (float)(dist_d / 1000.0);

        /* Color: brighter when closer */
        float brightness = 1.0f - (dist_km / 20000.0f);
        if (brightness < 0.3f) brightness = 0.3f;
        if (brightness > 1.0f) brightness = 1.0f;

        /* Center name horizontally */
        int name_len = (int)strlen(m->name);
        float name_cell_x = cell_x - (float)name_len * 0.5f;

        /* Moon name */
        sdtx_color3f(brightness * 0.8f, brightness * 0.9f, brightness * 1.0f);
        sdtx_pos(name_cell_x, cell_y - 1.0f);
        sdtx_puts(m->name);

        /* Distance below name */
        sdtx_color3f(brightness * 0.5f, brightness * 0.6f, brightness * 0.7f);
        sdtx_pos(name_cell_x, cell_y);
        if (dist_km >= 1000.0f) {
            sdtx_printf("%.0fkm", dist_km);
        } else {
            sdtx_printf("%.1fkm", dist_km);
        }
    }

    /* ---- Tenebris (main planet at origin) ---- */
    {
        float rel_x = (float)(-cam->pos_d[0]);
        float rel_y = (float)(-cam->pos_d[1]);
        float rel_z = (float)(-cam->pos_d[2]);

        HMM_Vec4 clip = HMM_MulM4V4(vp_rot, (HMM_Vec4){{rel_x, rel_y, rel_z, 1.0f}});
        if (clip.W > 0.0f) {
            float ndc_x = clip.X / clip.W;
            float ndc_y = clip.Y / clip.W;
            if (ndc_x >= -1.2f && ndc_x <= 1.2f && ndc_y >= -1.2f && ndc_y <= 1.2f) {
                float px = (ndc_x * 0.5f + 0.5f) * sw;
                float py = (-ndc_y * 0.5f + 0.5f) * sh;
                float cell_x = (px * 0.5f) / 8.0f - origin_x;
                float cell_y = (py * 0.5f) / 8.0f - origin_y;

                double dist_d = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                                     cam->pos_d[1]*cam->pos_d[1] +
                                     cam->pos_d[2]*cam->pos_d[2]);
                float dist_km = (float)(dist_d / 1000.0);

                float brightness = 1.0f - (dist_km / 20000.0f);
                if (brightness < 0.3f) brightness = 0.3f;
                if (brightness > 1.0f) brightness = 1.0f;

                float name_cell_x = cell_x - 4.5f;  /* "Tenebris" = 8 chars, center */
                sdtx_color3f(brightness * 1.0f, brightness * 0.9f, brightness * 0.7f);
                sdtx_pos(name_cell_x, cell_y - 1.0f);
                sdtx_puts("Tenebris");

                sdtx_color3f(brightness * 0.6f, brightness * 0.5f, brightness * 0.4f);
                sdtx_pos(name_cell_x, cell_y);
                if (dist_km >= 1000.0f) {
                    sdtx_printf("%.0fkm", dist_km);
                } else {
                    sdtx_printf("%.1fkm", dist_km);
                }
            }
        }
    }
}

void solar_system_draw_orbits(const SolarSystem* ss,
                              const Camera* cam,
                              HMM_Mat4 vp_rot) {
    if (!cam->space_mode) return;

    /* Set up sgl with camera-relative perspective + depth testing */
    sgl_defaults();
    sgl_load_pipeline(ss->orbit_pip);
    sgl_matrix_mode_projection();
    sgl_load_matrix((const float*)&vp_rot);
    sgl_matrix_mode_modelview();
    sgl_load_identity();

    #define ORBIT_SEGMENTS 128

    for (int i = 0; i < ss->moon_count; i++) {
        const CelestialBody* m = &ss->moons[i];

        /* Dim orbit line color — slightly brighter for the pinned body */
        float alpha = (i == ss->pinned_body) ? 0.5f : 0.2f;

        sgl_begin_line_strip();
        sgl_c4f(0.4f, 0.6f, 0.9f, alpha);

        for (int s = 0; s <= ORBIT_SEGMENTS; s++) {
            /* Sample the orbit at evenly-spaced mean anomaly */
            double t_sample = (double)s / (double)ORBIT_SEGMENTS * (double)m->orbit.period;
            double pos[3];
            kepler_position(&m->orbit, t_sample, pos);

            /* Camera-relative position */
            float rx = (float)(pos[0] - cam->pos_d[0]);
            float ry = (float)(pos[1] - cam->pos_d[1]);
            float rz = (float)(pos[2] - cam->pos_d[2]);
            sgl_v3f(rx, ry, rz);
        }
        sgl_end();
    }

    sgl_draw();
    #undef ORBIT_SEGMENTS
}

void solar_system_shutdown(SolarSystem* ss) {
    for (int i = 0; i < ss->moon_count; i++) {
        if (ss->moons[i].mesh_ready) {
            sg_destroy_buffer(ss->moons[i].gpu_buffer);
            ss->moons[i].mesh_ready = false;
        }
    }
    if (ss->planet_mesh_ready) {
        sg_destroy_buffer(ss->planet_buffer);
        ss->planet_mesh_ready = false;
    }
}

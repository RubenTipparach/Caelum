#include "test_framework.h"
#include "planet.h"
#include "mesh.h"
#include "math_utils.h"
#include <stdlib.h>

// ---- Mesh generation ----

TEST(mesh_produces_vertices) {
    Planet p;
    planet_init(&p, 2);  // Small planet
    PlanetMesh mesh = planet_mesh_generate(&p);
    ASSERT(mesh.vertex_count > 0);
    ASSERT(mesh.vertices != NULL);
    printf("(%d verts) ", mesh.vertex_count);
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_vertex_count_divisible_by_3) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    // Each triangle = 3 vertices, so total must be divisible by 3
    ASSERT_EQ_INT(mesh.vertex_count % 3, 0);
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_no_nan_positions) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    for (int i = 0; i < mesh.vertex_count; i++) {
        ASSERT(!isnan(mesh.vertices[i].pos[0]));
        ASSERT(!isnan(mesh.vertices[i].pos[1]));
        ASSERT(!isnan(mesh.vertices[i].pos[2]));
    }
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_no_nan_normals) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    for (int i = 0; i < mesh.vertex_count; i++) {
        ASSERT(!isnan(mesh.vertices[i].normal[0]));
        ASSERT(!isnan(mesh.vertices[i].normal[1]));
        ASSERT(!isnan(mesh.vertices[i].normal[2]));
    }
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_normals_are_unit_length) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    // Spot check every 100th normal
    for (int i = 0; i < mesh.vertex_count; i += 100) {
        float* n = mesh.vertices[i].normal;
        float len = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        ASSERT_NEAR(len, 1.0f, 0.1f);  // Approximate unit length
    }
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_colors_in_range) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    for (int i = 0; i < mesh.vertex_count; i++) {
        for (int c = 0; c < 3; c++) {
            ASSERT(mesh.vertices[i].color[c] >= 0.0f);
            ASSERT(mesh.vertices[i].color[c] <= 1.0f);
        }
    }
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_no_degenerate_triangles) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    int degenerate = 0;
    for (int i = 0; i < mesh.vertex_count; i += 3) {
        // Check that the 3 positions aren't all identical
        float* a = mesh.vertices[i].pos;
        float* b = mesh.vertices[i+1].pos;
        float* c = mesh.vertices[i+2].pos;
        float ab2 = (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]);
        float bc2 = (b[0]-c[0])*(b[0]-c[0]) + (b[1]-c[1])*(b[1]-c[1]) + (b[2]-c[2])*(b[2]-c[2]);
        float ca2 = (c[0]-a[0])*(c[0]-a[0]) + (c[1]-a[1])*(c[1]-a[1]) + (c[2]-a[2])*(c[2]-a[2]);
        if (ab2 < 1e-12f || bc2 < 1e-12f || ca2 < 1e-12f) {
            degenerate++;
        }
    }
    if (degenerate > 0) {
        printf("(%d degenerate tris) ", degenerate);
    }
    // Allow a small number of degenerate triangles (< 1%)
    int tri_count = mesh.vertex_count / 3;
    ASSERT(degenerate < tri_count / 100 + 1);
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_vertices_near_planet_surface) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate(&p);
    float min_r = p.radius;
    float max_r = p.radius + MAX_VOXEL_HEIGHT * p.layer_thickness;
    for (int i = 0; i < mesh.vertex_count; i += 50) {  // Spot check
        float* pos = mesh.vertices[i].pos;
        float r = sqrtf(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        float eps = min_r * 1e-5f;  // Scale tolerance with radius
        ASSERT_GT(r, min_r - eps);
        ASSERT(r < max_r + eps);
    }
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_updates_after_break) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh1 = planet_mesh_generate(&p);
    int count1 = mesh1.vertex_count;
    planet_mesh_destroy(&mesh1);

    // Break several cells to ensure visible difference
    for (int i = 0; i < 5 && i < p.cell_count; i++) {
        planet_break_cell(&p, i);
    }

    PlanetMesh mesh2 = planet_mesh_generate(&p);
    int count2 = mesh2.vertex_count;
    planet_mesh_destroy(&mesh2);

    // After breaking cells, vertex count should change
    ASSERT(count1 != count2);
    planet_destroy(&p);
    PASS();
}

// ---- Nearby mesh tests ----

TEST(mesh_nearby_has_geometry) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate_nearby(&p, 0, 5);
    ASSERT(mesh.vertex_count > 0);
    ASSERT(mesh.vertices != NULL);
    printf("(%d verts) ", mesh.vertex_count);
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_nearby_produces_fewer_vertices) {
    Planet p;
    planet_init(&p, 4);  // 162 cells
    PlanetMesh full = planet_mesh_generate(&p);
    PlanetMesh nearby = planet_mesh_generate_nearby(&p, 0, 3);  // Only ~3 hops = small subset
    ASSERT(nearby.vertex_count > 0);
    ASSERT(nearby.vertex_count < full.vertex_count);
    printf("(nearby=%d < full=%d) ", nearby.vertex_count, full.vertex_count);
    planet_mesh_destroy(&full);
    planet_mesh_destroy(&nearby);
    planet_destroy(&p);
    PASS();
}

TEST(mesh_nearby_vertex_count_divisible_by_3) {
    Planet p;
    planet_init(&p, 2);
    PlanetMesh mesh = planet_mesh_generate_nearby(&p, 0, 5);
    ASSERT_EQ_INT(mesh.vertex_count % 3, 0);
    planet_mesh_destroy(&mesh);
    planet_destroy(&p);
    PASS();
}

void run_mesh_tests(void) {
    TEST_SUITE("Mesh Generation");
    RUN_TEST(mesh_produces_vertices);
    RUN_TEST(mesh_vertex_count_divisible_by_3);
    RUN_TEST(mesh_no_nan_positions);
    RUN_TEST(mesh_no_nan_normals);
    RUN_TEST(mesh_normals_are_unit_length);
    RUN_TEST(mesh_colors_in_range);
    RUN_TEST(mesh_no_degenerate_triangles);
    RUN_TEST(mesh_vertices_near_planet_surface);
    RUN_TEST(mesh_updates_after_break);

    TEST_SUITE("Nearby Mesh");
    RUN_TEST(mesh_nearby_has_geometry);
    RUN_TEST(mesh_nearby_produces_fewer_vertices);
    RUN_TEST(mesh_nearby_vertex_count_divisible_by_3);
}

#include "test_framework.h"
#include "planet.h"
#include "math_utils.h"
#include <stdlib.h>

// Use a small subdivision for fast tests
static Planet test_planet;
static int planet_initialized = 0;

static void ensure_planet(void) {
    if (!planet_initialized) {
        planet_init(&test_planet, 4);  // GP(4,0) = small, fast
        planet_initialized = 1;
    }
}

// ---- Planet generation ----

TEST(planet_cell_count_gp4) {
    ensure_planet();
    // GP(4,0): T = 4^2 = 16, cells = 10*T + 2 = 162
    ASSERT_EQ_INT(test_planet.cell_count, 162);
    PASS();
}

TEST(planet_has_12_pentagons) {
    ensure_planet();
    int pent_count = 0;
    for (int i = 0; i < test_planet.cell_count; i++) {
        if (test_planet.cells[i].vertex_count == 5) pent_count++;
    }
    ASSERT_EQ_INT(pent_count, 12);
    PASS();
}

TEST(planet_remaining_are_hexagons) {
    ensure_planet();
    int hex_count = 0;
    for (int i = 0; i < test_planet.cell_count; i++) {
        if (test_planet.cells[i].vertex_count == 6) hex_count++;
    }
    ASSERT_EQ_INT(hex_count, test_planet.cell_count - 12);
    PASS();
}

TEST(planet_cells_on_unit_sphere) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        HMM_Vec3 c = test_planet.cells[i].center;
        float len = sqrtf(c.X*c.X + c.Y*c.Y + c.Z*c.Z);
        ASSERT_NEAR(len, 1.0f, 1e-3f);
    }
    PASS();
}

TEST(planet_normals_match_centers) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        HMM_Vec3 n = test_planet.cells[i].normal;
        HMM_Vec3 c = vec3_normalize(test_planet.cells[i].center);
        ASSERT_NEAR(n.X, c.X, 1e-3f);
        ASSERT_NEAR(n.Y, c.Y, 1e-3f);
        ASSERT_NEAR(n.Z, c.Z, 1e-3f);
    }
    PASS();
}

// ---- Neighbor adjacency ----

TEST(planet_all_cells_have_neighbors) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        int expected = test_planet.cells[i].vertex_count;
        ASSERT_EQ_INT(test_planet.cells[i].neighbor_count, expected);
    }
    PASS();
}

TEST(planet_neighbors_are_valid_indices) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        for (int k = 0; k < test_planet.cells[i].neighbor_count; k++) {
            int nb = test_planet.cells[i].neighbors[k];
            ASSERT(nb >= 0);
            ASSERT(nb < test_planet.cell_count);
            ASSERT(nb != i);  // No self-neighbor
        }
    }
    PASS();
}

TEST(planet_neighbor_symmetry) {
    // If A is a neighbor of B, then B should be a neighbor of A
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        for (int k = 0; k < test_planet.cells[i].neighbor_count; k++) {
            int nb = test_planet.cells[i].neighbors[k];
            // Check that i appears in nb's neighbor list
            int found = 0;
            for (int j = 0; j < test_planet.cells[nb].neighbor_count; j++) {
                if (test_planet.cells[nb].neighbors[j] == i) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf("FAIL\n    Cell %d has neighbor %d, but %d does not list %d\n", i, nb, nb, i);
                _tests_failed++;
                return;
            }
        }
    }
    PASS();
}

// ---- Terrain generation ----

TEST(planet_sea_level_set) {
    ensure_planet();
    ASSERT(test_planet.sea_level > 0);
    ASSERT(test_planet.sea_level < MAX_VOXEL_HEIGHT);
    PASS();
}

TEST(planet_terrain_heights_valid) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        int h = test_planet.cells[i].terrain_height;
        ASSERT(h >= 0);
        ASSERT(h < MAX_VOXEL_HEIGHT);
    }
    PASS();
}

TEST(planet_voxel_columns_have_solid) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        int h = test_planet.cells[i].terrain_height;
        uint8_t vtype = test_planet.cells[i].voxels[h];
        // Terrain height should point to a solid (non-air) voxel
        ASSERT(vtype != VOXEL_AIR);
    }
    PASS();
}

TEST(planet_above_terrain_is_air) {
    ensure_planet();
    for (int i = 0; i < test_planet.cell_count; i++) {
        int h = test_planet.cells[i].terrain_height;
        // If terrain is above sea level, layer above terrain should be air
        if (h >= test_planet.sea_level && h + 1 < MAX_VOXEL_HEIGHT) {
            ASSERT_EQ_INT(test_planet.cells[i].voxels[h + 1], VOXEL_AIR);
        }
    }
    PASS();
}

// ---- Spatial queries ----

TEST(find_cell_north_pole) {
    ensure_planet();
    HMM_Vec3 north = {{0, 5.0f, 0}};  // Way above north pole
    int cell = planet_find_cell(&test_planet, north);
    ASSERT(cell >= 0);
    ASSERT(cell < test_planet.cell_count);
    // Cell normal should point roughly upward
    ASSERT_GT(test_planet.cells[cell].normal.Y, 0.5f);
    PASS();
}

TEST(find_cell_south_pole) {
    ensure_planet();
    HMM_Vec3 south = {{0, -5.0f, 0}};
    int cell = planet_find_cell(&test_planet, south);
    ASSERT(cell >= 0);
    // Cell normal should point roughly downward
    ASSERT_GT(-test_planet.cells[cell].normal.Y, 0.5f);
    PASS();
}

TEST(find_cell_returns_nearest) {
    ensure_planet();
    // A position right at cell 0's center should return cell 0
    HMM_Vec3 pos = vec3_scale(test_planet.cells[0].center, 2.0f);
    int cell = planet_find_cell(&test_planet, pos);
    ASSERT_EQ_INT(cell, 0);
    PASS();
}

TEST(terrain_radius_above_base) {
    ensure_planet();
    for (int i = 0; i < 10; i++) {  // Spot check first 10 cells
        float r = planet_terrain_radius(&test_planet, i);
        ASSERT_GT(r, test_planet.radius);
    }
    PASS();
}

TEST(terrain_radius_consistent) {
    ensure_planet();
    // terrain_radius should be radius + (top_layer+1) * layer_thickness
    int cell = 0;
    int top = test_planet.cells[cell].terrain_height;
    if (top < test_planet.sea_level) top = test_planet.sea_level;
    float expected = test_planet.radius + (top + 1) * test_planet.layer_thickness;
    float actual = planet_terrain_radius(&test_planet, cell);
    ASSERT_NEAR(actual, expected, 1e-5f);
    PASS();
}

// ---- Block interactions ----

TEST(break_cell_reduces_height) {
    // Use a fresh planet for mutation tests
    Planet p;
    planet_init(&p, 2);  // Tiny planet for speed
    int cell = 0;
    int old_height = p.cells[cell].terrain_height;
    int ok = planet_break_cell(&p, cell);
    ASSERT(ok);
    ASSERT(p.cells[cell].terrain_height < old_height);
    ASSERT(p.dirty);
    planet_destroy(&p);
    PASS();
}

TEST(place_cell_increases_height) {
    Planet p;
    planet_init(&p, 2);
    int cell = 0;
    int old_height = p.cells[cell].terrain_height;
    HMM_Vec3 dummy = {{0, 0, 0}};
    int result = planet_place_cell(&p, cell, dummy);
    ASSERT_EQ_INT(result, cell);
    ASSERT_EQ_INT(p.cells[cell].terrain_height, old_height + 1);
    ASSERT(p.dirty);
    planet_destroy(&p);
    PASS();
}

TEST(break_and_place_roundtrip) {
    Planet p;
    planet_init(&p, 2);
    int cell = 5;
    int original_height = p.cells[cell].terrain_height;
    planet_break_cell(&p, cell);
    HMM_Vec3 dummy = {{0, 0, 0}};
    planet_place_cell(&p, cell, dummy);
    // Height should be back to original (placed stone on top)
    ASSERT_EQ_INT(p.cells[cell].terrain_height, original_height);
    planet_destroy(&p);
    PASS();
}

// ---- Greedy walk + cell distance ----

TEST(find_cell_greedy_matches_brute_force) {
    Planet p;
    planet_init(&p, 4);
    // Test directions using greedy walk from a NEARBY hint (the actual use case).
    // The greedy result may land on an adjacent cell with nearly identical dot product
    // (Voronoi boundary). We accept this as correct if the dot product is within epsilon.
    HMM_Vec3 test_dirs[] = {
        {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}},
        {{-1, 0, 0}}, {{0, -1, 0}}, {{0, 0, -1}},
        {{1, 1, 0}}, {{1, 0, 1}}, {{0, 1, 1}},
        {{-1, 1, 0}}, {{1, -1, 0}}, {{1, 1, 1}},
        {{-1, -1, -1}}, {{0.3f, 0.7f, 0.1f}}, {{-0.5f, 0.2f, 0.8f}},
    };
    for (int i = 0; i < 15; i++) {
        HMM_Vec3 dir = vec3_normalize(test_dirs[i]);
        HMM_Vec3 pos = vec3_scale(dir, 150.0f);
        // Brute force to find the correct cell
        int brute = planet_find_cell_from_hint(&p, pos, -1);
        // Greedy from a neighbor of the correct cell
        int hint = p.cells[brute].neighbors[0];
        if (hint < 0) hint = brute;
        int greedy = planet_find_cell_from_hint(&p, pos, hint);
        // Accept if same cell, or if dot products are nearly identical (boundary case)
        float dot_brute = vec3_dot(dir, p.cells[brute].normal);
        float dot_greedy = vec3_dot(dir, p.cells[greedy].normal);
        ASSERT_NEAR(dot_greedy, dot_brute, 0.01f);
    }
    planet_destroy(&p);
    PASS();
}

TEST(cell_distance_self_zero) {
    ensure_planet();
    int d = planet_cell_distance(&test_planet, 0, 0, 100);
    ASSERT_EQ_INT(d, 0);
    PASS();
}

TEST(cell_distance_neighbor_one) {
    ensure_planet();
    int nb = test_planet.cells[0].neighbors[0];
    ASSERT(nb >= 0);
    int d = planet_cell_distance(&test_planet, 0, nb, 100);
    ASSERT_EQ_INT(d, 1);
    PASS();
}

void run_planet_tests(void) {
    TEST_SUITE("Planet Generation");
    RUN_TEST(planet_cell_count_gp4);
    RUN_TEST(planet_has_12_pentagons);
    RUN_TEST(planet_remaining_are_hexagons);
    RUN_TEST(planet_cells_on_unit_sphere);
    RUN_TEST(planet_normals_match_centers);

    TEST_SUITE("Neighbor Adjacency");
    RUN_TEST(planet_all_cells_have_neighbors);
    RUN_TEST(planet_neighbors_are_valid_indices);
    RUN_TEST(planet_neighbor_symmetry);

    TEST_SUITE("Terrain Generation");
    RUN_TEST(planet_sea_level_set);
    RUN_TEST(planet_terrain_heights_valid);
    RUN_TEST(planet_voxel_columns_have_solid);
    RUN_TEST(planet_above_terrain_is_air);

    TEST_SUITE("Spatial Queries");
    RUN_TEST(find_cell_north_pole);
    RUN_TEST(find_cell_south_pole);
    RUN_TEST(find_cell_returns_nearest);
    RUN_TEST(terrain_radius_above_base);
    RUN_TEST(terrain_radius_consistent);
    RUN_TEST(find_cell_greedy_matches_brute_force);
    RUN_TEST(cell_distance_self_zero);
    RUN_TEST(cell_distance_neighbor_one);

    TEST_SUITE("Block Interactions");
    RUN_TEST(break_cell_reduces_height);
    RUN_TEST(place_cell_increases_height);
    RUN_TEST(break_and_place_roundtrip);
}

#include "test_framework.h"
#include "HandmadeMath.h"
#include "math_utils.h"

// ---- vec3_normalize ----

TEST(normalize_unit_x) {
    HMM_Vec3 v = vec3_normalize((HMM_Vec3){{3.0f, 0.0f, 0.0f}});
    ASSERT_NEAR(v.X, 1.0f, 1e-5f);
    ASSERT_NEAR(v.Y, 0.0f, 1e-5f);
    ASSERT_NEAR(v.Z, 0.0f, 1e-5f);
    PASS();
}

TEST(normalize_diagonal) {
    HMM_Vec3 v = vec3_normalize((HMM_Vec3){{1.0f, 1.0f, 1.0f}});
    float expected = 1.0f / sqrtf(3.0f);
    ASSERT_NEAR(v.X, expected, 1e-5f);
    ASSERT_NEAR(v.Y, expected, 1e-5f);
    ASSERT_NEAR(v.Z, expected, 1e-5f);
    PASS();
}

// ---- vec3_dot ----

TEST(dot_perpendicular) {
    float d = vec3_dot((HMM_Vec3){{1, 0, 0}}, (HMM_Vec3){{0, 1, 0}});
    ASSERT_NEAR(d, 0.0f, 1e-5f);
    PASS();
}

TEST(dot_parallel) {
    float d = vec3_dot((HMM_Vec3){{2, 0, 0}}, (HMM_Vec3){{3, 0, 0}});
    ASSERT_NEAR(d, 6.0f, 1e-5f);
    PASS();
}

TEST(dot_antiparallel) {
    float d = vec3_dot((HMM_Vec3){{1, 0, 0}}, (HMM_Vec3){{-1, 0, 0}});
    ASSERT_NEAR(d, -1.0f, 1e-5f);
    PASS();
}

// ---- vec3_cross ----

TEST(cross_xy_gives_z) {
    HMM_Vec3 c = vec3_cross((HMM_Vec3){{1, 0, 0}}, (HMM_Vec3){{0, 1, 0}});
    ASSERT_NEAR(c.X, 0.0f, 1e-5f);
    ASSERT_NEAR(c.Y, 0.0f, 1e-5f);
    ASSERT_NEAR(c.Z, 1.0f, 1e-5f);
    PASS();
}

TEST(cross_anticommutative) {
    HMM_Vec3 a = {{1, 2, 3}};
    HMM_Vec3 b = {{4, 5, 6}};
    HMM_Vec3 ab = vec3_cross(a, b);
    HMM_Vec3 ba = vec3_cross(b, a);
    ASSERT_NEAR(ab.X, -ba.X, 1e-5f);
    ASSERT_NEAR(ab.Y, -ba.Y, 1e-5f);
    ASSERT_NEAR(ab.Z, -ba.Z, 1e-5f);
    PASS();
}

// ---- vec3_lerp ----

TEST(lerp_endpoints) {
    HMM_Vec3 a = {{0, 0, 0}};
    HMM_Vec3 b = {{10, 20, 30}};
    HMM_Vec3 r0 = vec3_lerp(a, b, 0.0f);
    HMM_Vec3 r1 = vec3_lerp(a, b, 1.0f);
    ASSERT_NEAR(r0.X, 0.0f, 1e-5f);
    ASSERT_NEAR(r1.X, 10.0f, 1e-5f);
    ASSERT_NEAR(r1.Y, 20.0f, 1e-5f);
    PASS();
}

TEST(lerp_midpoint) {
    HMM_Vec3 mid = vec3_lerp((HMM_Vec3){{0, 0, 0}}, (HMM_Vec3){{4, 8, 12}}, 0.5f);
    ASSERT_NEAR(mid.X, 2.0f, 1e-5f);
    ASSERT_NEAR(mid.Y, 4.0f, 1e-5f);
    ASSERT_NEAR(mid.Z, 6.0f, 1e-5f);
    PASS();
}

// ---- Icosahedron data ----

TEST(ico_vertex_count) {
    ASSERT_EQ_INT(ICO_VERTEX_COUNT, 12);
    PASS();
}

TEST(ico_face_count) {
    ASSERT_EQ_INT(ICO_FACE_COUNT, 20);
    PASS();
}

TEST(ico_vertices_on_sphere) {
    for (int i = 0; i < ICO_VERTEX_COUNT; i++) {
        HMM_Vec3 v = vec3_normalize(ICO_VERTICES[i]);
        float len = sqrtf(v.X*v.X + v.Y*v.Y + v.Z*v.Z);
        ASSERT_NEAR(len, 1.0f, 1e-4f);
    }
    PASS();
}

TEST(ico_faces_valid_indices) {
    for (int i = 0; i < ICO_FACE_COUNT; i++) {
        for (int j = 0; j < 3; j++) {
            ASSERT(ICO_FACES[i][j] >= 0);
            ASSERT(ICO_FACES[i][j] < ICO_VERTEX_COUNT);
        }
    }
    PASS();
}

void run_math_tests(void) {
    TEST_SUITE("Math Utilities");
    RUN_TEST(normalize_unit_x);
    RUN_TEST(normalize_diagonal);
    RUN_TEST(dot_perpendicular);
    RUN_TEST(dot_parallel);
    RUN_TEST(dot_antiparallel);
    RUN_TEST(cross_xy_gives_z);
    RUN_TEST(cross_anticommutative);
    RUN_TEST(lerp_endpoints);
    RUN_TEST(lerp_midpoint);
    RUN_TEST(ico_vertex_count);
    RUN_TEST(ico_face_count);
    RUN_TEST(ico_vertices_on_sphere);
    RUN_TEST(ico_faces_valid_indices);
}

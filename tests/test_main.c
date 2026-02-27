// Hex-Planets Test Suite
// Runs without GPU — tests planet generation, mesh, math, and physics logic.
// Build: cmake --build build --target hex-tests
// Run:   build\Release\hex-tests.exe  (or Debug)

#define HANDMADE_MATH_IMPLEMENTATION
#include "HandmadeMath.h"

#include "test_framework.h"
#include <stdio.h>

// Define the shared test counters
int _tests_run = 0;
int _tests_passed = 0;
int _tests_failed = 0;

// Test suite declarations
extern void run_math_tests(void);
extern void run_planet_tests(void);
extern void run_mesh_tests(void);

int main(void) {
    printf("========================================\n");
    printf("  Hex-Planets Test Suite\n");
    printf("========================================\n");

    run_math_tests();
    run_planet_tests();
    run_mesh_tests();

    TEST_SUMMARY();
    return TEST_EXIT_CODE();
}

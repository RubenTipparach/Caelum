#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <stdio.h>
#include <math.h>

// Shared counters (defined in test_main.c)
extern int _tests_run;
extern int _tests_passed;
extern int _tests_failed;

#define TEST(name) \
    static void test_##name(void); \
    static void run_test_##name(void) { \
        _tests_run++; \
        printf("  %-50s ", #name); \
        fflush(stdout); \
        test_##name(); \
    } \
    static void test_##name(void)

#define ASSERT(cond) do { \
    if (!(cond)) { \
        printf("FAIL\n    %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        _tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_EQ_INT(a, b) do { \
    int _a = (a), _b = (b); \
    if (_a != _b) { \
        printf("FAIL\n    %s:%d: %s == %d, expected %d\n", __FILE__, __LINE__, #a, _a, _b); \
        _tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, eps) do { \
    float _a = (a), _b = (b), _e = (eps); \
    if (fabsf(_a - _b) > _e) { \
        printf("FAIL\n    %s:%d: %s == %.6f, expected %.6f (eps=%.6f)\n", __FILE__, __LINE__, #a, _a, _b, _e); \
        _tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_GT(a, b) do { \
    float _a = (a), _b = (b); \
    if (!(_a > _b)) { \
        printf("FAIL\n    %s:%d: %s == %.6f, expected > %.6f\n", __FILE__, __LINE__, #a, _a, _b); \
        _tests_failed++; \
        return; \
    } \
} while(0)

#define PASS() do { _tests_passed++; printf("OK\n"); } while(0)

#define RUN_TEST(name) run_test_##name()

#define TEST_SUITE(name) printf("\n=== %s ===\n", name)

#define TEST_SUMMARY() do { \
    printf("\n========================================\n"); \
    printf("Results: %d/%d passed", _tests_passed, _tests_run); \
    if (_tests_failed > 0) printf(", %d FAILED", _tests_failed); \
    printf("\n========================================\n"); \
} while(0)

#define TEST_EXIT_CODE() (_tests_failed > 0 ? 1 : 0)

#endif

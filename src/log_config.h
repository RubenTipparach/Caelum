#ifndef LOG_CONFIG_H
#define LOG_CONFIG_H

#include <stdio.h>
#include <stdbool.h>

// ---- Log category enables (compile-time, flip to 0 to silence a category) ----
#define LOG_ENABLE_GAME       1
#define LOG_ENABLE_PLAYER     1
#define LOG_ENABLE_LOD        1
#define LOG_ENABLE_RENDER     1
#define LOG_ENABLE_ATMOS      1
#define LOG_ENABLE_JOBS       1
#define LOG_ENABLE_SCREENSHOT 1

// ---- Verbose logging (runtime toggle, press V in-game) ----
// When false, periodic high-frequency logs (position spam, LOD stats table) are suppressed.
// One-time and event-driven logs still print if their category is enabled.
extern bool log_verbose;

// ---- LOG macro ----
// Category-gated printf with auto [TAG] prefix.
// Usage: LOG(PLAYER, "Jetpack %s\n", "ON");
#define LOG(category, fmt, ...) \
    do { \
        if (LOG_ENABLE_##category) { \
            printf("[" #category "] " fmt, ##__VA_ARGS__); \
            fflush(stdout); \
        } \
    } while (0)

// ---- LOG_VERBOSE macro ----
// Same as LOG but also checks the runtime verbose flag.
// Use for periodic/spammy logs (every-N-frames position dumps, stats tables).
#define LOG_VERBOSE(category, fmt, ...) \
    do { \
        if (LOG_ENABLE_##category && log_verbose) { \
            printf("[" #category "] " fmt, ##__VA_ARGS__); \
            fflush(stdout); \
        } \
    } while (0)

// ---- RAW_LOG_VERBOSE ----
// Verbose-gated, no [TAG] prefix. For continuation lines (table rows, etc).
#define RAW_LOG_VERBOSE(category, fmt, ...) \
    do { \
        if (LOG_ENABLE_##category && log_verbose) { \
            printf(fmt, ##__VA_ARGS__); \
        } \
    } while (0)

#endif

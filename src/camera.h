#ifndef CAMERA_H
#define CAMERA_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "sokol_app.h"
#include "planet.h"
#include "lod.h"
#include "hex_terrain.h"

/* Forward declaration */
typedef struct SolarSystem SolarSystem;

typedef struct Camera {
    HMM_Vec3 position;   // float copy for rendering (derived from pos_d each frame)
    double pos_d[3];     // double-precision position accumulator (sub-mm at any distance)
    float yaw;           // horizontal angle (radians) in local tangent frame
    float pitch;         // vertical angle (radians), clamped
    float speed;         // movement speed (units/sec)
    float sensitivity;   // mouse look sensitivity

    HMM_Vec3 forward;   // look direction (in world space)
    HMM_Vec3 right;     // right direction (in world space)
    HMM_Vec3 up;        // up direction (in world space, = local_up)
    HMM_Vec3 local_up;  // normalize(position) — radial "up" on sphere
    HMM_Vec3 prev_local_up;  // Previous frame's local_up (for parallel transport)
    HMM_Vec3 tangent_north;  // Parallel-transported tangent frame north direction
    bool tangent_initialized;

    HMM_Vec3 velocity;
    float eye_height;    // height of eyes above ground
    float gravity;       // gravitational acceleration toward center
    bool on_ground;

    HMM_Mat4 view;
    HMM_Mat4 proj;
    bool mouse_locked;

    // Jetpack
    bool jetpack_active;
    float jetpack_speed_mult;  // scroll-wheel speed multiplier (floor=1.0)
    uint64_t last_space_time;

    // Movement keys held
    bool key_w, key_s, key_a, key_d, key_space, key_shift, key_ctrl;

    // Space flight mode
    bool space_mode;              // true when >50km from any surface
    float roll;                   // accumulated roll angle (Q/E)
    HMM_Vec3 space_up;           // free-form up vector in space mode
    HMM_Vec3 space_forward;      // free-form forward vector in space mode
    float space_prev_yaw;        // yaw at last frame (for computing mouse deltas)
    float space_prev_pitch;      // pitch at last frame (for computing mouse deltas)
    bool key_q, key_e;           // roll keys
    float mouse_dx_accum;        // accumulated mouse delta (applied once per frame)
    float mouse_dy_accum;

    // Mouse diagnostics (toggled with Alt+M)
    bool mouse_diag_enabled;
    int mouse_events_this_frame;  // count of MOUSE_MOVE events since last camera_update
    uint64_t mouse_last_event_time; // timestamp of last mouse event
    float mouse_max_gap_ms;       // largest gap between consecutive events this frame
    float mouse_max_delta;        // largest single-event delta this frame
    // Rolling stats for periodic log
    int diag_frame_count;
    int diag_total_events;
    int diag_zero_event_frames;   // frames with no mouse events while locked
    float diag_max_gap_ms;
    float diag_max_delta;
    float diag_max_accum;         // largest accumulated delta applied in one frame
    float diag_max_frame_ms;

    int gravity_body;             // -1 = Tenebris, 0-9 = moon index
    double gravity_center_d[3];   // Fixed center of gravity body (latched on transition)
    float transition_alpha;       // 0.0 = space, 1.0 = grounded (lerps over 0.5s)
    float tenebris_gravity;       // Tenebris surface gravity (from config, default 10.0)
} Camera;

void camera_init(Camera* cam);
void camera_update(Camera* cam, Planet* planet, const LodTree* lod,
                   const HexTerrain* hex, float dt,
                   const SolarSystem* solar);
void camera_handle_event(Camera* cam, const sapp_event* ev);

#endif

#ifndef CAMERA_H
#define CAMERA_H

#include <stdbool.h>
#include <stdint.h>
#include "HandmadeMath.h"
#include "sokol_app.h"
#include "planet.h"
#include "lod.h"

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
    bool key_w, key_s, key_a, key_d, key_space, key_shift;
} Camera;

void camera_init(Camera* cam);
void camera_update(Camera* cam, Planet* planet, const LodTree* lod, float dt);
void camera_handle_event(Camera* cam, const sapp_event* ev);

#endif

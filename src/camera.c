#include "camera.h"
#include "math_utils.h"
#include "sokol_time.h"
#include "log_config.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Physics constants (inspired by tenebris PlayerConfig)
#define WALK_SPEED      4.0f    // Meters per second on tangent plane
#define SPRINT_SPEED    8.0f
#define JUMP_FORCE      8.0f    // Initial upward velocity on jump (m/s)
#define BASE_GRAVITY    20.0f   // Gravitational acceleration toward center (m/s^2)
#define PLAYER_EYE_HEIGHT 1.6f  // Eye height above feet (meters)
#define PLAYER_HEIGHT   1.8f    // Total player height (meters)
#define PLAYER_RADIUS   0.3f    // Horizontal collision radius (meters)
#define GROUND_SNAP_THRESHOLD 0.5f  // Snap to ground if within this distance
#define AUTO_STEP_HEIGHT 0.3f   // Can step up blocks this tall without jumping

// Jetpack constants
#define JETPACK_THRUST      30.0f   // m/s^2 upward acceleration
#define JETPACK_MAX_SPEED   100.0f  // m/s max radial velocity
#define JETPACK_FLY_SPEED   20.0f   // m/s lateral movement while flying
#define JETPACK_DRAG        0.97f   // velocity damping per frame (air resistance)
#define DOUBLE_TAP_MS       300.0   // ms window for double-tap detection

// Logging throttle
static int log_frame_counter = 0;

// Estimate the approximate radius of a hex cell on the surface
static float estimate_cell_radius(const Planet* planet) {
    float surface_r = planet->radius + planet->sea_level * planet->layer_thickness;
    float cell_area = 4.0f * (float)M_PI * surface_r * surface_r / (float)planet->cell_count;
    // Area of regular hexagon = (3*sqrt(3)/2) * r^2, solve for r
    return sqrtf(cell_area / (3.0f * 1.732f / 2.0f));
}

// Check if a position collides with walls on the planet.
// Returns true if the position is blocked.
// current_ground_r: the ground radius at the player's starting cell (before movement)
static bool check_wall_collision(Planet* planet, HMM_Vec3 pos,
                                  float current_ground_r, float cell_radius) {
    int dest_cell = planet_find_cell(planet, pos);
    if (dest_cell < 0) return false;

    HMM_Vec3 up = HMM_NormV3(pos);
    float player_r = HMM_LenV3(pos);
    float feet_r = player_r - PLAYER_EYE_HEIGHT;
    float head_r = feet_r + PLAYER_HEIGHT;
    float step_threshold_r = current_ground_r + AUTO_STEP_HEIGHT;

    // Check destination cell and all its neighbors
    const HexCell* dcell = &planet->cells[dest_cell];
    int cells_to_check[7];
    int check_count = 0;
    cells_to_check[check_count++] = dest_cell;
    for (int i = 0; i < dcell->neighbor_count && check_count < 7; i++) {
        if (dcell->neighbors[i] >= 0) {
            cells_to_check[check_count++] = dcell->neighbors[i];
        }
    }

    for (int ci = 0; ci < check_count; ci++) {
        int cell_idx = cells_to_check[ci];
        const HexCell* cell = &planet->cells[cell_idx];

        // Horizontal distance: tangent-plane distance from player to cell center
        HMM_Vec3 cell_world = HMM_MulV3F(cell->center, player_r);
        HMM_Vec3 diff = HMM_SubV3(cell_world, pos);
        float radial_comp = HMM_DotV3(diff, up);
        HMM_Vec3 tangent_diff = HMM_SubV3(diff, HMM_MulV3F(up, radial_comp));
        float horiz_dist = HMM_LenV3(tangent_diff);

        // If cell center is too far, skip
        if (horiz_dist > PLAYER_RADIUS + cell_radius * 1.5f) continue;

        // Check each solid voxel layer in this cell
        for (int layer = 0; layer < MAX_VOXEL_HEIGHT; layer++) {
            if (cell->voxels[layer] == VOXEL_AIR) continue;

            float block_bot_r = planet->radius + layer * planet->layer_thickness;
            float block_top_r = block_bot_r + planet->layer_thickness;

            // Skip if this block is ground (below step threshold)
            if (block_top_r <= step_threshold_r) continue;

            // Skip if no vertical overlap with player body
            if (block_bot_r >= head_r || block_top_r <= feet_r) continue;

            // This block overlaps the player vertically and is above step height
            // Check horizontal proximity
            if (horiz_dist < PLAYER_RADIUS + cell_radius * 0.5f) {
                return true;  // Blocked!
            }
        }
    }

    return false;
}

void camera_init(Camera* cam) {
    cam->position = (HMM_Vec3){{0.0f, 796050.0f, 0.0f}};
    cam->yaw = 0.0f;
    cam->pitch = 0.0f;
    cam->speed = WALK_SPEED;
    cam->sensitivity = 0.003f;
    cam->mouse_locked = false;
    cam->key_w = cam->key_s = cam->key_a = cam->key_d = false;
    cam->key_space = cam->key_shift = false;

    cam->local_up = (HMM_Vec3){{0, 1, 0}};
    cam->forward = (HMM_Vec3){{0, 0, -1}};
    cam->right = (HMM_Vec3){{1, 0, 0}};
    cam->up = (HMM_Vec3){{0, 1, 0}};

    cam->velocity = (HMM_Vec3){{0, 0, 0}};
    cam->eye_height = PLAYER_EYE_HEIGHT;
    cam->gravity = BASE_GRAVITY;
    cam->on_ground = false;
    cam->jetpack_active = false;
    cam->last_space_time = 0;

    cam->view = HMM_M4D(1.0f);
    cam->proj = HMM_M4D(1.0f);
}

void camera_update(Camera* cam, Planet* planet, float dt) {
    // Local "up" is the radial direction from planet center
    float pos_len = HMM_LenV3(cam->position);
    if (pos_len < 0.001f) pos_len = 0.001f;
    cam->local_up = HMM_NormV3(cam->position);

    // Build a local reference frame on the sphere's tangent plane.
    HMM_Vec3 world_ref = (HMM_Vec3){{0, 0, 1}};
    if (fabsf(HMM_DotV3(cam->local_up, world_ref)) > 0.99f) {
        world_ref = (HMM_Vec3){{1, 0, 0}};
    }

    // tangent_north = component of world_ref perpendicular to local_up
    HMM_Vec3 tangent_north = HMM_SubV3(world_ref, HMM_MulV3F(cam->local_up, HMM_DotV3(world_ref, cam->local_up)));
    tangent_north = HMM_NormV3(tangent_north);

    // tangent_east = local_up cross tangent_north
    HMM_Vec3 tangent_east = HMM_Cross(cam->local_up, tangent_north);

    // Apply yaw to get forward direction on tangent plane
    HMM_Vec3 flat_forward = HMM_AddV3(
        HMM_MulV3F(tangent_north, cosf(cam->yaw)),
        HMM_MulV3F(tangent_east, sinf(cam->yaw))
    );
    flat_forward = HMM_NormV3(flat_forward);

    // Right direction on tangent plane
    HMM_Vec3 flat_right = HMM_Cross(cam->local_up, flat_forward);
    flat_right = HMM_NormV3(flat_right);

    // Full forward with pitch (for look direction only, not movement)
    cam->forward = HMM_AddV3(
        HMM_MulV3F(flat_forward, cosf(cam->pitch)),
        HMM_MulV3F(cam->local_up, sinf(cam->pitch))
    );
    cam->forward = HMM_NormV3(cam->forward);
    cam->right = flat_right;
    cam->up = cam->local_up;

    // ---- Get current ground info for collision ----
    float current_ground_r = planet ? planet->radius : pos_len;
    float cell_radius = 1.0f;
    if (planet) {
        int standing_cell = planet_find_cell(planet, cam->position);
        if (standing_cell >= 0) {
            current_ground_r = planet_terrain_radius(planet, standing_cell);
        }
        cell_radius = estimate_cell_radius(planet);
    }

    // ---- Movement with wall collision (axis-decomposed) ----
    float move_speed;
    if (cam->jetpack_active) {
        move_speed = JETPACK_FLY_SPEED;
    } else {
        move_speed = cam->key_shift ? SPRINT_SPEED : WALK_SPEED;
    }

    HMM_Vec3 forward_input = (HMM_Vec3){{0, 0, 0}};
    HMM_Vec3 side_input = (HMM_Vec3){{0, 0, 0}};
    if (cam->key_w) forward_input = HMM_AddV3(forward_input, flat_forward);
    if (cam->key_s) forward_input = HMM_SubV3(forward_input, flat_forward);
    if (cam->key_a) side_input = HMM_AddV3(side_input, flat_right);
    if (cam->key_d) side_input = HMM_SubV3(side_input, flat_right);

    // Normalize each component independently
    float fwd_len = HMM_LenV3(forward_input);
    float side_len = HMM_LenV3(side_input);
    if (fwd_len > 0.001f) forward_input = HMM_NormV3(forward_input);
    if (side_len > 0.001f) side_input = HMM_NormV3(side_input);

    HMM_Vec3 full_input = HMM_AddV3(forward_input, side_input);
    float full_len = HMM_LenV3(full_input);
    bool has_input = full_len > 0.001f;

    if (has_input && planet && !cam->jetpack_active) {
        // Ground movement with wall collision
        HMM_Vec3 full_dir = HMM_NormV3(full_input);
        HMM_Vec3 full_move = HMM_MulV3F(full_dir, move_speed * dt);
        HMM_Vec3 new_pos = HMM_AddV3(cam->position, full_move);

        // Try full movement
        if (!check_wall_collision(planet, new_pos, current_ground_r, cell_radius)) {
            cam->position = new_pos;
        } else {
            // Blocked: try forward only
            bool moved = false;
            if (fwd_len > 0.001f) {
                HMM_Vec3 fwd_move = HMM_MulV3F(forward_input, move_speed * dt);
                HMM_Vec3 fwd_pos = HMM_AddV3(cam->position, fwd_move);
                if (!check_wall_collision(planet, fwd_pos, current_ground_r, cell_radius)) {
                    cam->position = fwd_pos;
                    moved = true;
                }
            }
            // Try sideways only
            if (side_len > 0.001f) {
                HMM_Vec3 side_move = HMM_MulV3F(side_input, move_speed * dt);
                HMM_Vec3 side_pos = HMM_AddV3(cam->position, side_move);
                if (!check_wall_collision(planet, side_pos, current_ground_r, cell_radius)) {
                    cam->position = side_pos;
                    moved = true;
                }
            }
            (void)moved;
        }
    } else if (has_input) {
        // Jetpack flight or no planet — move freely (no wall collision)
        HMM_Vec3 full_dir = HMM_NormV3(full_input);
        cam->position = HMM_AddV3(cam->position, HMM_MulV3F(full_dir, move_speed * dt));
    }

    // ---- Jetpack / Jump / Gravity ----
    if (cam->jetpack_active) {
        if (cam->key_space) {
            // Thrust upward
            float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
            if (radial_vel < JETPACK_MAX_SPEED) {
                cam->velocity = HMM_AddV3(cam->velocity,
                    HMM_MulV3F(cam->local_up, JETPACK_THRUST * dt));
            }
        } else {
            // No thrust: reduced gravity for slow descent
            cam->velocity = HMM_SubV3(cam->velocity,
                HMM_MulV3F(cam->local_up, cam->gravity * 0.3f * dt));
        }
        // Air drag
        cam->velocity = HMM_MulV3F(cam->velocity, JETPACK_DRAG);
        cam->on_ground = false;
    } else {
        // Normal jump
        if (cam->key_space && cam->on_ground) {
            cam->velocity = HMM_MulV3F(cam->local_up, JUMP_FORCE);
            cam->on_ground = false;
        }
        // Normal gravity when airborne
        if (!cam->on_ground) {
            cam->velocity = HMM_SubV3(cam->velocity, HMM_MulV3F(cam->local_up, cam->gravity * dt));
        }
    }

    // ---- Apply velocity ----
    if (HMM_LenV3(cam->velocity) > 0.001f) {
        HMM_Vec3 vel_move = HMM_MulV3F(cam->velocity, dt);
        HMM_Vec3 new_pos = HMM_AddV3(cam->position, vel_move);

        if (planet && !cam->jetpack_active &&
            check_wall_collision(planet, new_pos, current_ground_r, cell_radius)) {
            // Blocked by wall while airborne — zero horizontal velocity, keep vertical
            float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
            cam->velocity = HMM_MulV3F(cam->local_up, radial_vel);
            vel_move = HMM_MulV3F(cam->velocity, dt);
            cam->position = HMM_AddV3(cam->position, vel_move);
        } else {
            cam->position = new_pos;
        }
    }

    // ---- Ground collision ----
    if (planet) {
        int cell = planet_find_cell(planet, cam->position);
        if (cell >= 0) {
            float ground_r = planet_terrain_radius(planet, cell);
            float player_r = HMM_LenV3(cam->position);
            float feet_r = player_r - cam->eye_height;
            float radial_vel = HMM_DotV3(cam->velocity, HMM_NormV3(cam->position));

            if (feet_r <= ground_r + GROUND_SNAP_THRESHOLD) {
                if (radial_vel <= 0.0f) {
                    // Landing
                    float target_r = ground_r + cam->eye_height;
                    cam->position = HMM_MulV3F(HMM_NormV3(cam->position), target_r);
                    cam->velocity = (HMM_Vec3){{0, 0, 0}};
                    cam->on_ground = true;
                    // Auto-deactivate jetpack on landing
                    if (cam->jetpack_active) {
                        cam->jetpack_active = false;
                        LOG(PLAYER, "Jetpack OFF (landed)\n");
                    }
                }
            } else {
                cam->on_ground = false;
            }
        }
    }

    // ---- View matrix ----
    HMM_Vec3 target = HMM_AddV3(cam->position, cam->forward);
    cam->view = HMM_LookAt_RH(cam->position, target, cam->local_up);

    // ---- Projection ----
    float aspect = sapp_widthf() / sapp_heightf();
    cam->proj = HMM_Perspective_RH_NO(HMM_ToRad(70.0f), aspect, 0.1f, 10000000.0f);

    // ---- Periodic logging ----
    log_frame_counter++;
    if (log_frame_counter % 120 == 0) {
        if (LOG_ENABLE_PLAYER && log_verbose) {
            printf("[PLAYER] pos=(%.3f, %.3f, %.3f) r=%.3f vel=(%.3f, %.3f, %.3f) ground=%d jetpack=%d",
                cam->position.X, cam->position.Y, cam->position.Z,
                pos_len,
                cam->velocity.X, cam->velocity.Y, cam->velocity.Z,
                cam->on_ground, cam->jetpack_active);
            if (planet) {
                int cell = planet_find_cell(planet, cam->position);
                if (cell >= 0) {
                    float gr = planet_terrain_radius(planet, cell);
                    printf(" cell=%d ground_r=%.3f feet_r=%.3f",
                        cell, gr, pos_len - cam->eye_height);
                }
            }
            printf("\n");
            fflush(stdout);
        }
    }
}

void camera_handle_event(Camera* cam, const sapp_event* ev) {
    switch (ev->type) {
        case SAPP_EVENTTYPE_MOUSE_DOWN:
            if (!cam->mouse_locked &&
                (ev->mouse_button == SAPP_MOUSEBUTTON_LEFT ||
                 ev->mouse_button == SAPP_MOUSEBUTTON_MIDDLE)) {
                cam->mouse_locked = true;
                sapp_lock_mouse(true);
            }
            break;

        case SAPP_EVENTTYPE_MOUSE_MOVE:
            if (cam->mouse_locked) {
                cam->yaw -= ev->mouse_dx * cam->sensitivity;
                cam->pitch -= ev->mouse_dy * cam->sensitivity;
                float limit = (float)M_PI / 2.0f - 0.01f;
                if (cam->pitch > limit) cam->pitch = limit;
                if (cam->pitch < -limit) cam->pitch = -limit;
            }
            break;

        case SAPP_EVENTTYPE_MOUSE_SCROLL:
            break;

        case SAPP_EVENTTYPE_KEY_DOWN:
            if (ev->key_code == SAPP_KEYCODE_W) cam->key_w = true;
            if (ev->key_code == SAPP_KEYCODE_S) cam->key_s = true;
            if (ev->key_code == SAPP_KEYCODE_A) cam->key_a = true;
            if (ev->key_code == SAPP_KEYCODE_D) cam->key_d = true;
            if (ev->key_code == SAPP_KEYCODE_SPACE) {
                if (!ev->key_repeat) {
                    // Double-tap detection for jetpack (ignore key repeats)
                    uint64_t now = stm_now();
                    double elapsed_ms = stm_ms(stm_diff(now, cam->last_space_time));
                    if (cam->last_space_time != 0 && elapsed_ms < DOUBLE_TAP_MS) {
                        cam->jetpack_active = !cam->jetpack_active;
                        LOG(PLAYER, "Jetpack %s\n", cam->jetpack_active ? "ON" : "OFF");
                        cam->last_space_time = 0;  // Reset so held key can't re-toggle
                    } else {
                        cam->last_space_time = now;
                    }
                }
                cam->key_space = true;
            }
            if (ev->key_code == SAPP_KEYCODE_LEFT_SHIFT) cam->key_shift = true;
            if (ev->key_code == SAPP_KEYCODE_ESCAPE) {
                cam->mouse_locked = false;
                sapp_lock_mouse(false);
            }
            break;

        case SAPP_EVENTTYPE_KEY_UP:
            if (ev->key_code == SAPP_KEYCODE_W) cam->key_w = false;
            if (ev->key_code == SAPP_KEYCODE_S) cam->key_s = false;
            if (ev->key_code == SAPP_KEYCODE_A) cam->key_a = false;
            if (ev->key_code == SAPP_KEYCODE_D) cam->key_d = false;
            if (ev->key_code == SAPP_KEYCODE_SPACE) cam->key_space = false;
            if (ev->key_code == SAPP_KEYCODE_LEFT_SHIFT) cam->key_shift = false;
            break;

        default:
            break;
    }
}

#include "camera.h"
#include "celestial.h"
#include "math_utils.h"
#include "sokol_time.h"
#include "log_config.h"
#include <math.h>

#ifdef __APPLE__
extern void platform_macos_set_mouse_tap(bool active);
extern void platform_macos_poll_mouse(float* out_dx, float* out_dy);
extern bool platform_macos_has_tap(void);
extern void platform_macos_poll_events(void);
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Physics constants (inspired by tenebris PlayerConfig)
#define WALK_SPEED      4.0f    // Meters per second on tangent plane
#define SPRINT_SPEED    8.0f
#define JUMP_FORCE      8.0f    // Initial upward velocity on jump (m/s)
#define DEFAULT_GRAVITY 20.0f   // Fallback gravity if tenebris_gravity not set
#define PLAYER_EYE_HEIGHT 1.6f  // Eye height above feet (meters)
#define PLAYER_HEIGHT   1.8f    // Total player height (1.8m, fits in 2m caves)
#define PLAYER_RADIUS   0.3f    // Horizontal collision radius (meters)
#define GROUND_SNAP_THRESHOLD 0.5f  // Snap to ground if within this distance
#define AUTO_STEP_HEIGHT 1.0f   // Can step up blocks this tall without jumping (1 hex layer)

// Jetpack constants
#define JETPACK_THRUST      30.0f   // m/s^2 upward acceleration
#define JETPACK_MAX_SPEED   100.0f  // m/s max radial velocity
#define JETPACK_FLY_SPEED   20.0f   // m/s lateral movement while flying
#define JETPACK_DRAG        0.97f   // velocity damping per frame (air resistance)
#define DOUBLE_TAP_MS       300.0   // ms window for double-tap detection

// Rodrigues rotation: rotate vector v around unit axis k by angle theta
static HMM_Vec3 rodrigues(HMM_Vec3 v, HMM_Vec3 k, float theta) {
    float c = cosf(theta);
    float s = sinf(theta);
    float dot = HMM_DotV3(k, v);
    HMM_Vec3 cross = HMM_Cross(k, v);
    return HMM_AddV3(
        HMM_AddV3(HMM_MulV3F(v, c), HMM_MulV3F(cross, s)),
        HMM_MulV3F(k, dot * (1.0f - c))
    );
}

// Compute the column center surface point on the ellipsoid for moon hex collision.
// col_r is the ellipsoid radius at the column center direction (from hex_terrain_col_base_r).
static HMM_Vec3 hex_col_surface(const HexTerrain* hex, int gcol, int grow, float col_r) {
    float lx = gcol * HEX_COL_SPACING;
    float lz = ((gcol & 1) ? (grow + 0.5f) : (float)grow) * HEX_ROW_SPACING;
    HMM_Vec3 col_dir = HMM_NormV3(HMM_AddV3(hex->tangent_origin,
        HMM_AddV3(HMM_MulV3F(hex->tangent_east, lx),
                   HMM_MulV3F(hex->tangent_north, lz))));
    return HMM_MulV3F(col_dir, col_r);
}

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
            if (cell->voxels[layer] == VOXEL_AIR || cell->voxels[layer] == VOXEL_TORCH) continue;

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
    cam->pos_d[0] = 0.0;
    cam->pos_d[1] = 796050.0;
    cam->pos_d[2] = 0.0;
    cam->yaw = 0.0f;
    cam->pitch = 0.0f;
    cam->speed = WALK_SPEED;
    cam->sensitivity = 0.003f;
    cam->mouse_locked = false;
    cam->key_w = cam->key_s = cam->key_a = cam->key_d = false;
    cam->key_space = cam->key_shift = cam->key_ctrl = false;

    cam->local_up = (HMM_Vec3){{0, 1, 0}};
    cam->prev_local_up = (HMM_Vec3){{0, 1, 0}};
    cam->tangent_north = (HMM_Vec3){{0, 0, -1}};
    cam->tangent_initialized = false;
    cam->forward = (HMM_Vec3){{0, 0, -1}};
    cam->right = (HMM_Vec3){{1, 0, 0}};
    cam->up = (HMM_Vec3){{0, 1, 0}};

    cam->velocity = (HMM_Vec3){{0, 0, 0}};
    cam->eye_height = PLAYER_EYE_HEIGHT;
    cam->gravity = DEFAULT_GRAVITY;
    cam->on_ground = false;
    cam->jetpack_active = false;
    cam->jetpack_speed_mult = 1.0f;
    cam->last_space_time = 0;

    cam->space_mode = false;
    cam->roll = 0.0f;
    cam->space_up = (HMM_Vec3){{0, 1, 0}};
    cam->space_forward = (HMM_Vec3){{0, 0, -1}};
    cam->space_prev_yaw = 0.0f;
    cam->space_prev_pitch = 0.0f;
    cam->key_q = cam->key_e = false;
    cam->gravity_body = -1;
    cam->gravity_center_d[0] = 0.0;
    cam->gravity_center_d[1] = 0.0;
    cam->gravity_center_d[2] = 0.0;
    cam->transition_alpha = 1.0f;

    cam->view = HMM_M4D(1.0f);
    cam->proj = HMM_M4D(1.0f);
}

void camera_update(Camera* cam, Planet* planet, const LodTree* lod,
                   const HexTerrain* hex, float dt,
                   const SolarSystem* solar) {
    double dd = (double)dt;

    // ---- Apply mouse deltas (once per frame) ----
    // On macOS: drain pending mouse events from the run loop to avoid
    // display-link starvation, then prefer CGEventTap if available.
    #ifdef __APPLE__
    platform_macos_poll_events();
    #endif

    float mouse_dx = cam->mouse_dx_accum;
    float mouse_dy = cam->mouse_dy_accum;
    cam->mouse_dx_accum = 0.0f;
    cam->mouse_dy_accum = 0.0f;

    #ifdef __APPLE__
    if (cam->mouse_locked && platform_macos_has_tap()) {
        float tap_dx, tap_dy;
        platform_macos_poll_mouse(&tap_dx, &tap_dy);
        // Use tap deltas if available (more reliable than sokol NSEvents)
        if (tap_dx != 0.0f || tap_dy != 0.0f) {
            mouse_dx = tap_dx;
            mouse_dy = tap_dy;
        }
        // If tap returned zero but sokol had deltas, use sokol deltas
    }
    #endif

    float applied_accum = 0.0f;
    if (mouse_dx != 0.0f || mouse_dy != 0.0f) {
        applied_accum = sqrtf(mouse_dx * mouse_dx + mouse_dy * mouse_dy);
        cam->yaw -= mouse_dx * cam->sensitivity;
        cam->pitch -= mouse_dy * cam->sensitivity;
        if (!cam->space_mode) {
            float limit = (float)M_PI / 2.0f - 0.01f;
            if (cam->pitch > limit) cam->pitch = limit;
            if (cam->pitch < -limit) cam->pitch = -limit;
        }
    }

    // ---- Mouse diagnostics (Alt+M to toggle) ----
    if (cam->mouse_diag_enabled) {
        float frame_ms = dt * 1000.0f;

        // Accumulate rolling stats
        cam->diag_frame_count++;
        cam->diag_total_events += cam->mouse_events_this_frame;
        if (cam->mouse_locked && cam->mouse_events_this_frame == 0)
            cam->diag_zero_event_frames++;
        if (cam->mouse_max_gap_ms > cam->diag_max_gap_ms)
            cam->diag_max_gap_ms = cam->mouse_max_gap_ms;
        if (cam->mouse_max_delta > cam->diag_max_delta)
            cam->diag_max_delta = cam->mouse_max_delta;
        if (applied_accum > cam->diag_max_accum)
            cam->diag_max_accum = applied_accum;
        if (frame_ms > cam->diag_max_frame_ms)
            cam->diag_max_frame_ms = frame_ms;

        // Log every 120 frames (~2 seconds at 60fps)
        if (cam->diag_frame_count >= 120) {
            float avg_events = (float)cam->diag_total_events / (float)cam->diag_frame_count;
            printf("[MOUSE DIAG] %d frames: avg_events/frame=%.1f, "
                   "zero_event_frames=%d, max_gap=%.1fms, "
                   "max_single_delta=%.0f, max_accum=%.0f, "
                   "max_frame=%.1fms\n",
                   cam->diag_frame_count, avg_events,
                   cam->diag_zero_event_frames,
                   cam->diag_max_gap_ms,
                   cam->diag_max_delta,
                   cam->diag_max_accum,
                   cam->diag_max_frame_ms);
            fflush(stdout);

            // Reset
            cam->diag_frame_count = 0;
            cam->diag_total_events = 0;
            cam->diag_zero_event_frames = 0;
            cam->diag_max_gap_ms = 0.0f;
            cam->diag_max_delta = 0.0f;
            cam->diag_max_accum = 0.0f;
            cam->diag_max_frame_ms = 0.0f;
        }

        // Reset per-frame counters
        cam->mouse_events_this_frame = 0;
        cam->mouse_max_gap_ms = 0.0f;
        cam->mouse_max_delta = 0.0f;
    }

    // ---- Teleport detection: capture radius before movement ----
    double radius_before = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                                cam->pos_d[1]*cam->pos_d[1] +
                                cam->pos_d[2]*cam->pos_d[2]);

    // ---- Sync float position from double (authoritative) ----
    cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

    // ---- Determine gravity body (planet vs moon) ----
    float planet_surface_r = 0.0f;
    if (planet) {
        planet_surface_r = planet->radius + planet->sea_level * planet->layer_thickness;
    }
    int prev_gravity = cam->gravity_body;

    // When on a moon, cam->pos_d is moon-local. Convert to world-space for SOI checks.
    double world_pos_d[3];
    if (cam->gravity_body >= 0 && solar && solar->pinned_body >= 0) {
        world_pos_d[0] = cam->pos_d[0] + solar->pinned_center_d[0];
        world_pos_d[1] = cam->pos_d[1] + solar->pinned_center_d[1];
        world_pos_d[2] = cam->pos_d[2] + solar->pinned_center_d[2];
    } else {
        world_pos_d[0] = cam->pos_d[0];
        world_pos_d[1] = cam->pos_d[1];
        world_pos_d[2] = cam->pos_d[2];
    }

    if (solar) {
        cam->gravity_body = solar_system_find_gravity_body(solar, world_pos_d, planet_surface_r);
    } else {
        cam->gravity_body = -1;
    }

    // ---- Gravity body center tracking ----
    if (prev_gravity != cam->gravity_body) {
        // Transition: gravity_center is always {0,0,0} — either planet origin
        // (world-space) or moon origin (moon-local frame, set up by main.c).
        cam->gravity_center_d[0] = 0.0;
        cam->gravity_center_d[1] = 0.0;
        cam->gravity_center_d[2] = 0.0;
        cam->tangent_initialized = false;
    }
    // No orbital delta needed — when on a moon, cam->pos_d is already moon-local
    // and the moon doesn't move in its own reference frame.

    // Set gravity based on current body
    if (cam->gravity_body >= 0 && solar) {
        cam->gravity = solar->moons[cam->gravity_body].surface_gravity;
    } else {
        cam->gravity = cam->tenebris_gravity > 0.0f ? cam->tenebris_gravity : DEFAULT_GRAVITY;
    }

    // ---- Compute body-relative distance and local_up ----
    // body_dist: distance from gravity body center (for altitude, collision, ground snap)
    // pos_len:   distance from world origin (for space mode planet check)
    double body_dx = cam->pos_d[0] - cam->gravity_center_d[0];
    double body_dy = cam->pos_d[1] - cam->gravity_center_d[1];
    double body_dz = cam->pos_d[2] - cam->gravity_center_d[2];
    double body_dist_d = sqrt(body_dx*body_dx + body_dy*body_dy + body_dz*body_dz);
    if (body_dist_d < 0.001) body_dist_d = 0.001;
    float body_dist = (float)body_dist_d;

    // Radial direction (sphere normal)
    HMM_Vec3 radial_up = {{
        (float)(body_dx / body_dist_d),
        (float)(body_dy / body_dist_d),
        (float)(body_dz / body_dist_d)
    }};

    // For moons: ellipsoid surface normal (like eNormal in hex-shell.html).
    // For planet: radial direction (sphere normal).
    if (cam->gravity_body >= 0 && solar) {
        cam->local_up = moon_ellipsoid_normal(
            &solar->moons[cam->gravity_body].shape, radial_up);
    } else {
        cam->local_up = radial_up;
    }

    // pos_len: distance from Tenebris center (world-space, for planet altitude check)
    double pos_len_d = sqrt(world_pos_d[0]*world_pos_d[0] +
                            world_pos_d[1]*world_pos_d[1] +
                            world_pos_d[2]*world_pos_d[2]);
    if (pos_len_d < 0.001) pos_len_d = 0.001;
    float pos_len = (float)pos_len_d;

    // ---- Space mode detection ----
    // Check if we're more than 50km from ALL surfaces (using world-space positions)
    {
        double planet_alt = pos_len_d - (double)planet_surface_r;
        bool within_any_soi = (planet_alt < (double)MOON_SOI_RADIUS);

        if (!within_any_soi && solar) {
            for (int i = 0; i < solar->moon_count; i++) {
                const CelestialBody* m = &solar->moons[i];
                double mx = m->pos_d[0];
                double my = m->pos_d[1];
                double mz = m->pos_d[2];
                double dx = world_pos_d[0] - mx;
                double dy = world_pos_d[1] - my;
                double dz = world_pos_d[2] - mz;
                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                if (dist - (double)m->radius < (double)MOON_SOI_RADIUS) {
                    within_any_soi = true;
                    break;
                }
            }
        }

        bool was_space = cam->space_mode;
        cam->space_mode = !within_any_soi;

        // Transition: entering space
        if (cam->space_mode && !was_space) {
            cam->space_up = cam->local_up;      // Snapshot current orientation
            cam->space_forward = cam->forward;   // Snapshot current look direction
            cam->space_prev_yaw = cam->yaw;      // Snapshot yaw/pitch for delta tracking
            cam->space_prev_pitch = cam->pitch;
            cam->transition_alpha = 0.0f;
            cam->roll = 0.0f;
        }
        // Transition: leaving space (entering SOI)
        if (!cam->space_mode && was_space) {
            cam->transition_alpha = 0.0f;  // Will lerp to 1.0
        }

        // Transition gravity body change (moon to moon or moon to planet)
        if (!cam->space_mode && prev_gravity != cam->gravity_body) {
            cam->transition_alpha = 0.0f;
        }
    }

    // ---- Update transition alpha ----
    if (cam->space_mode) {
        cam->transition_alpha = 0.0f;
    } else {
        if (cam->transition_alpha < 1.0f) {
            cam->transition_alpha += dt * 2.0f;  // 0.5s transition
            if (cam->transition_alpha > 1.0f) cam->transition_alpha = 1.0f;
            // Lerp space_up toward radial local_up during transition
            cam->space_up = HMM_NormV3(HMM_LerpV3(cam->space_up, cam->transition_alpha, cam->local_up));
        }
    }

    // ---- Orientation: space mode uses incremental rotations, grounded uses tangent frame ----
    HMM_Vec3 flat_forward, flat_right;

    if (cam->space_mode) {
        // === SPACE MODE: incremental rotation from mouse deltas ===
        // Compute how much yaw/pitch changed this frame
        float delta_yaw = cam->yaw - cam->space_prev_yaw;
        float delta_pitch = cam->pitch - cam->space_prev_pitch;
        cam->space_prev_yaw = cam->yaw;
        cam->space_prev_pitch = cam->pitch;

        // Current right vector (before this frame's rotations)
        HMM_Vec3 right = HMM_NormV3(HMM_Cross(cam->space_forward, cam->space_up));

        // Yaw: rotate forward around up
        if (fabsf(delta_yaw) > 1e-6f) {
            cam->space_forward = rodrigues(cam->space_forward, cam->space_up, delta_yaw);
        }
        // Pitch: rotate forward and up around right
        if (fabsf(delta_pitch) > 1e-6f) {
            cam->space_forward = rodrigues(cam->space_forward, right, delta_pitch);
            cam->space_up = rodrigues(cam->space_up, right, delta_pitch);
        }
        // Roll (Q/E): rotate up around forward
        float roll_amount = 0.0f;
        if (cam->key_q) roll_amount -= 2.0f * dt;
        if (cam->key_e) roll_amount += 2.0f * dt;
        if (fabsf(roll_amount) > 1e-6f) {
            cam->space_up = rodrigues(cam->space_up, cam->space_forward, roll_amount);
        }

        // Orthonormalize to prevent drift
        cam->space_forward = HMM_NormV3(cam->space_forward);
        cam->space_up = HMM_SubV3(cam->space_up,
            HMM_MulV3F(cam->space_forward, HMM_DotV3(cam->space_up, cam->space_forward)));
        cam->space_up = HMM_NormV3(cam->space_up);

        // Set camera vectors (cross order matches grounded: cross(up, forward) for movement)
        cam->forward = cam->space_forward;
        cam->up = cam->space_up;
        cam->local_up = cam->space_up;
        cam->right = HMM_NormV3(HMM_Cross(cam->space_up, cam->space_forward));
        cam->prev_local_up = cam->local_up;

        // For movement code below
        flat_forward = cam->forward;  // in space mode, movement follows look direction
        flat_right = cam->right;
    } else {
        // === GROUNDED MODE: tangent frame + absolute yaw/pitch ===

        // Build a local reference frame on the sphere's tangent plane
        // using parallel transport to avoid rotation snapping.
        if (!cam->tangent_initialized) {
            // First frame: initialize tangent from global reference
            HMM_Vec3 world_ref = (HMM_Vec3){{0, 0, 1}};
            if (fabsf(HMM_DotV3(cam->local_up, world_ref)) > 0.99f) {
                world_ref = (HMM_Vec3){{1, 0, 0}};
            }
            cam->tangent_north = HMM_SubV3(world_ref,
                HMM_MulV3F(cam->local_up, HMM_DotV3(world_ref, cam->local_up)));
            cam->tangent_north = HMM_NormV3(cam->tangent_north);
            cam->tangent_initialized = true;
        } else {
            // Parallel transport: rotate tangent_north from prev_up to current local_up
            HMM_Vec3 prev_up = cam->prev_local_up;
            float cos_a = HMM_DotV3(prev_up, cam->local_up);
            if (cos_a < 0.999999f) {
                HMM_Vec3 axis = HMM_NormV3(HMM_Cross(prev_up, cam->local_up));
                float angle = acosf(fminf(fmaxf(cos_a, -1.0f), 1.0f));
                float c = cosf(angle), s = sinf(angle);
                HMM_Vec3 tn = cam->tangent_north;
                float dot_ax = HMM_DotV3(axis, tn);
                HMM_Vec3 cross_ax = HMM_Cross(axis, tn);
                cam->tangent_north = HMM_AddV3(
                    HMM_AddV3(HMM_MulV3F(tn, c), HMM_MulV3F(cross_ax, s)),
                    HMM_MulV3F(axis, dot_ax * (1.0f - c))
                );
            }
            // Re-orthogonalize against new local_up (remove accumulated drift)
            cam->tangent_north = HMM_SubV3(cam->tangent_north,
                HMM_MulV3F(cam->local_up, HMM_DotV3(cam->tangent_north, cam->local_up)));
            float tn_len = HMM_LenV3(cam->tangent_north);
            if (tn_len < 0.001f) {
                HMM_Vec3 world_ref = (HMM_Vec3){{0, 0, 1}};
                if (fabsf(HMM_DotV3(cam->local_up, world_ref)) > 0.99f) {
                    world_ref = (HMM_Vec3){{1, 0, 0}};
                }
                cam->tangent_north = HMM_SubV3(world_ref,
                    HMM_MulV3F(cam->local_up, HMM_DotV3(world_ref, cam->local_up)));
            }
            cam->tangent_north = HMM_NormV3(cam->tangent_north);
        }
        cam->prev_local_up = cam->local_up;

        HMM_Vec3 tangent_north = cam->tangent_north;
        HMM_Vec3 tangent_east = HMM_Cross(cam->local_up, tangent_north);

        flat_forward = HMM_AddV3(
            HMM_MulV3F(tangent_north, cosf(cam->yaw)),
            HMM_MulV3F(tangent_east, sinf(cam->yaw))
        );
        flat_forward = HMM_NormV3(flat_forward);

        flat_right = HMM_Cross(cam->local_up, flat_forward);
        flat_right = HMM_NormV3(flat_right);

        cam->forward = HMM_AddV3(
            HMM_MulV3F(flat_forward, cosf(cam->pitch)),
            HMM_MulV3F(cam->local_up, sinf(cam->pitch))
        );
        cam->forward = HMM_NormV3(cam->forward);
        cam->right = flat_right;
        cam->up = cam->local_up;
    }

    // ---- Get current ground info for collision ----
    // current_ground_r is a body-relative radius (distance from gravity body center)
    float current_ground_r = body_dist;
    if (cam->gravity_body >= 0 && solar) {
        // On a moon: use the moon's shape model (NOT the Tenebris LOD tree)
        // Use radial direction (not ellipsoid normal) for surface radius lookup
        current_ground_r = moon_surface_radius(
            &solar->moons[cam->gravity_body].shape, radial_up);
    } else if (lod) {
        current_ground_r = (float)lod_tree_terrain_height(lod, cam->position);
    } else if (planet) {
        current_ground_r = planet->radius + planet->sea_level * planet->layer_thickness;
    }

    // ---- Cap jetpack speed by altitude ----
    // <1km: max 25x, 1-5km: max 200x, 5-50km: max 500x, >50km: unlimited
    // Once you drop below 50km, speed is clamped to the tier cap.
    {
        float altitude = body_dist - current_ground_r;
        float max_mult;
        if (altitude >= 50000.0f) {
            max_mult = 1e9f;  // unlimited above 50km
        } else if (altitude >= 5000.0f) {
            max_mult = 500.0f;
        } else if (altitude >= 1000.0f) {
            max_mult = 200.0f;
        } else {
            max_mult = 25.0f;
        }
        if (cam->jetpack_speed_mult > max_mult) {
            cam->jetpack_speed_mult = max_mult;
        }
    }

    // ---- Movement (accumulated in double precision) ----
    float move_speed;
    if (cam->space_mode) {
        move_speed = JETPACK_FLY_SPEED * cam->jetpack_speed_mult;
    } else if (cam->jetpack_active) {
        move_speed = JETPACK_FLY_SPEED * cam->jetpack_speed_mult;
    } else {
        move_speed = cam->key_shift ? SPRINT_SPEED : WALK_SPEED;
    }

    HMM_Vec3 forward_input = (HMM_Vec3){{0, 0, 0}};
    HMM_Vec3 side_input = (HMM_Vec3){{0, 0, 0}};
    if (cam->space_mode) {
        // Space mode: WASD moves in look direction (full 3D, not projected)
        if (cam->key_w) forward_input = HMM_AddV3(forward_input, cam->forward);
        if (cam->key_s) forward_input = HMM_SubV3(forward_input, cam->forward);
        if (cam->key_a) side_input = HMM_AddV3(side_input, flat_right);
        if (cam->key_d) side_input = HMM_SubV3(side_input, flat_right);
    } else {
        if (cam->key_w) forward_input = HMM_AddV3(forward_input, flat_forward);
        if (cam->key_s) forward_input = HMM_SubV3(forward_input, flat_forward);
        if (cam->key_a) side_input = HMM_AddV3(side_input, flat_right);
        if (cam->key_d) side_input = HMM_SubV3(side_input, flat_right);
    }

    float fwd_len = HMM_LenV3(forward_input);
    float side_len = HMM_LenV3(side_input);
    if (fwd_len > 0.001f) forward_input = HMM_NormV3(forward_input);
    if (side_len > 0.001f) side_input = HMM_NormV3(side_input);

    // Check if hex terrain covers our position (for voxel-based collision)
    bool use_hex_collision = (hex && hex->frame_valid && !cam->jetpack_active);
    if (use_hex_collision) {
        // Get current hex ground height to determine if we're within hex range
        int cur_gcol, cur_grow, cur_layer;
        hex_terrain_world_to_hex(hex, cam->position, &cur_gcol, &cur_grow, &cur_layer);
        float cur_ground_r = hex_terrain_ground_height(hex, cur_gcol, cur_grow);
        if (cur_ground_r < 1.0f) use_hex_collision = false; // no chunk loaded here
    }

    // Save position before movement (used by both hex and jetpack collision)
    double old_pos_d[3] = { cam->pos_d[0], cam->pos_d[1], cam->pos_d[2] };

    if (use_hex_collision) {
        // ---- Hex voxel collision: try movement, check walls, slide ----
        // See docs/hex-collision-design.md for physics case documentation.

        // Get current ground info using player's feet altitude as search ceiling.
        // This correctly finds ground inside arches/caves (not the roof above).
        bool moon_hex = (cam->gravity_body >= 0 && solar);
        // Use local_up for moons: the actual ellipsoid normal at the player's position.
        // tangent_up (fixed at grid center) diverges by up to 20°+ on small moons when
        // the player is far from the grid center, causing 50m+ altitude errors.
        // Using local_up for BOTH altitude AND snap is self-consistent — any per-frame
        // drift from the normal changing is corrected by the ground snap each frame.
        HMM_Vec3 hex_up = cam->local_up;
        int cur_gcol, cur_grow, cur_layer;
        hex_terrain_world_to_hex(hex, cam->position, &cur_gcol, &cur_grow, &cur_layer);
        float cur_col_r = hex_terrain_col_base_r(hex, cur_gcol, cur_grow);
        if (cur_col_r < 1.0f) cur_col_r = hex->planet_radius;

        // Feet layer: normal-aligned for moons, radial for planet
        int cur_feet_layer;
        if (moon_hex) {
            HMM_Vec3 surf = hex_col_surface(hex, cur_gcol, cur_grow, cur_col_r);
            float normal_alt = HMM_DotV3(HMM_SubV3(cam->position, surf), hex_up);
            cur_feet_layer = (int)floorf((normal_alt - cam->eye_height) / HEX_HEIGHT);
        } else {
            cur_feet_layer = (int)floor((pos_len_d - (double)cam->eye_height - (double)cur_col_r) / (double)HEX_HEIGHT);
        }
        int search_ceiling = cur_feet_layer + (int)ceilf(AUTO_STEP_HEIGHT / HEX_HEIGHT) + 1;
        float cur_ground_r = hex_terrain_ground_height_below(hex, cur_gcol, cur_grow, search_ceiling);
        if (cur_ground_r < 1.0f) cur_ground_r = hex_terrain_ground_height(hex, cur_gcol, cur_grow);
        // Normal-aligned ground height (offset from ellipsoid surface)
        float cur_ground_h = cur_ground_r - cur_col_r;

        double step = (double)move_speed * dd;

        // Try full movement, then component-wise if blocked
        HMM_Vec3 move_dirs[3];
        int num_dirs = 0;
        HMM_Vec3 full_input = HMM_AddV3(forward_input, side_input);
        if (HMM_LenV3(full_input) > 0.001f) {
            move_dirs[num_dirs++] = HMM_NormV3(full_input);
            // Also prepare component fallbacks for wall sliding
            if (fwd_len > 0.001f && side_len > 0.001f) {
                move_dirs[num_dirs++] = forward_input;
                move_dirs[num_dirs++] = side_input;
            }
        }

        // --- Continuous movement diagnostic (F key toggle) ---
        // Tracks actual position delta per frame to detect snap eating movement.
        if (cam->move_debug && moon_hex) {
            static int diag_throttle = 0;
            static double prev_pos[3] = {0, 0, 0};
            bool diag_print = (++diag_throttle % 30 == 0);

            // Compute actual frame delta
            double dx = cam->pos_d[0] - prev_pos[0];
            double dy = cam->pos_d[1] - prev_pos[1];
            double dz = cam->pos_d[2] - prev_pos[2];
            HMM_Vec3 delta = {{(float)dx, (float)dy, (float)dz}};

            // Decompose into tangent (horizontal) and normal (vertical) components
            float vert_delta = HMM_DotV3(delta, hex_up);
            HMM_Vec3 vert_vec = HMM_MulV3F(hex_up, vert_delta);
            HMM_Vec3 horiz_vec = HMM_SubV3(delta, vert_vec);
            float horiz_speed = HMM_LenV3(horiz_vec) / (dt > 0.0001f ? dt : 0.016f);

            // Check if any movement key is pressed
            bool any_key = cam->key_w || cam->key_s || cam->key_a || cam->key_d;

            // Compute expected step this frame
            float expected_step = move_speed * dt;

            // Check flat_forward/flat_right dot with hex_up (should be ~0 for tangent movement)
            float fwd_up_dot = HMM_DotV3(flat_forward, hex_up);
            float right_up_dot = HMM_DotV3(flat_right, hex_up);

            if (diag_print) {
                printf("[MOVE] keys=%c%c%c%c horiz_spd=%.3f vert_delta=%.4f expected_step=%.4f\n",
                    cam->key_w ? 'W' : '.', cam->key_s ? 'S' : '.',
                    cam->key_a ? 'A' : '.', cam->key_d ? 'D' : '.',
                    horiz_speed, vert_delta, expected_step);
                printf("  fwd.hex_up=%.4f right.hex_up=%.4f num_dirs=%d moved_last=on_ground=%d\n",
                    fwd_up_dot, right_up_dot, num_dirs, cam->on_ground);
                printf("  local_up=(%.4f,%.4f,%.4f) hex_up=(%.4f,%.4f,%.4f) dot=%.6f\n",
                    cam->local_up.X, cam->local_up.Y, cam->local_up.Z,
                    hex_up.X, hex_up.Y, hex_up.Z,
                    HMM_DotV3(cam->local_up, hex_up));
                printf("  pos_delta=(%.5f,%.5f,%.5f) gcol=%d grow=%d ground_h=%.2f\n",
                    (float)dx, (float)dy, (float)dz, cur_gcol, cur_grow, cur_ground_h);
                fflush(stdout);
            }

            prev_pos[0] = cam->pos_d[0];
            prev_pos[1] = cam->pos_d[1];
            prev_pos[2] = cam->pos_d[2];
        }

        bool moved = false;
        static int ml_throttle = 0;
        bool ml_print = (cam->move_debug && moon_hex && (++ml_throttle % 30 == 0));
        for (int attempt = 0; attempt < num_dirs && !moved; attempt++) {
            HMM_Vec3 dir = move_dirs[attempt];

            // Proposed new position
            double test_pos[3] = {
                old_pos_d[0] + (double)dir.X * step,
                old_pos_d[1] + (double)dir.Y * step,
                old_pos_d[2] + (double)dir.Z * step,
            };

            // Probe ahead: check collision at position + margin in movement direction.
            // This prevents the player from reaching the wall column boundary before
            // being blocked (otherwise they visually overlap the wall by one frame's step).
            float probe_margin = 0.3f;
            HMM_Vec3 probe_pos_f = (HMM_Vec3){{
                (float)(test_pos[0] + (double)dir.X * (double)probe_margin),
                (float)(test_pos[1] + (double)dir.Y * (double)probe_margin),
                (float)(test_pos[2] + (double)dir.Z * (double)probe_margin)
            }};
            int dest_gcol, dest_grow, dest_layer;
            hex_terrain_world_to_hex(hex, probe_pos_f, &dest_gcol, &dest_grow, &dest_layer);

            // Find ground at destination, searching down from player's current altitude
            // + step margin. This lets the player walk under arches/caves.
            float dest_ground_r = hex_terrain_ground_height_below(
                hex, dest_gcol, dest_grow, search_ceiling);

            // If no chunk at destination, allow movement (LOD handles it)
            if (dest_ground_r < 1.0f) {
                cam->pos_d[0] = test_pos[0];
                cam->pos_d[1] = test_pos[1];
                cam->pos_d[2] = test_pos[2];
                moved = true;
                if (ml_print) { printf("[MLOOP] ALLOW(no_chunk) a=%d dest(%d,%d) ceil=%d\n", attempt, dest_gcol, dest_grow, search_ceiling); fflush(stdout); }
                break;
            }

            // Step height check: compare normal-aligned ground offsets
            float dest_col_r = hex_terrain_col_base_r(hex, dest_gcol, dest_grow);
            if (dest_col_r < 1.0f) dest_col_r = hex->planet_radius;
            float dest_ground_h = dest_ground_r - dest_col_r;
            float height_diff = dest_ground_h - cur_ground_h;
            if (height_diff > AUTO_STEP_HEIGHT) {
                if (ml_print) { printf("[MLOOP] BLOCK(step) a=%d dest(%d,%d) dgh=%.2f cgh=%.2f diff=%.2f limit=%.2f dgr=%.2f dcr=%.2f ceil=%d\n",
                    attempt, dest_gcol, dest_grow, dest_ground_h, cur_ground_h, height_diff, AUTO_STEP_HEIGHT, dest_ground_r, dest_col_r, search_ceiling); fflush(stdout); }
                continue;
            }

            // Headroom check: PLAYER_HEIGHT layers of air above dest ground
            int dest_ground_layer = (int)floorf(dest_ground_h / HEX_HEIGHT);
            if (!hex_terrain_has_headroom(hex, dest_gcol, dest_grow,
                                          dest_ground_layer + 1, (int)ceilf(PLAYER_HEIGHT))) {
                if (ml_print) { printf("[MLOOP] BLOCK(headroom) a=%d dest(%d,%d) dgl=%d\n", attempt, dest_gcol, dest_grow, dest_ground_layer); fflush(stdout); }
                continue;
            }

            // Movement allowed
            cam->pos_d[0] = test_pos[0];
            cam->pos_d[1] = test_pos[1];
            cam->pos_d[2] = test_pos[2];
            moved = true;
            if (ml_print) { printf("[MLOOP] ALLOW a=%d dest(%d,%d) dgh=%.2f cgh=%.2f diff=%.2f\n", attempt, dest_gcol, dest_grow, dest_ground_h, cur_ground_h, height_diff); fflush(stdout); }
        }

        // Jump / gravity: use local_up (actual surface normal at player position)
        // rather than hex_up (fixed patch normal). hex_up is for altitude/snap only.
        if (cam->key_space && cam->on_ground) {
            cam->velocity = HMM_MulV3F(cam->local_up, JUMP_FORCE);
            cam->on_ground = false;
        }
        if (!cam->on_ground) {
            cam->velocity = HMM_SubV3(cam->velocity, HMM_MulV3F(cam->local_up, cam->gravity * dt));
        }

        // Apply velocity (double precision)
        if (HMM_LenV3(cam->velocity) > 0.001f) {
            cam->pos_d[0] += (double)cam->velocity.X * dd;
            cam->pos_d[1] += (double)cam->velocity.Y * dd;
            cam->pos_d[2] += (double)cam->velocity.Z * dd;
        }

        // Sync float position
        cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

        // Ground collision using hex terrain voxel data (altitude-aware for arches)
        hex_terrain_world_to_hex(hex, cam->position, &cur_gcol, &cur_grow, &cur_layer);
        float snap_col_r = hex_terrain_col_base_r(hex, cur_gcol, cur_grow);
        if (snap_col_r < 1.0f) snap_col_r = hex->planet_radius;

        // Compute feet layer (normal-aligned for moons, radial for planet)
        float snap_normal_alt;  // player altitude above ellipsoid surface (normal-aligned)
        int snap_feet_layer;
        if (moon_hex) {
            HMM_Vec3 surf = hex_col_surface(hex, cur_gcol, cur_grow, snap_col_r);
            snap_normal_alt = HMM_DotV3(HMM_SubV3(cam->position, surf), hex_up);
            snap_feet_layer = (int)floorf((snap_normal_alt - cam->eye_height) / HEX_HEIGHT);
        } else {
            double pr_snap = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                                  cam->pos_d[1]*cam->pos_d[1] +
                                  cam->pos_d[2]*cam->pos_d[2]);
            snap_normal_alt = (float)(pr_snap - (double)snap_col_r);
            snap_feet_layer = (int)floor((pr_snap - (double)cam->eye_height - (double)snap_col_r) / (double)HEX_HEIGHT);
        }

        // When on ground, allow step-up margin so walking auto-steps.
        // When airborne, search only at/below feet — don't find surfaces above the player.
        int snap_ceiling = cam->on_ground
            ? snap_feet_layer + (int)ceilf(AUTO_STEP_HEIGHT / HEX_HEIGHT) + 1
            : snap_feet_layer;

        float ground_r = hex_terrain_ground_height_below(hex, cur_gcol, cur_grow, snap_ceiling);
        if (ground_r < 1.0f) {
            // Fallback to LOD if no hex chunk
            if (lod) ground_r = (float)lod_tree_terrain_height(lod, cam->position);
            else ground_r = pos_len;
        }
        // Normal-aligned ground height above surface
        float ground_h = ground_r - snap_col_r;

        float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);

        // Ceiling collision: if moving upward, check for solid blocks above head
        if (radial_vel > 0.0f) {
            int head_layer = snap_feet_layer + (int)ceilf(PLAYER_HEIGHT / HEX_HEIGHT);
            if (hex_terrain_get_voxel(hex, cur_gcol, cur_grow, head_layer) != 0) {
                // Hit ceiling — kill upward velocity, start falling
                HMM_Vec3 up_component = HMM_MulV3F(cam->local_up, radial_vel);
                cam->velocity = HMM_SubV3(cam->velocity, up_component);
            }
        }

        float feet_alt = snap_normal_alt - cam->eye_height;
        if (feet_alt <= ground_h) {
            // Feet at or below ground: snap altitude along normal (preserve tangent position).
            // Only adjust the component along hex_up, like planet path only scales radial.
            // Setting absolute position would snap to column center, killing lateral movement.
            if (moon_hex) {
                float correction = (ground_h + cam->eye_height) - snap_normal_alt;
                cam->pos_d[0] += (double)(hex_up.X * correction);
                cam->pos_d[1] += (double)(hex_up.Y * correction);
                cam->pos_d[2] += (double)(hex_up.Z * correction);
            } else {
                double pr = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                                 cam->pos_d[1]*cam->pos_d[1] +
                                 cam->pos_d[2]*cam->pos_d[2]);
                double target_r = (double)ground_r + (double)cam->eye_height;
                double scale = target_r / pr;
                cam->pos_d[0] *= scale;
                cam->pos_d[1] *= scale;
                cam->pos_d[2] *= scale;
            }
            if (radial_vel <= 0.0f) {
                cam->velocity = (HMM_Vec3){{0, 0, 0}};
            }
            cam->on_ground = true;
        } else {
            // Feet above ground: freefall (gravity pulls player down naturally)
            cam->on_ground = false;
        }
    } else if (cam->space_mode) {
        // ---- Space flight mode: no gravity, free movement ----
        HMM_Vec3 full_input = HMM_AddV3(forward_input, side_input);
        bool has_input = HMM_LenV3(full_input) > 0.001f;
        if (has_input) {
            HMM_Vec3 full_dir = HMM_NormV3(full_input);
            double step = (double)move_speed * dd;
            cam->pos_d[0] += (double)full_dir.X * step;
            cam->pos_d[1] += (double)full_dir.Y * step;
            cam->pos_d[2] += (double)full_dir.Z * step;
        }

        // Space/Ctrl for vertical thrust (relative to space_up)
        float sm = cam->jetpack_speed_mult;
        if (cam->key_space) {
            float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
            if (radial_vel < JETPACK_MAX_SPEED * sm) {
                cam->velocity = HMM_AddV3(cam->velocity,
                    HMM_MulV3F(cam->local_up, JETPACK_THRUST * sm * dt));
            }
        } else if (cam->key_ctrl) {
            float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
            if (radial_vel > -JETPACK_MAX_SPEED * sm) {
                cam->velocity = HMM_SubV3(cam->velocity,
                    HMM_MulV3F(cam->local_up, JETPACK_THRUST * sm * dt));
            }
        }
        cam->velocity = HMM_MulV3F(cam->velocity, JETPACK_DRAG);
        cam->on_ground = false;

        // Apply velocity
        if (HMM_LenV3(cam->velocity) > 0.001f) {
            cam->pos_d[0] += (double)cam->velocity.X * dd;
            cam->pos_d[1] += (double)cam->velocity.Y * dd;
            cam->pos_d[2] += (double)cam->velocity.Z * dd;
        }
        cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

    } else if (cam->gravity_body >= 0 && solar) {
        // ---- On a moon: simplified gravity + ground collision ----
        HMM_Vec3 full_input = HMM_AddV3(forward_input, side_input);
        bool has_input = HMM_LenV3(full_input) > 0.001f;
        if (has_input) {
            HMM_Vec3 full_dir = HMM_NormV3(full_input);
            double step = (double)move_speed * dd;
            cam->pos_d[0] += (double)full_dir.X * step;
            cam->pos_d[1] += (double)full_dir.Y * step;
            cam->pos_d[2] += (double)full_dir.Z * step;
        }

        // Jetpack / Jump / Gravity (same as planet but toward moon center)
        if (cam->jetpack_active) {
            float sm = cam->jetpack_speed_mult;
            if (cam->key_space) {
                float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
                if (radial_vel < JETPACK_MAX_SPEED * sm) {
                    cam->velocity = HMM_AddV3(cam->velocity,
                        HMM_MulV3F(cam->local_up, JETPACK_THRUST * sm * dt));
                }
            } else if (cam->key_ctrl) {
                float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
                if (radial_vel > -JETPACK_MAX_SPEED * sm) {
                    cam->velocity = HMM_SubV3(cam->velocity,
                        HMM_MulV3F(cam->local_up, JETPACK_THRUST * sm * dt));
                }
            } else {
                cam->velocity = HMM_SubV3(cam->velocity,
                    HMM_MulV3F(cam->local_up, cam->gravity * 0.3f * dt));
            }
            cam->velocity = HMM_MulV3F(cam->velocity, JETPACK_DRAG);
            cam->on_ground = false;
        } else {
            if (cam->key_space && cam->on_ground) {
                cam->velocity = HMM_MulV3F(cam->local_up, JUMP_FORCE);
                cam->on_ground = false;
            }
            if (!cam->on_ground) {
                cam->velocity = HMM_SubV3(cam->velocity,
                    HMM_MulV3F(cam->local_up, cam->gravity * dt));
            }
        }

        // Apply velocity
        if (HMM_LenV3(cam->velocity) > 0.001f) {
            cam->pos_d[0] += (double)cam->velocity.X * dd;
            cam->pos_d[1] += (double)cam->velocity.Y * dd;
            cam->pos_d[2] += (double)cam->velocity.Z * dd;
        }
        cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

        // Moon ground collision: normal-aligned altitude measurement.
        // Surface point is along radial ray at ground_r; altitude measured along
        // ellipsoid normal (local_up) so the player aligns with the ellipsoid shape.
        const CelestialBody* moon = &solar->moons[cam->gravity_body];
        double dx = cam->pos_d[0] - cam->gravity_center_d[0];
        double dy = cam->pos_d[1] - cam->gravity_center_d[1];
        double dz = cam->pos_d[2] - cam->gravity_center_d[2];
        double dist_from_moon = sqrt(dx*dx + dy*dy + dz*dz);
        if (dist_from_moon < 0.001) dist_from_moon = 0.001;
        HMM_Vec3 dir_from_moon = {{(float)(dx/dist_from_moon), (float)(dy/dist_from_moon), (float)(dz/dist_from_moon)}};
        float ground_r = current_ground_r;
        if (ground_r < 1.0f) ground_r = moon_surface_radius(&moon->shape, dir_from_moon);

        // Surface contact point (where radial ray hits noisy surface)
        HMM_Vec3 surface_pt = HMM_MulV3F(dir_from_moon, ground_r);
        // Player offset from moon center
        HMM_Vec3 player_offset = {{(float)dx, (float)dy, (float)dz}};
        // Normal-aligned altitude (distance from surface along ellipsoid normal)
        float normal_alt = HMM_DotV3(HMM_SubV3(player_offset, surface_pt), cam->local_up);
        float feet_alt = normal_alt - cam->eye_height;
        float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);

        if (feet_alt <= GROUND_SNAP_THRESHOLD) {
            if (radial_vel <= 0.0f) {
                // Snap: place feet on surface, body extends along ellipsoid normal
                HMM_Vec3 target = HMM_AddV3(surface_pt,
                    HMM_MulV3F(cam->local_up, cam->eye_height));
                cam->pos_d[0] = cam->gravity_center_d[0] + (double)target.X;
                cam->pos_d[1] = cam->gravity_center_d[1] + (double)target.Y;
                cam->pos_d[2] = cam->gravity_center_d[2] + (double)target.Z;
                cam->velocity = (HMM_Vec3){{0, 0, 0}};
                cam->on_ground = true;
                if (cam->jetpack_active) {
                    cam->jetpack_active = false;
                    cam->jetpack_speed_mult = 1.0f;
                    LOG(PLAYER, "Jetpack OFF (landed on moon)\n");
                }
            }
        } else {
            cam->on_ground = false;
        }

    } else {
        // ---- No hex collision: jetpack or outside hex range (on Tenebris) ----
        HMM_Vec3 full_input = HMM_AddV3(forward_input, side_input);
        bool has_input = HMM_LenV3(full_input) > 0.001f;
        if (has_input) {
            HMM_Vec3 full_dir = HMM_NormV3(full_input);
            double step = (double)move_speed * dd;
            cam->pos_d[0] += (double)full_dir.X * step;
            cam->pos_d[1] += (double)full_dir.Y * step;
            cam->pos_d[2] += (double)full_dir.Z * step;
        }

        // Jetpack / Jump / Gravity
        if (cam->jetpack_active) {
            float sm = cam->jetpack_speed_mult;
            if (cam->key_space) {
                float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
                if (radial_vel < JETPACK_MAX_SPEED * sm) {
                    cam->velocity = HMM_AddV3(cam->velocity,
                        HMM_MulV3F(cam->local_up, JETPACK_THRUST * sm * dt));
                }
            } else if (cam->key_ctrl) {
                float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);
                if (radial_vel > -JETPACK_MAX_SPEED * sm) {
                    cam->velocity = HMM_SubV3(cam->velocity,
                        HMM_MulV3F(cam->local_up, JETPACK_THRUST * sm * dt));
                }
            } else {
                cam->velocity = HMM_SubV3(cam->velocity,
                    HMM_MulV3F(cam->local_up, cam->gravity * 0.3f * dt));
            }
            cam->velocity = HMM_MulV3F(cam->velocity, JETPACK_DRAG);
            cam->on_ground = false;
        } else {
            if (cam->key_space && cam->on_ground) {
                cam->velocity = HMM_MulV3F(cam->local_up, JUMP_FORCE);
                cam->on_ground = false;
            }
            if (!cam->on_ground) {
                cam->velocity = HMM_SubV3(cam->velocity, HMM_MulV3F(cam->local_up, cam->gravity * dt));
            }
        }

        // Apply velocity (double precision)
        if (HMM_LenV3(cam->velocity) > 0.001f) {
            cam->pos_d[0] += (double)cam->velocity.X * dd;
            cam->pos_d[1] += (double)cam->velocity.Y * dd;
            cam->pos_d[2] += (double)cam->velocity.Z * dd;
        }

        // Sync float position for collision queries
        cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

        // Ground collision using LOD terrain height (body-relative)
        if (lod) {
            static double smoothed_ground_r = 0.0;
            double ground_r = lod_tree_terrain_height(lod, cam->position);
            double bdx = cam->pos_d[0] - cam->gravity_center_d[0];
            double bdy = cam->pos_d[1] - cam->gravity_center_d[1];
            double bdz = cam->pos_d[2] - cam->gravity_center_d[2];
            double pr = sqrt(bdx*bdx + bdy*bdy + bdz*bdz);
            double feet_r = pr - (double)cam->eye_height;
            float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);

            if (smoothed_ground_r == 0.0) smoothed_ground_r = ground_r;
            if (ground_r < smoothed_ground_r) {
                smoothed_ground_r = ground_r;
            } else {
                double alpha = 1.0 - exp(-3.0 * dd);
                smoothed_ground_r += (ground_r - smoothed_ground_r) * alpha;
            }

            if (feet_r <= smoothed_ground_r + (double)GROUND_SNAP_THRESHOLD) {
                if (radial_vel <= 0.0f) {
                    double target_r = smoothed_ground_r + (double)cam->eye_height;
                    double scale = target_r / pr;
                    // Scale relative to gravity body center
                    cam->pos_d[0] = cam->gravity_center_d[0] + bdx * scale;
                    cam->pos_d[1] = cam->gravity_center_d[1] + bdy * scale;
                    cam->pos_d[2] = cam->gravity_center_d[2] + bdz * scale;
                    cam->velocity = (HMM_Vec3){{0, 0, 0}};
                    cam->on_ground = true;
                    if (cam->jetpack_active) {
                        cam->jetpack_active = false;
                        cam->jetpack_speed_mult = 1.0f;
                        LOG(PLAYER, "Jetpack OFF (landed)\n");
                    }
                }
            } else {
                cam->on_ground = false;
            }
        }
    }

    // ---- Hex solid voxel rejection (applies to jetpack AND walking) ----
    // Prevent the player from moving INTO solid hex voxels regardless of movement mode.
    // For moons (ellipsoids), use normal-aligned altitude via hex_col_surface + hex_up
    // instead of radial distance, which disagrees with the hex collision block.
    if (hex && hex->frame_valid) {
        HMM_Vec3 test_pos = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};
        int gcol, grow, layer;
        hex_terrain_world_to_hex(hex, test_pos, &gcol, &grow, &layer);

        float vox_col_r = hex_terrain_col_base_r(hex, gcol, grow);
        if (vox_col_r < 1.0f) vox_col_r = hex->planet_radius;

        // Compute feet altitude: normal-aligned for moons, radial for planet
        double feet_alt;
        bool is_moon_vox = (cam->gravity_body >= 0);
        if (is_moon_vox) {
            // Use local_up (same as hex collision block) — tangent_up diverges on small moons
            HMM_Vec3 surf = hex_col_surface(hex, gcol, grow, vox_col_r);
            float normal_alt = HMM_DotV3(HMM_SubV3(test_pos, surf), cam->local_up);
            feet_alt = (double)(normal_alt - cam->eye_height);
        } else {
            double pr_d = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                               cam->pos_d[1]*cam->pos_d[1] +
                               cam->pos_d[2]*cam->pos_d[2]);
            feet_alt = pr_d - (double)cam->eye_height - (double)vox_col_r;
        }
        // Add small upward tolerance (5cm) to prevent false-positive rejection when
        // standing exactly at a layer boundary (snap targets ground_h exactly, but
        // re-computing altitude with different float ops can yield ground_h - epsilon,
        // causing floor() to round down into the solid layer below).
        int feet_layer = (int)floor((feet_alt + 0.05) / (double)HEX_HEIGHT);
        int head_layer = feet_layer + (int)ceil((double)PLAYER_HEIGHT / (double)HEX_HEIGHT);

        bool blocked = false;
        for (int check_l = feet_layer; check_l <= head_layer; check_l++) {
            uint8_t voxel = hex_terrain_get_voxel(hex, gcol, grow, check_l);
            if (voxel != VOXEL_AIR && voxel != VOXEL_TORCH) {
                blocked = true;
                break;
            }
        }

        if (blocked) {
            // Only revert if old position was NOT also blocked (prevents permanent trapping)
            HMM_Vec3 old_pos_f = (HMM_Vec3){{(float)old_pos_d[0], (float)old_pos_d[1], (float)old_pos_d[2]}};
            int old_gcol, old_grow, old_layer;
            hex_terrain_world_to_hex(hex, old_pos_f, &old_gcol, &old_grow, &old_layer);

            float old_col_r = hex_terrain_col_base_r(hex, old_gcol, old_grow);
            if (old_col_r < 1.0f) old_col_r = hex->planet_radius;

            double old_feet_alt;
            if (is_moon_vox) {
                HMM_Vec3 old_surf = hex_col_surface(hex, old_gcol, old_grow, old_col_r);
                float old_normal_alt = HMM_DotV3(HMM_SubV3(old_pos_f, old_surf), cam->local_up);
                old_feet_alt = (double)(old_normal_alt - cam->eye_height);
            } else {
                double old_pr_d = sqrt(old_pos_d[0]*old_pos_d[0] +
                                       old_pos_d[1]*old_pos_d[1] +
                                       old_pos_d[2]*old_pos_d[2]);
                old_feet_alt = old_pr_d - (double)cam->eye_height - (double)old_col_r;
            }
            int old_feet_layer = (int)floor((old_feet_alt + 0.05) / (double)HEX_HEIGHT);
            int old_head_layer = old_feet_layer + (int)ceil((double)PLAYER_HEIGHT / (double)HEX_HEIGHT);

            bool old_blocked = false;
            for (int check_l = old_feet_layer; check_l <= old_head_layer; check_l++) {
                if (hex_terrain_get_voxel(hex, old_gcol, old_grow, check_l) != 0) {
                    old_blocked = true;
                    break;
                }
            }

            if (!old_blocked) {
                cam->pos_d[0] = old_pos_d[0];
                cam->pos_d[1] = old_pos_d[1];
                cam->pos_d[2] = old_pos_d[2];
                cam->velocity = (HMM_Vec3){{0, 0, 0}};
            }
            // If old position was also blocked, let the player through — they need
            // to escape (e.g., jetpack out). Don't trap them forever.
        }
    }

    // ---- Final sync: double → float for rendering ----
    cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

    // ---- Teleport detection: check for unexpected upward position jumps ----
    // Only active when jetpack is OFF and not in space/moon mode.
    if (!cam->jetpack_active && !cam->space_mode && cam->gravity_body < 0) {
        double radius_after = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                                   cam->pos_d[1]*cam->pos_d[1] +
                                   cam->pos_d[2]*cam->pos_d[2]);
        double radius_change = radius_after - radius_before;
        if (radius_change > 0.5) {
            LOG(PLAYER, "========== UNEXPECTED UPWARD TELEPORT ==========\n");
            LOG(PLAYER, "Radius change: %.3f -> %.3f (+%.3f)\n",
                radius_before, radius_after, radius_change);
            LOG(PLAYER, "on_ground=%d key_w=%d key_space=%d dt=%.4f\n",
                cam->on_ground, cam->key_w, cam->key_space, dt);
            LOG(PLAYER, "pos=(%.3f, %.3f, %.3f) vel=(%.3f, %.3f, %.3f)\n",
                cam->position.X, cam->position.Y, cam->position.Z,
                cam->velocity.X, cam->velocity.Y, cam->velocity.Z);
            if (hex && hex->frame_valid) {
                int gcol, grow, layer;
                hex_terrain_world_to_hex(hex, cam->position, &gcol, &grow, &layer);
                float ground_r = hex_terrain_ground_height(hex, gcol, grow);
                LOG(PLAYER, "hex=(%d,%d) layer=%d ground_r=%.3f\n", gcol, grow, layer, ground_r);
            }
            LOG(PLAYER, "================================================\n");
        }
    }

    // ---- View matrix (built directly from basis vectors) ----
    // AVOID HMM_LookAt_RH(pos, pos+forward, up) — it does normalize(target-eye)
    // internally, which loses the pitch component of forward due to float
    // rounding at 800km (e.g., 800000+0.12 rounds to 800000.125 → 4% pitch error).
    // Instead, build the rotation matrix directly from the already-precise unit vectors.
    {
        HMM_Vec3 F = cam->forward;
        HMM_Vec3 S = HMM_NormV3(HMM_Cross(F, cam->local_up));
        HMM_Vec3 U = HMM_Cross(S, F);

        cam->view.Elements[0][0] = S.X;
        cam->view.Elements[0][1] = U.X;
        cam->view.Elements[0][2] = -F.X;
        cam->view.Elements[0][3] = 0.0f;

        cam->view.Elements[1][0] = S.Y;
        cam->view.Elements[1][1] = U.Y;
        cam->view.Elements[1][2] = -F.Y;
        cam->view.Elements[1][3] = 0.0f;

        cam->view.Elements[2][0] = S.Z;
        cam->view.Elements[2][1] = U.Z;
        cam->view.Elements[2][2] = -F.Z;
        cam->view.Elements[2][3] = 0.0f;

        // Translation row (stripped by camera-relative rendering, but kept for completeness)
        cam->view.Elements[3][0] = -HMM_DotV3(S, cam->position);
        cam->view.Elements[3][1] = -HMM_DotV3(U, cam->position);
        cam->view.Elements[3][2] = HMM_DotV3(F, cam->position);
        cam->view.Elements[3][3] = 1.0f;
    }

    // ---- Projection ----
    float aspect = sapp_widthf() / sapp_heightf();
    cam->proj = HMM_Perspective_RH_NO(HMM_ToRad(70.0f), aspect, 0.01f, 10000000.0f);

    // ---- Periodic logging ----
    log_frame_counter++;
    if (log_frame_counter % 120 == 0) {
        if (LOG_ENABLE_PLAYER && log_verbose) {
            printf("[PLAYER] pos=(%.3f, %.3f, %.3f) r=%.3f vel=(%.3f, %.3f, %.3f) ground=%d jetpack=%d",
                cam->position.X, cam->position.Y, cam->position.Z,
                pos_len,
                cam->velocity.X, cam->velocity.Y, cam->velocity.Z,
                cam->on_ground, cam->jetpack_active);
            if (lod) {
                double gr = lod_tree_terrain_height(lod, cam->position);
                printf(" ground_r=%.3f feet_r=%.3f",
                    gr, (float)(pos_len_d - (double)cam->eye_height));
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
                #ifdef __APPLE__
                platform_macos_set_mouse_tap(true);
                #endif
            }
            break;

        case SAPP_EVENTTYPE_MOUSE_MOVE:
            if (cam->mouse_locked) {
                // Always accumulate sokol events. On macOS with CGEventTap,
                // tap deltas override these in camera_update if available.
                cam->mouse_dx_accum += ev->mouse_dx;
                cam->mouse_dy_accum += ev->mouse_dy;

                // Diagnostics (works on all platforms via sokol events)
                if (cam->mouse_diag_enabled) {
                    cam->mouse_events_this_frame++;
                    uint64_t now = stm_now();
                    if (cam->mouse_last_event_time != 0) {
                        float gap = (float)stm_ms(stm_diff(now, cam->mouse_last_event_time));
                        if (gap > cam->mouse_max_gap_ms) cam->mouse_max_gap_ms = gap;
                    }
                    cam->mouse_last_event_time = now;
                    float delta = sqrtf(ev->mouse_dx * ev->mouse_dx + ev->mouse_dy * ev->mouse_dy);
                    if (delta > cam->mouse_max_delta) cam->mouse_max_delta = delta;
                }
            }
            break;

        case SAPP_EVENTTYPE_MOUSE_SCROLL:
            if (cam->jetpack_active || cam->space_mode) {
                // Scroll wheel adjusts jetpack speed multiplier (floor = 1.0)
                float scroll = ev->scroll_y;
                if (scroll > 0.0f) {
                    cam->jetpack_speed_mult *= 1.25f;  // 4 clicks = ~2.4x
                } else if (scroll < 0.0f) {
                    cam->jetpack_speed_mult /= 1.25f;
                    if (cam->jetpack_speed_mult < 1.0f) cam->jetpack_speed_mult = 1.0f;
                }
            }
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
                        if (!cam->jetpack_active) cam->jetpack_speed_mult = 1.0f;
                        LOG(PLAYER, "Jetpack %s\n", cam->jetpack_active ? "ON" : "OFF");
                        cam->last_space_time = 0;  // Reset so held key can't re-toggle
                    } else {
                        cam->last_space_time = now;
                    }
                }
                cam->key_space = true;
            }
            if (ev->key_code == SAPP_KEYCODE_Q) cam->key_q = true;
            if (ev->key_code == SAPP_KEYCODE_E) cam->key_e = true;
            if (ev->key_code == SAPP_KEYCODE_LEFT_SHIFT) cam->key_shift = true;
            if (ev->key_code == SAPP_KEYCODE_LEFT_CONTROL || ev->key_code == SAPP_KEYCODE_C) cam->key_ctrl = true;
            if (ev->key_code == SAPP_KEYCODE_ESCAPE) {
                cam->mouse_locked = false;
                sapp_lock_mouse(false);
                #ifdef __APPLE__
                platform_macos_set_mouse_tap(false);
                #endif
            }
            break;

        case SAPP_EVENTTYPE_KEY_UP:
            if (ev->key_code == SAPP_KEYCODE_W) cam->key_w = false;
            if (ev->key_code == SAPP_KEYCODE_S) cam->key_s = false;
            if (ev->key_code == SAPP_KEYCODE_A) cam->key_a = false;
            if (ev->key_code == SAPP_KEYCODE_D) cam->key_d = false;
            if (ev->key_code == SAPP_KEYCODE_SPACE) cam->key_space = false;
            if (ev->key_code == SAPP_KEYCODE_Q) cam->key_q = false;
            if (ev->key_code == SAPP_KEYCODE_E) cam->key_e = false;
            if (ev->key_code == SAPP_KEYCODE_LEFT_SHIFT) cam->key_shift = false;
            if (ev->key_code == SAPP_KEYCODE_LEFT_CONTROL || ev->key_code == SAPP_KEYCODE_C) cam->key_ctrl = false;
            break;

        default:
            break;
    }
}

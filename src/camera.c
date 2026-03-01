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
#define PLAYER_EYE_HEIGHT 1.7f  // Eye height above feet (meters)
#define PLAYER_HEIGHT   2.0f    // Total player height (2m capsule)
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
    cam->gravity = BASE_GRAVITY;
    cam->on_ground = false;
    cam->jetpack_active = false;
    cam->jetpack_speed_mult = 1.0f;
    cam->last_space_time = 0;

    cam->view = HMM_M4D(1.0f);
    cam->proj = HMM_M4D(1.0f);
}

void camera_update(Camera* cam, Planet* planet, const LodTree* lod,
                   const HexTerrain* hex, float dt) {
    double dd = (double)dt;

    // ---- Sync float position from double (authoritative) ----
    cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

    // Local "up" computed from double-precision position (smooth even at 800 km)
    double pos_len_d = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                            cam->pos_d[1]*cam->pos_d[1] +
                            cam->pos_d[2]*cam->pos_d[2]);
    if (pos_len_d < 0.001) pos_len_d = 0.001;
    float pos_len = (float)pos_len_d;
    cam->local_up = (HMM_Vec3){{
        (float)(cam->pos_d[0] / pos_len_d),
        (float)(cam->pos_d[1] / pos_len_d),
        (float)(cam->pos_d[2] / pos_len_d)
    }};

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
            // Rodrigues rotation around axis = prev_up × local_up
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
            // Degenerate: fall back to global reference
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
    float current_ground_r = pos_len;
    if (lod) {
        current_ground_r = (float)lod_tree_terrain_height(lod, cam->position);
    } else if (planet) {
        current_ground_r = planet->radius + planet->sea_level * planet->layer_thickness;
    }

    // ---- Cap jetpack speed by altitude ----
    // <1km: max 25x, 1-5km: max 200x, 5-50km: max 500x, >50km: unlimited
    // Once you drop below 50km, speed is clamped to the tier cap.
    {
        float altitude = pos_len - current_ground_r;
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
    if (cam->jetpack_active) {
        move_speed = JETPACK_FLY_SPEED * cam->jetpack_speed_mult;
    } else {
        move_speed = cam->key_shift ? SPRINT_SPEED : WALK_SPEED;
    }

    HMM_Vec3 forward_input = (HMM_Vec3){{0, 0, 0}};
    HMM_Vec3 side_input = (HMM_Vec3){{0, 0, 0}};
    if (cam->key_w) forward_input = HMM_AddV3(forward_input, flat_forward);
    if (cam->key_s) forward_input = HMM_SubV3(forward_input, flat_forward);
    if (cam->key_a) side_input = HMM_AddV3(side_input, flat_right);
    if (cam->key_d) side_input = HMM_SubV3(side_input, flat_right);

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

    if (use_hex_collision) {
        // ---- Hex voxel collision: try movement, check walls, slide ----
        double old_pos_d[3] = { cam->pos_d[0], cam->pos_d[1], cam->pos_d[2] };

        // Get current ground info
        int cur_gcol, cur_grow, cur_layer;
        hex_terrain_world_to_hex(hex, cam->position, &cur_gcol, &cur_grow, &cur_layer);
        float cur_ground_r = hex_terrain_ground_height(hex, cur_gcol, cur_grow);

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

        bool moved = false;
        for (int attempt = 0; attempt < num_dirs && !moved; attempt++) {
            HMM_Vec3 dir = move_dirs[attempt];

            // Proposed new position
            double test_pos[3] = {
                old_pos_d[0] + (double)dir.X * step,
                old_pos_d[1] + (double)dir.Y * step,
                old_pos_d[2] + (double)dir.Z * step,
            };

            HMM_Vec3 test_pos_f = (HMM_Vec3){{(float)test_pos[0], (float)test_pos[1], (float)test_pos[2]}};
            int dest_gcol, dest_grow, dest_layer;
            hex_terrain_world_to_hex(hex, test_pos_f, &dest_gcol, &dest_grow, &dest_layer);
            float dest_ground_r = hex_terrain_ground_height(hex, dest_gcol, dest_grow);

            // If no chunk at destination, allow movement (LOD handles it)
            if (dest_ground_r < 1.0f) {
                cam->pos_d[0] = test_pos[0];
                cam->pos_d[1] = test_pos[1];
                cam->pos_d[2] = test_pos[2];
                moved = true;
                break;
            }

            // Step height check: can we step onto the destination ground?
            float height_diff = dest_ground_r - cur_ground_r;
            if (height_diff > AUTO_STEP_HEIGHT) {
                // Wall too tall to step over — BLOCKED
                continue;
            }

            // Headroom check: 2 layers of air above ground at destination
            int dest_ground_layer = (int)floorf((dest_ground_r - hex->planet_radius) / HEX_HEIGHT);
            if (!hex_terrain_has_headroom(hex, dest_gcol, dest_grow,
                                          dest_ground_layer + 1, (int)ceilf(PLAYER_HEIGHT))) {
                // Ceiling too low — BLOCKED
                continue;
            }

            // Movement allowed
            cam->pos_d[0] = test_pos[0];
            cam->pos_d[1] = test_pos[1];
            cam->pos_d[2] = test_pos[2];
            moved = true;
        }

        // Jump / gravity (same as before)
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

        // Ground collision using hex terrain voxel data
        hex_terrain_world_to_hex(hex, cam->position, &cur_gcol, &cur_grow, &cur_layer);
        float ground_r = hex_terrain_ground_height(hex, cur_gcol, cur_grow);
        if (ground_r < 1.0f) {
            // Fallback to LOD if no hex chunk
            if (lod) ground_r = (float)lod_tree_terrain_height(lod, cam->position);
            else ground_r = pos_len;
        }

        double pr = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                         cam->pos_d[1]*cam->pos_d[1] +
                         cam->pos_d[2]*cam->pos_d[2]);
        double feet_r = pr - (double)cam->eye_height;
        float radial_vel = HMM_DotV3(cam->velocity, cam->local_up);

        if (feet_r <= (double)ground_r + (double)GROUND_SNAP_THRESHOLD) {
            if (radial_vel <= 0.0f) {
                double target_r = (double)ground_r + (double)cam->eye_height;
                double scale = target_r / pr;
                cam->pos_d[0] *= scale;
                cam->pos_d[1] *= scale;
                cam->pos_d[2] *= scale;
                cam->velocity = (HMM_Vec3){{0, 0, 0}};
                cam->on_ground = true;
            }
        } else {
            cam->on_ground = false;
        }
    } else {
        // ---- No hex collision: jetpack or outside hex range ----
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

        // Ground collision using LOD terrain height
        if (lod) {
            static double smoothed_ground_r = 0.0;
            double ground_r = lod_tree_terrain_height(lod, cam->position);
            double pr = sqrt(cam->pos_d[0]*cam->pos_d[0] +
                             cam->pos_d[1]*cam->pos_d[1] +
                             cam->pos_d[2]*cam->pos_d[2]);
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
                    cam->pos_d[0] *= scale;
                    cam->pos_d[1] *= scale;
                    cam->pos_d[2] *= scale;
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

    // ---- Final sync: double → float for rendering ----
    cam->position = (HMM_Vec3){{(float)cam->pos_d[0], (float)cam->pos_d[1], (float)cam->pos_d[2]}};

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
            if (cam->jetpack_active) {
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
            if (ev->key_code == SAPP_KEYCODE_LEFT_SHIFT) cam->key_shift = true;
            if (ev->key_code == SAPP_KEYCODE_LEFT_CONTROL) cam->key_ctrl = true;
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
            if (ev->key_code == SAPP_KEYCODE_LEFT_CONTROL) cam->key_ctrl = false;
            break;

        default:
            break;
    }
}

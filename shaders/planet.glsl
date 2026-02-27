// planet.glsl -- Hex planet terrain shader
// Features: tenebris-style lighting, aerial perspective fog, camera-relative rendering,
// logarithmic depth buffer for massive draw distances.
// Double-float emulation: camera offset split into high+low for sub-cm precision at any distance.
@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs planet_vs
layout(binding=0) uniform planet_vs_params {
    mat4 mvp;
    vec4 camera_offset;      // xyz = camera position (float part)
    vec4 camera_offset_low;  // xyz = residual from double precision (sub-ULP correction)
    vec4 log_depth;          // x = Fcoef (0 = disabled), y = far_plane, z = z_bias (-1 for GL, 0 for D3D11/Metal/WebGPU)
};

layout(location=0) in vec3 a_position;
layout(location=1) in vec3 a_normal;
layout(location=2) in vec3 a_color;

out vec3 fs_normal;
out vec3 fs_color;
out vec3 fs_cam_rel_pos;

void main() {
    // Double-float camera-relative subtraction:
    // (a_position - offset_high) is precise for nearby verts (similar magnitude cancels)
    // then subtracting offset_low adds back the sub-ULP residual from the double split
    vec3 cam_rel_pos = (a_position - camera_offset.xyz) - camera_offset_low.xyz;
    gl_Position = mvp * vec4(cam_rel_pos, 1.0);

    // Logarithmic depth buffer (GPU Gems style)
    float Fcoef = log_depth.x;
    if (Fcoef > 0.0) {
        float w = gl_Position.w;
        gl_Position.z = (log2(max(1e-6, 1.0 + w)) * Fcoef + log_depth.z) * w;
    }

    fs_normal = a_normal;
    fs_color = a_color;
    fs_cam_rel_pos = cam_rel_pos;  // Already camera-relative — precise for fog
}
@end

@fs planet_fs
layout(binding=1) uniform planet_fs_params {
    vec4 sun_direction;     // xyz = normalized sun dir
    vec4 camera_pos;        // xyz = camera world position (for surface_dir computation)
    vec4 fog_params;        // x = fog_density, y = fog_max_distance, z = inscatter_strength
    vec4 lod_debug;         // x = LOD depth (0 = off/normal), y = max_depth for color mapping
};

in vec3 fs_normal;
in vec3 fs_color;
in vec3 fs_cam_rel_pos;

out vec4 frag_color;

void main() {
    vec3 N = normalize(fs_normal);
    vec3 L = normalize(sun_direction.xyz);
    // Reconstruct world direction from camera-relative position
    vec3 world_pos_approx = fs_cam_rel_pos + camera_pos.xyz;
    vec3 surface_dir = normalize(world_pos_approx);

    // Smooth terminator (from tenebris): gradual day/night transition
    float sun_facing = dot(surface_dir, L);
    float sun_brightness = smoothstep(-0.1, 0.3, sun_facing);

    // Ambient: transitions from dim starlight (night) to fill light (day)
    float ambient = mix(0.06, 0.35, sun_brightness);

    // Directional: Lambert diffuse, only on day side
    float ndotl = max(0.0, dot(N, L));
    float directional = ndotl * sun_brightness * 0.7;

    float brightness = ambient + directional;

    // LOD debug mode: color by depth level
    vec3 base_color = fs_color;
    if (lod_debug.x > 0.0) {
        float t = clamp(lod_debug.x / max(lod_debug.y, 1.0), 0.0, 1.0);
        base_color.r = clamp(1.0 - t * 2.0, 0.0, 1.0) + clamp(t * 4.0 - 3.0, 0.0, 1.0);
        base_color.g = clamp(t * 2.0, 0.0, 1.0) - clamp(t * 2.0 - 1.0, 0.0, 1.0);
        base_color.b = clamp(t * 2.0 - 1.0, 0.0, 1.0);
    }
    vec3 terrain_color = base_color * brightness;

    // ---- Aerial perspective fog ----
    // Use camera-relative position directly (already precise, no large-float subtraction)
    float dist = length(fs_cam_rel_pos);
    float fog_density = fog_params.x;
    float fog_max_dist = fog_params.y;
    float inscatter_strength = fog_params.z;

    float normalized_dist = dist / max(fog_max_dist, 1.0);
    float fog_amount = 1.0 - exp(-(normalized_dist * fog_density) * (normalized_dist * fog_density));
    fog_amount = clamp(fog_amount, 0.0, 0.85);

    vec3 day_fog = vec3(0.35, 0.55, 0.85);
    vec3 night_fog = vec3(0.01, 0.01, 0.03);
    vec3 fog_color = mix(night_fog, day_fog, sun_brightness) * inscatter_strength;

    frag_color = vec4(mix(terrain_color, fog_color, fog_amount), 1.0);
}
@end

@program planet planet_vs planet_fs

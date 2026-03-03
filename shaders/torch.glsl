// torch.glsl -- Torch model shader with flame vertex animation
// Camera-relative rendering + log depth (same pattern as hex_terrain)
@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs torch_vs
layout(binding=0) uniform torch_vs_params {
    mat4 vp;                // View-projection matrix
    vec4 camera_offset;     // xyz = camera position (float part)
    vec4 camera_offset_low; // xyz = residual from double precision
    vec4 log_depth;         // x = Fcoef, y = far_plane, z = z_bias
    vec4 torch_pos;         // xyz = world position of this torch, w = time
    vec4 torch_up;          // xyz = surface normal (local up), w = flicker_amount
};

layout(location=0) in vec3 a_position;    // Model-space position (Y-up)
layout(location=1) in vec3 a_normal;      // Model-space normal
layout(location=2) in vec3 a_color;       // Vertex color
layout(location=3) in float a_anim_weight; // 0 = static, 1 = flame

out vec3 fs_color;
out vec3 fs_normal;
out vec3 fs_cam_rel_pos;

void main() {
    float t = torch_pos.w;
    float flickerAmount = torch_up.w;

    // Flame animation (matching tenebris approach)
    vec3 pos = a_position;
    if (a_anim_weight > 0.0) {
        float flicker = sin(t) * sin(t * 2.3) * sin(t * 0.7);
        float flameBaseY = 0.5;
        float heightAboveBase = pos.y - flameBaseY;
        float scaleY = 1.0 + flicker * flickerAmount * a_anim_weight;
        pos.y = flameBaseY + heightAboveBase * scaleY;

        float heightFactor = max(0.0, heightAboveBase * 4.0);
        pos.x += sin(t * 1.7) * 0.008 * heightFactor * a_anim_weight;
        pos.z += cos(t * 1.3) * 0.008 * heightFactor * a_anim_weight;
    }

    // Build orientation matrix from surface normal (torch_up)
    // torch stands perpendicular to the planet surface
    vec3 up = normalize(torch_up.xyz);
    vec3 world_y = vec3(0.0, 1.0, 0.0);
    if (abs(dot(up, world_y)) > 0.99) world_y = vec3(1.0, 0.0, 0.0);
    vec3 right = normalize(cross(world_y, up));
    vec3 forward = cross(up, right);

    // Transform model position to world (model Y-up → surface normal up)
    vec3 world_pos = torch_pos.xyz + right * pos.x + up * pos.y + forward * pos.z;
    vec3 world_normal = right * a_normal.x + up * a_normal.y + forward * a_normal.z;

    // Camera-relative rendering
    vec3 cam_rel_pos = (world_pos - camera_offset.xyz) - camera_offset_low.xyz;
    gl_Position = vp * vec4(cam_rel_pos, 1.0);

    // Logarithmic depth
    float Fcoef = log_depth.x;
    if (Fcoef > 0.0) {
        float w = gl_Position.w;
        gl_Position.z = (log2(max(1e-6, 1.0 + w)) * Fcoef + log_depth.z) * w;
    }

    fs_color = a_color;
    fs_normal = world_normal;
    fs_cam_rel_pos = cam_rel_pos;
}
@end

@fs torch_fs
layout(binding=1) uniform torch_fs_params {
    vec4 sun_direction;  // xyz = normalized sun dir
};

in vec3 fs_color;
in vec3 fs_normal;
in vec3 fs_cam_rel_pos;

out vec4 frag_color;

void main() {
    vec3 N = normalize(fs_normal);
    vec3 L = normalize(sun_direction.xyz);

    // Simple Lambert diffuse + ambient
    float ndotl = max(0.0, dot(N, L));
    vec3 ambient = vec3(0.3);
    vec3 lit_color = fs_color * (ambient + vec3(0.7) * ndotl);

    frag_color = vec4(lit_color, 1.0);
}
@end

@program torch torch_vs torch_fs

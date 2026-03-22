// agent.glsl -- Simple shaded mesh for AI agents
// Same camera-relative rendering + log depth as planet.glsl
// No textures, no atmosphere, just vertex colors + sun lighting.
@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs agent_vs
layout(binding=0) uniform agent_vs_params {
    mat4 mvp;
    vec4 camera_offset;
    vec4 camera_offset_low;
    vec4 log_depth;
};

layout(location=0) in vec3 a_position;
layout(location=1) in vec3 a_normal;
layout(location=2) in vec3 a_color;

out vec3 fs_normal;
out vec3 fs_color;

void main() {
    vec3 cam_rel_pos = (a_position - camera_offset.xyz) - camera_offset_low.xyz;
    gl_Position = mvp * vec4(cam_rel_pos, 1.0);

    float Fcoef = log_depth.x;
    if (Fcoef > 0.0) {
        float w = gl_Position.w;
        gl_Position.z = (log2(max(1e-6, 1.0 + w)) * Fcoef + log_depth.z) * w;
    }

    fs_normal = a_normal;
    fs_color = a_color;
}
@end

@fs agent_fs
layout(binding=1) uniform agent_fs_params {
    vec4 sun_direction;
};

in vec3 fs_normal;
in vec3 fs_color;

out vec4 frag_color;

void main() {
    vec3 N = normalize(fs_normal);
    vec3 L = normalize(sun_direction.xyz);

    float ndotl = max(0.0, dot(N, L));
    vec3 ambient = vec3(0.3);
    vec3 diffuse = vec3(0.7) * ndotl;

    frag_color = vec4(fs_color * (ambient + diffuse), 1.0);
}
@end

@program agent agent_vs agent_fs

// highlight.glsl -- Wireframe highlight for hex selection
// Minimal shader: camera-relative position + log depth, flat color output.
// Used with SG_PRIMITIVETYPE_LINES for hex outline rendering.
@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs highlight_vs
layout(binding=0) uniform highlight_vs_params {
    mat4 mvp;
    vec4 camera_offset;      // xyz = camera position (float part)
    vec4 camera_offset_low;  // xyz = residual from double precision
    vec4 log_depth;          // x = Fcoef, y = far_plane, z = z_bias
};

layout(location=0) in vec3 a_position;

void main() {
    vec3 cam_rel_pos = (a_position - camera_offset.xyz) - camera_offset_low.xyz;
    gl_Position = mvp * vec4(cam_rel_pos, 1.0);

    // Logarithmic depth buffer (same as planet.glsl)
    float Fcoef = log_depth.x;
    if (Fcoef > 0.0) {
        float w = gl_Position.w;
        gl_Position.z = (log2(max(1e-6, 1.0 + w)) * Fcoef + log_depth.z) * w;
    }
}
@end

@fs highlight_fs
layout(binding=1) uniform highlight_fs_params {
    vec4 color;  // rgba highlight color
};

out vec4 frag_color;

void main() {
    frag_color = color;
}
@end

@program highlight highlight_vs highlight_fs

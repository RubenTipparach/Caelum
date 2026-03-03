// hotbar.glsl -- Screen-space textured quad shader for block selection UI
@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs hotbar_vs
layout(binding=0) uniform hotbar_vs_params {
    mat4 ortho;
};

layout(location=0) in vec2 a_position;
layout(location=1) in vec2 a_uv;

out vec2 fs_uv;

void main() {
    gl_Position = ortho * vec4(a_position, 0.0, 1.0);
    fs_uv = a_uv;
}
@end

@fs hotbar_fs
layout(binding=1) uniform hotbar_fs_params {
    vec4 tint;  // xyz = color tint, w = alpha
};

layout(binding=0) uniform texture2D hotbar_tex;
layout(binding=0) uniform sampler hotbar_smp;

in vec2 fs_uv;
out vec4 frag_color;

void main() {
    vec4 tex = texture(sampler2D(hotbar_tex, hotbar_smp), fs_uv);
    frag_color = tex * tint;
}
@end

@program hotbar hotbar_vs hotbar_fs

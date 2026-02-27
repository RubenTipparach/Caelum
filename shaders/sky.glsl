// sky.glsl — Fullscreen sky shader with Voronoi starfield + nebula
// Adapted from tenebris skybox.frag for sokol-shdc

@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs sky_vs

layout(binding=0) uniform sky_vs_params {
    mat4 inv_vp;
    vec4 camera_pos;  // xyz = camera world position
};

layout(location=0) in vec2 a_pos;

out vec3 view_dir;

void main() {
    gl_Position = vec4(a_pos, 0.0, 1.0);

    // Reconstruct world-space view direction from clip coordinates
    vec4 far_pt = inv_vp * vec4(a_pos, 1.0, 1.0);
    view_dir = far_pt.xyz / far_pt.w - camera_pos.xyz;
}
@end

@fs sky_fs

in vec3 view_dir;
out vec4 frag_color;

// ---- Star parameters (from tenebris) ----
const float STAR_SHARPNESS = 0.02;
const float STAR_BRIGHTNESS_THRESHOLD = 0.7;
const vec3 SPACE_COLOR = vec3(0.0, 0.0, 0.02);
const vec3 STAR_COLOR = vec3(1.0, 0.98, 0.95);

// ---- Nebula parameters ----
const vec3 NEBULA_COLOR_1 = vec3(0.4, 0.1, 0.6);   // Purple
const vec3 NEBULA_COLOR_2 = vec3(0.1, 0.3, 0.5);   // Blue
const vec3 NEBULA_COLOR_3 = vec3(0.6, 0.2, 0.3);   // Red/pink
const float NEBULA_INTENSITY = 0.15;
const float NEBULA_SCALE = 2.0;

// ---- Hash functions (from tenebris) ----

float hash_f(float x) {
    return fract(x + 1.3215 * 1.8152);
}

float hash_3f(vec3 a) {
    return fract((hash_f(a.z * 42.8883) + hash_f(a.y * 36.9125) + hash_f(a.x * 65.4321)) * 291.1257);
}

vec3 rehash_3f(float x) {
    return vec3(
        hash_f(((x + 0.5283) * 59.3829) * 274.3487),
        hash_f(((x + 0.8192) * 83.6621) * 345.3871),
        hash_f(((x + 0.2157) * 36.6521) * 458.3971)
    );
}

// ---- Value noise for nebula ----

float hash31(vec3 p) {
    p = fract(p * vec3(443.897, 441.423, 437.195));
    p += dot(p, p.yzx + 19.19);
    return fract((p.x + p.y) * p.z);
}

float noise3d(vec3 p) {
    vec3 ip = floor(p);
    vec3 fp = fract(p);
    fp = fp * fp * (3.0 - 2.0 * fp);

    float n000 = hash31(ip);
    float n100 = hash31(ip + vec3(1.0, 0.0, 0.0));
    float n010 = hash31(ip + vec3(0.0, 1.0, 0.0));
    float n110 = hash31(ip + vec3(1.0, 1.0, 0.0));
    float n001 = hash31(ip + vec3(0.0, 0.0, 1.0));
    float n101 = hash31(ip + vec3(1.0, 0.0, 1.0));
    float n011 = hash31(ip + vec3(0.0, 1.0, 1.0));
    float n111 = hash31(ip + vec3(1.0, 1.0, 1.0));

    float n00 = mix(n000, n100, fp.x);
    float n01 = mix(n001, n101, fp.x);
    float n10 = mix(n010, n110, fp.x);
    float n11 = mix(n011, n111, fp.x);

    float n0 = mix(n00, n10, fp.y);
    float n1 = mix(n01, n11, fp.y);

    return mix(n0, n1, fp.z);
}

float fbm3(vec3 p) {
    float v = 0.0;
    v += noise3d(p) * 0.5;
    v += noise3d(p * 2.0) * 0.25;
    v += noise3d(p * 4.0) * 0.125;
    v += noise3d(p * 8.0) * 0.0625;
    return v;
}

// ---- Voronoi starfield (from tenebris) ----

float starfield(vec3 dir, float density) {
    vec3 pos = dir * density;

    float minDist = 9999.0;
    float cellId = 0.0;

    for (int x = -1; x < 2; x++) {
        for (int y = -1; y < 2; y++) {
            for (int z = -1; z < 2; z++) {
                vec3 cellPos = floor(pos) + vec3(float(x), float(y), float(z));
                float h = hash_3f(cellPos);
                vec3 starPos = rehash_3f(h) + cellPos;

                vec3 diff = pos - starPos;
                float d = dot(diff, diff);  // squared distance
                if (d < minDist) {
                    minDist = d;
                    cellId = h;
                }
            }
        }
    }

    float star = smoothstep(STAR_SHARPNESS, 0.0, minDist);
    float brightness = fract(cellId * 123.456);
    star *= step(STAR_BRIGHTNESS_THRESHOLD, brightness);

    return star * brightness * 1.2;
}

// ---- Nebula ----

vec3 nebula(vec3 dir) {
    vec3 p = dir * NEBULA_SCALE;

    float n1 = fbm3(p);
    float n2 = fbm3(p * 1.5 + vec3(100.0, 0.0, 0.0));
    float n3 = fbm3(p * 2.0 + vec3(0.0, 200.0, 0.0));

    n1 = smoothstep(0.1, 0.6, n1 * 0.5 + 0.5);
    n2 = smoothstep(0.2, 0.7, n2 * 0.5 + 0.5);
    n3 = smoothstep(0.15, 0.5, n3 * 0.5 + 0.5);

    vec3 nebColor = vec3(0.0);
    nebColor += NEBULA_COLOR_1 * n1 * 0.4;
    nebColor += NEBULA_COLOR_2 * n2 * 0.3;
    nebColor += NEBULA_COLOR_3 * n3 * 0.3;

    return nebColor * NEBULA_INTENSITY;
}

// ---- Main ----

void main() {
    vec3 dir = normalize(view_dir);

    // Starfield with three density layers (80 / 160 / 280)
    float stars = 0.0;
    stars += starfield(dir, 80.0) * 1.0;
    stars += starfield(dir, 160.0) * 0.7;
    stars += starfield(dir, 280.0) * 0.5;

    // Nebula background
    vec3 nebColor = nebula(dir);

    // Composite: space + nebula + stars
    vec3 color = SPACE_COLOR + nebColor + stars * STAR_COLOR;

    frag_color = vec4(color, 1.0);
}
@end

@program sky sky_vs sky_fs

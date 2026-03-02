// hex_terrain.glsl -- Textured hex terrain shader
// Based on planet.glsl with texture atlas sampling replacing vertex colors.
// Uses same lighting model (sun, ambient, AO, fog, specular, rim, shadow desat).
@ctype mat4 HMM_Mat4
@ctype vec4 HMM_Vec4

@vs hex_terrain_vs
layout(binding=0) uniform hex_terrain_vs_params {
    mat4 mvp;
    vec4 camera_offset;      // xyz = camera position (float part)
    vec4 camera_offset_low;  // xyz = residual from double precision
    vec4 log_depth;          // x = Fcoef, y = far_plane, z = z_bias
};

layout(location=0) in vec3 a_position;
layout(location=1) in vec3 a_normal;
layout(location=2) in vec2 a_uv;
layout(location=3) in vec3 a_color;
layout(location=4) in float a_sky_light;

out vec3 fs_normal;
out vec2 fs_uv;
out vec3 fs_cam_rel_pos;
out vec3 fs_color;
out float fs_sky_light;

void main() {
    vec3 cam_rel_pos = (a_position - camera_offset.xyz) - camera_offset_low.xyz;
    gl_Position = mvp * vec4(cam_rel_pos, 1.0);

    float Fcoef = log_depth.x;
    if (Fcoef > 0.0) {
        float w = gl_Position.w;
        gl_Position.z = (log2(max(1e-6, 1.0 + w)) * Fcoef + log_depth.z) * w;
    }

    fs_normal = a_normal;
    fs_uv = a_uv;
    fs_cam_rel_pos = cam_rel_pos;
    fs_color = a_color;
    fs_sky_light = a_sky_light;
}
@end

@fs hex_terrain_fs
layout(binding=1) uniform hex_terrain_fs_params {
    vec4 sun_direction;     // xyz = normalized sun dir, w = fog_scale_height
    vec4 camera_pos;        // xyz = camera world position
    vec4 atmos_params;      // x = planet_radius, y = atmos_radius, z = rayleigh_scale, w = sun_intensity
    vec4 lod_debug;         // x = LOD depth (0 = off), y = max_depth
    vec4 dusk_sun_color;    // xyz = warm sunset color
    vec4 day_sun_color;     // xyz = neutral daylight color
    vec4 hex_fade;          // x = fade_start, y = fade_end
};

layout(binding=0) uniform texture2D hex_atlas_tex;
layout(binding=0) uniform sampler hex_atlas_smp;

in vec3 fs_normal;
in vec2 fs_uv;
in vec3 fs_cam_rel_pos;
in vec3 fs_color;
in float fs_sky_light;

out vec4 frag_color;

void main() {
    vec3 N = normalize(fs_normal);
    vec3 L = normalize(sun_direction.xyz);
    vec3 V = -normalize(fs_cam_rel_pos);

    vec3 world_pos_approx = fs_cam_rel_pos + camera_pos.xyz;
    vec3 surface_dir = normalize(world_pos_approx);

    // Sample texture atlas
    vec3 base_color = texture(sampler2D(hex_atlas_tex, hex_atlas_smp), fs_uv).rgb;

    // Distance fade: blend textured hex into vertex terrain color
    {
        float fade_dist = length(fs_cam_rel_pos);
        float fade = smoothstep(hex_fade.x, hex_fade.y, fade_dist);
        base_color = mix(base_color, fs_color, fade);
    }

    // Apply tenebris depth-based sky light (darkens underground blocks uniformly)
    // Applied AFTER texture/vertex blend so both texture and vertex color get darkened
    base_color *= fs_sky_light;

    // Smooth terminator
    float sun_facing = dot(surface_dir, L);
    float sun_brightness = smoothstep(-0.1, 0.3, sun_facing);

    // Atmospheric sun tinting
    float dusk_factor = smoothstep(0.0, 0.35, sun_facing);
    vec3 sunColor = mix(dusk_sun_color.xyz, day_sun_color.xyz, dusk_factor);

    // Normal-based ambient occlusion
    float ao_dot = dot(N, surface_dir);
    float ao_factor = smoothstep(-0.2, 0.8, ao_dot);
    ao_factor = mix(0.45, 1.0, ao_factor);

    // Ambient
    vec3 ambientColor = mix(vec3(0.06), vec3(0.35) * mix(vec3(0.9, 0.7, 0.5), vec3(1.0), dusk_factor), sun_brightness) * ao_factor;

    // Directional Lambert
    float ndotl = max(0.0, dot(N, L));
    vec3 directionalColor = sunColor * ndotl * sun_brightness * 0.7;

    vec3 terrain_color = base_color * (ambientColor + directionalColor);

    // Ocean specular
    float ocean_mask = step(0.5, base_color.b)
                     * step(base_color.r * 1.5, base_color.b)
                     * step(base_color.g * 1.3, base_color.b);
    vec3 H = normalize(L + V);
    float spec = pow(max(0.0, dot(N, H)), 32.0);
    terrain_color += sunColor * spec * ocean_mask * 0.4 * sun_brightness;

    // Rim lighting disabled for hex terrain — voxel faces should shade uniformly

    // Shadow desaturation
    float shadow_amount = 1.0 - ndotl * sun_brightness;
    float saturation = mix(1.0, 0.7, shadow_amount * 0.4);
    vec3 lum = vec3(dot(terrain_color, vec3(0.299, 0.587, 0.114)));
    terrain_color = mix(lum, terrain_color, saturation);

    // Atmospheric fog (Rayleigh aerial perspective)
    const vec3 wl4 = vec3(5.602, 9.473, 19.644);
    float fogScaleH = sun_direction.w;

    float pR = atmos_params.x;
    float aR = atmos_params.y;
    float rScale = atmos_params.z;
    float sInt = atmos_params.w;
    float aThick = aR - pR;

    float camAlt = clamp((length(camera_pos.xyz) - pR) / aThick, 0.0, 1.0);
    float fragAlt = clamp((length(world_pos_approx) - pR) / aThick, 0.0, 1.0);
    float avgDensity = exp(-((camAlt + fragAlt) * 0.5) / fogScaleH);

    float dist = length(fs_cam_rel_pos);
    vec3 tau = rScale * wl4 * avgDensity * (dist / aThick);
    vec3 transmittance = exp(-tau);

    float cosTheta = dot(normalize(fs_cam_rel_pos), L);
    float phaseR = 0.75 * (1.0 + cosTheta * cosTheta);
    vec3 fogEquilibrium = sunColor * sInt * phaseR * sun_brightness;
    fogEquilibrium = vec3(1.0) - exp(-fogEquilibrium);

    terrain_color = terrain_color * transmittance + fogEquilibrium * (vec3(1.0) - transmittance);

    frag_color = vec4(terrain_color, 1.0);
}
@end

@program hex_terrain hex_terrain_vs hex_terrain_fs

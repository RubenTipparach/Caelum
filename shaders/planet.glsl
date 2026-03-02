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
    vec4 sun_direction;     // xyz = normalized sun dir, w = fog_scale_height
    vec4 camera_pos;        // xyz = camera world position, w = hex_suppress_range (0 = disabled)
    vec4 atmos_params;      // x = planet_radius, y = atmos_radius, z = rayleigh_scale, w = sun_intensity
    vec4 lod_debug;         // x = LOD depth (0 = off/normal), y = max_depth for color mapping
    vec4 dusk_sun_color;    // xyz = warm sunset color (from config.yaml)
    vec4 day_sun_color;     // xyz = neutral daylight color (from config.yaml)
};

in vec3 fs_normal;
in vec3 fs_color;
in vec3 fs_cam_rel_pos;

out vec4 frag_color;

void main() {
    // Discard LOD fragments inside hex terrain range (hex terrain renders its own voxels there)
    float hex_suppress = camera_pos.w;
    if (hex_suppress > 0.0 && length(fs_cam_rel_pos) < hex_suppress) {
        discard;
    }

    vec3 N = normalize(fs_normal);
    vec3 L = normalize(sun_direction.xyz);
    vec3 V = -normalize(fs_cam_rel_pos);  // View direction: fragment toward camera

    // Reconstruct world direction from camera-relative position
    vec3 world_pos_approx = fs_cam_rel_pos + camera_pos.xyz;
    vec3 surface_dir = normalize(world_pos_approx);

    // Smooth terminator (from tenebris): gradual day/night transition
    float sun_facing = dot(surface_dir, L);
    float sun_brightness = smoothstep(-0.1, 0.3, sun_facing);

    // Atmospheric sun tinting: when the sun is low on the local horizon,
    // sunlight passes through more atmosphere, scattering out blue → warm orange.
    float dusk_factor = smoothstep(0.0, 0.35, sun_facing);
    vec3 sunColor = mix(dusk_sun_color.xyz, day_sun_color.xyz, dusk_factor);

    // ---- Normal-based ambient occlusion (tenebris crevice darkening) ----
    // Faces pointing outward (aligned with surface_dir) get full ambient.
    // Side/inward faces (hex walls, cliff sides) get reduced ambient.
    float ao_dot = dot(N, surface_dir);
    float ao_factor = smoothstep(-0.2, 0.8, ao_dot);  // 0 at inward → 1 at outward
    ao_factor = mix(0.45, 1.0, ao_factor);             // Range: 45% to 100%

    // Ambient: transitions from dim starlight (night) to fill light (day)
    // Dusk ambient picks up warm tint from the sky
    vec3 ambientColor = mix(vec3(0.06), vec3(0.35) * mix(vec3(0.9, 0.7, 0.5), vec3(1.0), dusk_factor), sun_brightness) * ao_factor;

    // Directional: Lambert diffuse with atmospheric sun tint, only on day side
    float ndotl = max(0.0, dot(N, L));
    vec3 directionalColor = sunColor * ndotl * sun_brightness * 0.7;

    // LOD debug mode: color by depth level
    vec3 base_color = fs_color;
    if (lod_debug.x > 0.0) {
        float t = clamp(lod_debug.x / max(lod_debug.y, 1.0), 0.0, 1.0);
        base_color.r = clamp(1.0 - t * 2.0, 0.0, 1.0) + clamp(t * 4.0 - 3.0, 0.0, 1.0);
        base_color.g = clamp(t * 2.0, 0.0, 1.0) - clamp(t * 2.0 - 1.0, 0.0, 1.0);
        base_color.b = clamp(t * 2.0 - 1.0, 0.0, 1.0);
    }
    vec3 terrain_color = base_color * (ambientColor + directionalColor);

    // ---- Ocean specular (tenebris Blinn-Phong) ----
    // Detect water: blue-dominant vertex color (b > 0.5, b > r*1.5, b > g*1.3)
    float ocean_mask = step(0.5, base_color.b)
                     * step(base_color.r * 1.5, base_color.b)
                     * step(base_color.g * 1.3, base_color.b);
    vec3 H = normalize(L + V);
    float spec = pow(max(0.0, dot(N, H)), 32.0);
    terrain_color += sunColor * spec * ocean_mask * 0.4 * sun_brightness;

    // ---- Rim lighting (tenebris silhouette backlight) ----
    // Subtle glow at terrain edges, emphasizes topography from any angle
    float rim_dot = 1.0 - max(0.0, dot(V, N));
    float rim = pow(rim_dot, 3.0) * 0.12 * sun_brightness;
    terrain_color += vec3(0.35, 0.45, 0.6) * rim;

    // ---- Shadow desaturation (tenebris) ----
    // Colors become more gray in shadow for a more realistic look
    float shadow_amount = 1.0 - ndotl * sun_brightness;
    float saturation = mix(1.0, 0.7, shadow_amount * 0.4);
    vec3 lum = vec3(dot(terrain_color, vec3(0.299, 0.587, 0.114)));
    terrain_color = mix(lum, terrain_color, saturation);

    // ---- Atmospheric fog (Rayleigh aerial perspective) ----
    // Physically-based: extinction dims terrain, inscatter adds sky-colored haze.
    // Naturally blue by day, warm orange at dusk, dark at night.
    const vec3 wl4 = vec3(5.602, 9.473, 19.644); // (1/wavelength)^4 for Rayleigh
    float fogScaleH = sun_direction.w;  // from config.yaml via uniform

    float pR = atmos_params.x;
    float aR = atmos_params.y;
    float rScale = atmos_params.z;
    float sInt = atmos_params.w;
    float aThick = aR - pR;

    // Average density along the camera-to-fragment ray
    float camAlt = clamp((length(camera_pos.xyz) - pR) / aThick, 0.0, 1.0);
    float fragAlt = clamp((length(world_pos_approx) - pR) / aThick, 0.0, 1.0);
    float avgDensity = exp(-((camAlt + fragAlt) * 0.5) / fogScaleH);

    // Optical depth: wavelength-dependent (blue extinguishes faster → warm at distance)
    float dist = length(fs_cam_rel_pos);
    vec3 tau = rScale * wl4 * avgDensity * (dist / aThick);
    vec3 transmittance = exp(-tau);

    // Inscatter: single-scatter Rayleigh approximation
    // Phase function adds view-direction dependence (brighter toward sun)
    float cosTheta = dot(normalize(fs_cam_rel_pos), L);
    float phaseR = 0.75 * (1.0 + cosTheta * cosTheta);
    vec3 fogEquilibrium = sunColor * sInt * phaseR * sun_brightness;
    fogEquilibrium = vec3(1.0) - exp(-fogEquilibrium);  // tone map (can exceed 1.0)

    // Aerial perspective blend: attenuate terrain + add inscattered light
    terrain_color = terrain_color * transmittance + fogEquilibrium * (vec3(1.0) - transmittance);

    frag_color = vec4(terrain_color, 1.0);
}
@end

@program planet planet_vs planet_fs

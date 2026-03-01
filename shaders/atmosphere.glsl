// atmosphere.glsl -- GPU Gems 2 Ch.16 atmospheric scattering (fullscreen ray march)
// Ported from tenebris atmosphere.frag to sokol-shdc format.
// Renders as fullscreen pass with additive blending on top of sky.

@ctype vec4 HMM_Vec4

@vs atmosphere_vs
layout(binding=0) uniform atmosphere_vs_params {
    vec4 cam_right;    // xyz = right vector, w = tan(fov/2)
    vec4 cam_up;       // xyz = up vector,    w = aspect ratio
    vec4 cam_forward;  // xyz = forward vector
    vec4 camera_pos;   // xyz = camera world position (for ray-sphere intersection)
};

layout(location=0) in vec2 a_pos;

out vec3 view_dir;
out vec3 v_camera_pos;

void main() {
    gl_Position = vec4(a_pos, 0.0, 1.0);

    // Reconstruct view direction from camera basis vectors.
    float tan_half_fov = cam_right.w;
    float aspect = cam_up.w;
    view_dir = cam_forward.xyz
             + cam_right.xyz * (a_pos.x * tan_half_fov * aspect)
             + cam_up.xyz    * (a_pos.y * tan_half_fov);
    v_camera_pos = camera_pos.xyz;
}
@end

@fs atmosphere_fs
layout(binding=1) uniform atmosphere_fs_params {
    vec4 sun_direction;     // xyz = normalized sun direction
    vec4 radii;             // x = planet_radius, y = atmosphere_radius, z = scale_height
    vec4 scatter_coeffs;    // x = rayleigh_scale, y = mie_scale, z = mie_g, w = sun_intensity
};

in vec3 view_dir;
in vec3 v_camera_pos;

out vec4 frag_color;

#define PI 3.14159265359
#define NUM_SAMPLES 8
#define NUM_LIGHT_SAMPLES 4

// (1/wavelength)^4 for Rayleigh scattering
// 680nm (red), 550nm (green), 440nm (blue)
const vec3 wavelengthsInv4 = vec3(5.602, 9.473, 19.644);

// Scale height: set from uniform radii.z at start of main() (configurable via config.yaml)
float scaleHeight;

// Simple hash for dithering to reduce banding artifacts
float hash(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453);
}

// Rayleigh phase function: symmetric scattering
float rayleighPhase(float cosTheta) {
    return 0.75 * (1.0 + cosTheta * cosTheta);
}

// Mie phase function (Henyey-Greenstein): forward-peaked scattering
float miePhase(float cosTheta, float g) {
    float g2 = g * g;
    float denom = 1.0 + g2 - 2.0 * g * cosTheta;
    return (1.0 - g2) / (4.0 * PI * pow(max(denom, 0.0001), 1.5));
}

// Ray-sphere intersection (sphere centered at origin)
// Returns vec2(near, far) distances, or vec2(-1) if no hit
vec2 raySphereIntersect(vec3 ro, vec3 rd, float radius) {
    float b = dot(ro, rd);
    float c = dot(ro, ro) - radius * radius;
    float disc = b * b - c;
    if (disc < 0.0) return vec2(-1.0);
    float sqrtDisc = sqrt(disc);
    return vec2(-b - sqrtDisc, -b + sqrtDisc);
}

// Atmospheric density at a given normalized altitude (0=surface, 1=edge)
float getDensity(float altitude) {
    return exp(-altitude / scaleHeight);
}

// Optical depth along a ray segment (sun ray toward atmosphere edge)
// Returns depth in normalized units (divided by atmosThickness) for scale-independence
float computeOpticalDepth(vec3 rayOrigin, vec3 rayDir, float rayLength,
                          float pRadius, float aRadius) {
    float atmosThickness = aRadius - pRadius;
    float stepSize = rayLength / float(NUM_LIGHT_SAMPLES);
    float normalizedStep = stepSize / atmosThickness;
    float depth = 0.0;
    for (int i = 0; i < NUM_LIGHT_SAMPLES; i++) {
        vec3 samplePos = rayOrigin + rayDir * (stepSize * (float(i) + 0.5));
        float altitude = (length(samplePos) - pRadius) / atmosThickness;
        altitude = clamp(altitude, 0.0, 1.0);
        depth += getDensity(altitude) * normalizedStep;
    }
    return depth;
}

void main() {
    scaleHeight = radii.z;
    vec3 rayDir = normalize(view_dir);
    vec3 rayOrigin = v_camera_pos;
    float pRadius = radii.x;
    float aRadius = radii.y;
    float rayleighScale = scatter_coeffs.x;
    float mieScale = scatter_coeffs.y;
    float mieG = scatter_coeffs.z;
    float sunIntensity = scatter_coeffs.w;
    vec3 sunDir = normalize(sun_direction.xyz);

    // Find atmosphere intersection (sphere at origin)
    vec2 atmosHit = raySphereIntersect(rayOrigin, rayDir, aRadius);
    if (atmosHit.y < 0.0) {
        frag_color = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }

    // Ray segment within atmosphere
    float rayStart = max(0.0, atmosHit.x);
    float rayEnd = atmosHit.y;

    // Check planet occlusion (terminate ray at planet surface)
    vec2 planetHit = raySphereIntersect(rayOrigin, rayDir, pRadius);
    if (planetHit.x > 0.0) {
        rayEnd = min(rayEnd, planetHit.x);
    }

    if (rayStart >= rayEnd) {
        frag_color = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }

    float rayLength = rayEnd - rayStart;
    float stepSize = rayLength / float(NUM_SAMPLES);
    float atmosThickness = aRadius - pRadius;
    // Normalize step by atmosphere thickness for scale-independent scattering.
    // Without this, world-unit step sizes (tens of km) cause optical depth to explode.
    float normalizedStep = stepSize / atmosThickness;

    // Dither offset to reduce banding (varies per pixel)
    float dither = hash(gl_FragCoord.xy) * 0.5;

    // Accumulate Rayleigh and Mie scattering
    vec3 rayleighSum = vec3(0.0);
    vec3 mieSum = vec3(0.0);
    float optDepthR = 0.0;
    float optDepthM = 0.0;

    for (int i = 0; i < NUM_SAMPLES; i++) {
        // Sample point along view ray (still in world space for position)
        vec3 samplePos = rayOrigin + rayDir * (rayStart + stepSize * (float(i) + dither));
        float height = length(samplePos);

        // Skip if inside planet body
        if (height < pRadius) continue;

        // Normalized altitude (0 at surface, 1 at atmosphere edge)
        float altitude = (height - pRadius) / atmosThickness;
        altitude = clamp(altitude, 0.0, 1.0);

        // Local atmospheric density
        float localDensity = getDensity(altitude);

        // Accumulate optical depth (normalized units)
        float segmentDepth = localDensity * normalizedStep;
        optDepthR += segmentDepth;
        optDepthM += segmentDepth;

        // Planet shadow test: skip if sun is occluded by planet
        vec2 sunPlanetHit = raySphereIntersect(samplePos, sunDir, pRadius);
        if (sunPlanetHit.x > 0.0) continue;

        // Optical depth from this point toward the sun (through atmosphere)
        vec2 sunAtmosHit = raySphereIntersect(samplePos, sunDir, aRadius);
        float sunRayLength = max(0.0, sunAtmosHit.y);
        float sunOptDepth = computeOpticalDepth(samplePos, sunDir, sunRayLength,
                                                 pRadius, aRadius);

        // Total attenuation: sun -> sample point -> camera
        // Rayleigh: wavelength-dependent (blue scatters more)
        vec3 tauR = rayleighScale * wavelengthsInv4 * (optDepthR + sunOptDepth);
        // Mie: wavelength-independent
        vec3 tauM = vec3(mieScale * (optDepthM + sunOptDepth));
        vec3 attenuation = exp(-(tauR + tauM));

        // Accumulate scattered light (normalized units)
        rayleighSum += localDensity * attenuation * normalizedStep;
        mieSum += localDensity * attenuation * normalizedStep;
    }

    // Apply scattering coefficients
    rayleighSum *= rayleighScale * wavelengthsInv4;
    mieSum *= mieScale;

    // Phase functions (angle between view and sun)
    float cosTheta = dot(rayDir, sunDir);
    float phaseR = rayleighPhase(cosTheta);
    float phaseM = miePhase(cosTheta, mieG);

    // Final scattered light
    vec3 color = sunIntensity * (rayleighSum * phaseR + mieSum * phaseM);

    // Tone mapping (HDR -> LDR)
    color = 1.0 - exp(-color);

    // Alpha based on scattering intensity (controls additive blend amount)
    float alpha = clamp(length(color) * 1.5, 0.0, 0.95);

    frag_color = vec4(color, alpha);
}
@end

@program atmosphere atmosphere_vs atmosphere_fs

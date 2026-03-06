# Moon Reference

## System Constants

| Parameter | Value |
|---|---|
| SOI radius | 50,000 m (50 km) |
| Planet gravity | 10.0 m/s² (Tenebris, ~Earth) |
| Noise amplitude | 0.01 (1% of base_radius, all moons) |
| LOD max depth | 13 |
| Hex terrain range | 500 m |
| Hex voxel size | 1.0 m radius, 1.0 m height |
| Icosphere subdivisions | 4 (~5120 tris) |

## Moon Table

| # | Name | Class | Base Radius | Gravity (m/s²) | Ellipsoid Scale (X,Y,Z) | Axis Ratio | Shape Description |
|---|---|---|---|---|---|---|---|
| 0 | Gorrath | Giant | 80,000 m | 6.0 (60%) | 1.00, 0.90, 1.05 | 1.17:1 | Slightly oblate, wide |
| 1 | Atheron | Giant | 55,000 m | 4.5 (45%) | 1.05, 0.95, 1.00 | 1.11:1 | Nearly spherical |
| 2 | Kelthos | Large | 7,000 m | 2.0 (20%) | 0.95, 1.10, 1.00 | 1.16:1 | Tall/prolate |
| 3 | Dravok | Large | 6,000 m | 1.8 (18%) | 1.20, 0.80, 0.95 | 1.50:1 | **Capsule** (very elongated X) |
| 4 | Serath | Large | 5,000 m | 1.5 (15%) | 1.00, 1.00, 1.15 | 1.15:1 | Stretched Z-axis |
| 5 | Cryx | Small | 3,000 m | 1.0 (10%) | 1.30, 0.70, 1.00 | 1.86:1 | **Capsule** (very elongated X) |
| 6 | Nyctra | Small | 2,000 m | 0.8 (8%) | 0.80, 1.20, 1.10 | 1.50:1 | Tall Y + wide Z |
| 7 | Thalwen | Small | 4,000 m | 1.2 (12%) | 1.10, 0.90, 1.20 | 1.33:1 | Stretched Z-axis |
| 8 | Vexis | Small | 1,500 m | 0.6 (6%) | 1.40, 0.60, 1.00 | 2.33:1 | **Capsule** (most elongated) |
| 9 | Zephyros | Small | 1,000 m | 0.5 (5%) | 1.00, 1.30, 0.80 | 1.63:1 | Tall Y-axis |

*Axis Ratio = max(scale) / min(scale). Higher = more asymmetric.*

## Effective Dimensions (base_radius * scale per axis)

| # | Name | X extent (m) | Y extent (m) | Z extent (m) | Approx Surface Area (km²) |
|---|---|---|---|---|---|
| 0 | Gorrath | 80,000 | 72,000 | 84,000 | ~79,000 |
| 1 | Atheron | 57,750 | 52,250 | 55,000 | ~37,000 |
| 2 | Kelthos | 6,650 | 7,700 | 7,000 | ~600 |
| 3 | Dravok | 7,200 | 4,800 | 5,700 | ~430 |
| 4 | Serath | 5,000 | 5,000 | 5,750 | ~310 |
| 5 | Cryx | 3,900 | 2,100 | 3,000 | ~110 |
| 6 | Nyctra | 1,600 | 2,400 | 2,200 | ~50 |
| 7 | Thalwen | 4,400 | 3,600 | 4,800 | ~200 |
| 8 | Vexis | 2,100 | 900 | 1,500 | ~16 |
| 9 | Zephyros | 1,000 | 1,300 | 800 | ~12 |

## Noise Parameters

| # | Name | Seed | Frequency | Amplitude | Octaves | Max Displacement (m) |
|---|---|---|---|---|---|---|
| 0 | Gorrath | 100 | 0.3 | 0.01 | 4 | ~800 |
| 1 | Atheron | 150 | 0.4 | 0.01 | 4 | ~550 |
| 2 | Kelthos | 300 | 0.5 | 0.01 | 4 | ~70 |
| 3 | Dravok | 400 | 0.6 | 0.01 | 3 | ~60 |
| 4 | Serath | 500 | 0.4 | 0.01 | 3 | ~50 |
| 5 | Cryx | 600 | 0.5 | 0.01 | 2 | ~30 |
| 6 | Nyctra | 700 | 0.6 | 0.01 | 2 | ~20 |
| 7 | Thalwen | 800 | 0.5 | 0.01 | 3 | ~40 |
| 8 | Vexis | 900 | 0.7 | 0.01 | 2 | ~15 |
| 9 | Zephyros | 1000 | 0.7 | 0.01 | 2 | ~10 |

*Max Displacement = base_radius * amplitude = base_radius * 0.01*

## Orbital Parameters

| # | Name | Semi-Major (km) | Period (s) | Period (h:m) | Eccentricity | Inclination (°) |
|---|---|---|---|---|---|---|
| 0 | Gorrath | 5,000 | 14,230 | 3h57 | 0.02 | 5.0 |
| 1 | Atheron | 10,000 | 40,250 | 11h11 | 0.03 | 10.0 |
| 2 | Kelthos | 3,000 | 6,600 | 1h50 | 0.04 | 12.0 |
| 3 | Dravok | 6,500 | 21,100 | 5h51 | 0.06 | 18.0 |
| 4 | Serath | 8,500 | 31,560 | 8h46 | 0.03 | 8.0 |
| 5 | Cryx | 2,000 | 3,600 | 1h00 | 0.08 | 25.0 |
| 6 | Nyctra | 4,500 | 12,150 | 3h22 | 0.05 | 30.0 |
| 7 | Thalwen | 7,500 | 26,200 | 7h17 | 0.07 | 15.0 |
| 8 | Vexis | 11,000 | 46,440 | 12h54 | 0.09 | 40.0 |
| 9 | Zephyros | 13,000 | 59,720 | 16h35 | 0.02 | 10.0 |

## Color Palettes (blue-gray tones matching moon.png)

| # | Name | Base RGB | Highlight RGB | Shadow RGB |
|---|---|---|---|---|
| 0 | Gorrath | 0.38, 0.38, 0.42 | 0.52, 0.52, 0.57 | 0.18, 0.18, 0.22 |
| 1 | Atheron | 0.40, 0.40, 0.44 | 0.55, 0.55, 0.60 | 0.20, 0.20, 0.24 |
| 2 | Kelthos | 0.36, 0.37, 0.40 | 0.50, 0.51, 0.55 | 0.17, 0.17, 0.20 |
| 3 | Dravok | 0.34, 0.34, 0.38 | 0.48, 0.48, 0.53 | 0.16, 0.16, 0.19 |
| 4 | Serath | 0.42, 0.42, 0.46 | 0.57, 0.57, 0.62 | 0.22, 0.22, 0.25 |
| 5 | Cryx | 0.32, 0.32, 0.36 | 0.46, 0.46, 0.51 | 0.15, 0.15, 0.18 |
| 6 | Nyctra | 0.37, 0.38, 0.42 | 0.52, 0.53, 0.58 | 0.18, 0.18, 0.22 |
| 7 | Thalwen | 0.35, 0.35, 0.40 | 0.50, 0.50, 0.56 | 0.17, 0.17, 0.21 |
| 8 | Vexis | 0.33, 0.33, 0.38 | 0.47, 0.47, 0.53 | 0.15, 0.15, 0.19 |
| 9 | Zephyros | 0.39, 0.40, 0.44 | 0.54, 0.55, 0.60 | 0.19, 0.20, 0.23 |

## Gravity / Up Direction

Surface gravity varies per moon (set in `CelestialBody.surface_gravity`). Tenebris = 10.0 m/s² (~Earth).

| # | Name | Gravity | vs Planet | Jump Height | Feel |
|---|---|---|---|---|---|
| - | Tenebris | 10.0 m/s² | 100% | ~3.2 m | Earth-like |
| 0 | Gorrath | 6.0 m/s² | 60% | ~5.3 m | Mars-like |
| 1 | Atheron | 4.5 m/s² | 45% | ~7.1 m | Light, floaty |
| 2 | Kelthos | 2.0 m/s² | 20% | ~16.0 m | Bouncy |
| 3 | Dravok | 1.8 m/s² | 18% | ~17.8 m | Bouncy |
| 4 | Serath | 1.5 m/s² | 15% | ~21.3 m | Very bouncy |
| 5 | Cryx | 1.0 m/s² | 10% | ~32.0 m | Low-grav playground |
| 6 | Nyctra | 0.8 m/s² | 8% | ~40.0 m | Floating |
| 7 | Thalwen | 1.2 m/s² | 12% | ~26.7 m | Low-grav |
| 8 | Vexis | 0.6 m/s² | 6% | ~53.3 m | Near-weightless |
| 9 | Zephyros | 0.5 m/s² | 5% | ~64.0 m | Near-weightless |

*Jump height = JUMP_FORCE² / (2 * gravity) with JUMP_FORCE = 8.0 m/s. Tenebris = 10.0 m/s² (~Earth).*

### Up Direction

Planet uses radial direction (sphere). Moon up direction depends on shape type:

**Ellipsoid** (implemented in `camera.c:332`):
```
radial   = normalize(position)
local_up = normalize(radial.x/sx², radial.y/sy², radial.z/sz²)
```

**Capsule** (TODO — for Dravok, Cryx, Vexis):
```
P        = player position (moon-local)
p_axis   = dot(P, capsule_axis)
clamped  = clamp(p_axis, -half_h, +half_h)
closest  = clamped * capsule_axis

if |p_axis| <= half_h:
    local_up = normalize(P - closest)          // cylinder: perpendicular to axis
else:
    cap_center = sign(p_axis) * half_h * axis
    local_up = normalize(P - cap_center)       // endcap: radial from cap center
```
See [Capsule geometry](#capsule-todo--for-high-axis-ratio-moons) for full derivation.

### Ellipsoid Normal Deviation by Moon

| # | Name | Axis Ratio | Max Normal Deviation | Noticeable? |
|---|---|---|---|---|
| 0 | Gorrath | 1.17:1 | ~5° | Slight |
| 1 | Atheron | 1.11:1 | ~3° | Barely |
| 2 | Kelthos | 1.16:1 | ~5° | Slight |
| 3 | Dravok | 1.50:1 | ~14° | **Yes** |
| 4 | Serath | 1.15:1 | ~4° | Barely |
| 5 | Cryx | 1.86:1 | ~22° | **Very** |
| 6 | Nyctra | 1.50:1 | ~14° | **Yes** |
| 7 | Thalwen | 1.33:1 | ~10° | Moderate |
| 8 | Vexis | 2.33:1 | ~30° | **Extreme** |
| 9 | Zephyros | 1.63:1 | ~17° | **Yes** |

*Max Normal Deviation = approximate max angle between radial direction and ellipsoid normal, at 45° latitude on the most stretched axis.*

## Surface Radius Functions

### Ellipsoid (current, all moons)

```c
float moon_surface_radius(const MoonShapeParams* p, HMM_Vec3 dir) {
    scaled = { dir.X * sx, dir.Y * sy, dir.Z * sz }
    ellipsoid_r = base_radius / length(scaled)
    noise = fbm_opensimplex2(dir * 100.0)
    return ellipsoid_r * (1.0 + noise * amplitude)
}
```

### Capsule (TODO — for high axis ratio moons)

A capsule = cylinder + two hemispherical endcaps. Better than ellipsoid for highly eccentric shapes because:
- Ellipsoids get very "pointy" at the poles when axis ratio > 1.5
- Capsules have uniform curvature at the endcaps — more natural asteroid/moon look
- Flat cylinder sides give interesting terrain to walk on

**Capsule candidates** (axis ratio > 1.5):

| # | Name | Ratio | Capsule Axis | Cap Radius (m) | Half-Height (m) | Cylinder Length (m) |
|---|---|---|---|---|---|---|
| 3 | Dravok | 1.50:1 | X | 5,250 | 1,950 | 3,900 |
| 5 | Cryx | 1.86:1 | X | 2,550 | 1,350 | 2,700 |
| 8 | Vexis | 2.33:1 | X | 1,200 | 900 | 1,800 |

*Cap radius = base_radius * average(two shorter scales). Half-height = base_radius * max_scale - cap_radius.*

**Capsule geometry:**
```
Capsule: axis unit vector `a`, half-height `h` (cylinder half-length), cap radius `r`
         Two cap centers at ±h*a, total extent = 2*(h+r) along axis

Given direction d (unit vector), surface radius from origin:

  da = dot(d, axis)                       // projection onto capsule axis
  threshold = h / sqrt(r² + h²)

  if |da| <= threshold:
      // CYLINDER region — ray hits the flat sides
      surface_r = r / sqrt(1 - da²)

  else:
      // CAP region — ray hits a hemispherical endcap
      surface_r = h * |da| + sqrt(r² - h² * (1 - da²))
```

**Capsule surface normal:**
```
Given point P on capsule surface:
  p_axis = dot(P, axis)                   // project onto axis
  clamped = clamp(p_axis, -h, +h)         // nearest point on cylinder axis
  closest = clamped * axis                 // nearest axis point

  if |p_axis| <= h:
      // CYLINDER region — normal is perpendicular to axis
      normal = normalize(P - closest)

  else:
      // CAP region — normal radiates from cap center
      cap_center = sign(p_axis) * h * axis
      normal = normalize(P - cap_center)
```

**Capsule up direction (gravity):**
```
Same as capsule normal — gravity is perpendicular to the local surface.
On the cylinder sides, "up" points straight out from the axis.
On the endcaps, "up" radiates from the cap center (like a sphere).
```

**Deriving capsule params from MoonShapeParams:**
```c
// Find longest axis
int axis_idx = (sx >= sy && sx >= sz) ? 0 : (sy >= sz) ? 1 : 2;
float max_s = ellipsoid_scale[axis_idx];

// Average of the two shorter scales
float other1 = ellipsoid_scale[(axis_idx+1) % 3];
float other2 = ellipsoid_scale[(axis_idx+2) % 3];
float avg_short = (other1 + other2) * 0.5f;

float cap_radius  = base_radius * avg_short;     // hemisphere radius
float total_extent = base_radius * max_s;         // tip-to-center distance
float half_height  = total_extent - cap_radius;   // cylinder half-length
```

**Implementation plan:**
1. Add `shape_type` field to `MoonShapeParams` (0=ellipsoid, 1=capsule)
2. Add precomputed capsule fields: `capsule_axis[3]`, `capsule_half_h`, `capsule_cap_r`
3. Branch `moon_surface_radius()` for capsule math
4. Branch camera.c `local_up` computation for capsule normal
5. Update `moon_generate_mesh()` to generate capsule icosphere
6. Update hex_terrain and LOD mesh gen to use new surface radius

## Key Code Locations

| System | File | Key Lines |
|---|---|---|
| Moon definitions | src/celestial.c | 448-492 |
| Surface radius | src/celestial.c | 109-135 |
| Icosphere mesh gen | src/celestial.c | 140-230 |
| Orbit update (Kepler) | src/celestial.c | 560-620 |
| Gravity detection | src/camera.c | 520-555 |
| Local up computation | src/camera.c | ~555 |
| Hex collision | src/camera.c | 565-720 |
| Moon gravity path | src/camera.c | 765-810 |
| SOI transitions | src/main.c | 1050-1120 |
| LOD retarget | src/lod.c | lod_tree_retarget() |
| Hex retarget | src/hex_terrain.c | hex_terrain_retarget() |
| Moon hex mesh (LOD) | src/lod.c | 1146-1420 |
| Moon voxel fill | src/hex_terrain.c | 1567-1592 |

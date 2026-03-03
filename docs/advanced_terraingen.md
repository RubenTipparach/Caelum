# Advanced Planet Terrain Generation

A curated collection of techniques, articles, and repositories for making procedural planet terrain more visually stunning — prioritizing noise-based and procedural approaches over expensive simulations.

---

## 1. Noise Patterns & Layering

### Swiss Turbulence
A variant of FBM that multiplies each octave by the previous one, creating terrain with interconnected ridges and valleys reminiscent of Swiss Alps.

```
value = 0, weight = 1
for each octave:
    signal = (1 - abs(noise(p))) * weight
    weight = clamp(signal * gain, 0, 1)
    value += signal * amplitude
    p *= lacunarity; amplitude *= persistence
```

- Creates correlated ridge networks (unlike standard ridged which has independent ridges)
- Weight feedback makes valleys smooth and ridges sharp
- **Key params:** gain (0.5-2.0), persistence (0.5), lacunarity (2.0)
- Reference: [Giliam de Carpentier — "Procedural Swiss Eroded Alps"](https://www.decarpentier.nl/scape-procedural-extensions)

### Ridged Multifractal (Enhanced)
Beyond FastNoiseLite's built-in `FNL_FRACTAL_RIDGED`:

```c
float ridged(vec3 p) {
    float signal = offset - abs(noise(p));
    signal *= signal;  // Sharpen ridges
    float result = signal;
    float weight = 1.0;
    for (int i = 1; i < octaves; i++) {
        p *= lacunarity;
        weight = clamp(signal * gain, 0, 1);
        signal = (offset - abs(noise(p))) * weight;
        signal *= signal;
        result += signal * pow(persistence, i);
    }
    return result;
}
```

- The `offset` parameter (0.9-1.0) controls ridge prominence
- Weight feedback creates dependent octaves — higher octaves only appear on ridges
- Great for mountain chains with realistic sub-ridges

### Domain Warping (Cascaded)
Already in our system at strength 0.3. For more dramatic terrain:

```c
// Single warp (current)
vec3 q = p + warp_noise(p) * strength;
float h = terrain_noise(q);

// Double warp (more organic)
vec3 q = p + noise_A(p) * strength;
vec3 r = p + noise_B(q) * strength;
float h = terrain_noise(r);
```

- Strength 0.5-0.8 for dramatic warping
- Double-cascade creates turbulent, chaotic coastlines and mountain shapes
- Use different seeds for each warp layer
- Reference: [Inigo Quilez — "Domain Warping"](https://iquilezles.org/articles/warp/)

### Multiplicative Noise Layers
Instead of adding octaves, multiply independent noise fields:

```c
float continental = fbm(p * 0.6);         // 0..1 range
float mountain = ridged(p * 1.5);          // 0..1 range
float detail = (fbm(p * 16.0) + 1) * 0.5; // 0..1 range

float height = continental * (0.3 + 0.7 * mountain) * (0.8 + 0.2 * detail);
```

- Multiplication naturally masks features: mountains only appear on continents
- Detail only roughens existing features, doesn't create floating bumps
- Much more realistic than pure additive FBM

### Power/Exponential Redistribution
Simple post-processing that dramatically changes terrain character:

```c
// Sharpen continents (flatten oceans, steepen land)
height = pow(max(height, 0), 1.5);

// Terrace effect
height = round(height * num_terraces) / num_terraces;

// Smooth terraces (Inigo Quilez)
float k = num_terraces;
height = floor(height * k + 0.5) / k + smoothstep(0.45, 0.55, fract(height * k)) / k;
```

---

## 2. Continental & Large-Scale Structure

### Plate Tectonics Simulation (Simplified)
Not full simulation — use Voronoi cells to approximate:

1. Generate N Voronoi cells on the sphere (plate boundaries)
2. Assign each plate a drift direction + speed (random unit vector)
3. At boundaries: convergent → mountain ridge noise, divergent → rift valley
4. Distance from boundary → falloff function
5. Overlay with FBM for local detail

- Creates geologically plausible continental shapes
- Much cheaper than real simulation
- Reference: [Primer on Plate Tectonics for Procedural Generation](https://undiscoveredworlds.blogspot.com/)

### Continent Masks via Threshold + Hysteresis
```c
float continent_noise = fbm(p * 0.5, 3 octaves);
float coast = smoothstep(-0.05, 0.05, continent_noise - sea_threshold);
// coast is 0 in ocean, 1 on land, smooth transition at shoreline
```

- Threshold creates sharp land/water distinction
- Smoothstep width controls coastline sharpness
- Multiply all land features by `coast` to prevent underwater mountains

---

## 3. Biome Systems

### Whittaker Diagram (Temperature + Moisture)
Two-axis biome selection based on real-world ecology:

```c
float temperature = base_temp * cos(latitude) - altitude_lapse * elevation;
float moisture = moisture_noise(p) + ocean_proximity_bonus;

// Biome lookup from 2D gradient
BiomeType biome = whittaker_lookup(temperature, moisture);
```

| | Arid | Moderate | Humid |
|---|---|---|---|
| **Hot** | Desert | Savanna | Tropical Rainforest |
| **Warm** | Scrubland | Grassland | Temperate Forest |
| **Cool** | Steppe | Taiga | Boreal Forest |
| **Cold** | Polar Desert | Tundra | Ice Sheet |

- Latitude → base temperature (equator hot, poles cold)
- Altitude lapse rate: -6.5°C per 1000m (real atmosphere)
- Ocean proximity adds moisture (distance-from-coast falloff)
- Noise perturbation on both axes prevents grid-aligned boundaries

### Slope-Based Material Selection
Cheap and highly effective visual improvement:

```c
float slope = 1.0 - dot(surface_normal, unit_up);  // 0=flat, 1=cliff

if (slope > 0.7) material = ROCK;        // Steep cliffs → exposed rock
else if (slope > 0.4) material = DIRT;    // Moderate slopes → dirt/gravel
else material = biome_ground_material;     // Flat → grass/sand/snow per biome
```

- Automatically creates realistic cliff faces and erosion-like patterns
- Zero additional noise sampling — uses existing surface normals
- Combine with height: high + steep = snow-capped rock, low + steep = mud

### Elevation Band Blending (Enhanced)
```c
// Noise-perturbed height thresholds
float threshold_noise = fbm(p * 20.0) * 100.0;  // ±100m variation
float effective_height = height + threshold_noise;

// Smooth 3-way blend between adjacent biomes
float blend = smoothstep(low_threshold, high_threshold, effective_height);
color = mix(biome_below_color, biome_above_color, blend);
```

---

## 4. Detail Enrichment (Without Erosion Simulation)

### Analytical Erosion Patterns
Fake erosion using noise, not particle simulation:

```c
// Erosion-like valleys using abs(noise)
float erosion = 1.0 - abs(fbm(p * 8.0));
erosion = pow(erosion, 3.0);  // Sharpen into channels

// Apply only on slopes
float slope_factor = saturate(slope * 2.0);
height -= erosion * slope_factor * erosion_depth;
```

- Creates convincing river-like channels in valleys
- Runs at noise sampling speed (not particle simulation)
- Combine with domain warping for meandering channels

### Crater Generation
For moons or rocky planets:

```c
float crater(vec3 p, vec3 center, float radius) {
    float d = distance(p, center) / radius;
    if (d > 1.5) return 0;

    // Bowl shape
    float bowl = smoothstep(0, 1, d) * depth;

    // Raised rim
    float rim = exp(-pow((d - 1.0) * 4.0, 2)) * rim_height;

    // Central peak (for large craters)
    float peak = exp(-d * d * 16.0) * peak_height;

    return -bowl + rim + peak;
}
```

- Poisson-disc distribute crater centers on sphere
- Vary radius with power-law distribution (many small, few large)
- Older craters partially eroded (reduce depth, smooth rim)

### River Carving (Noise-Based)
```c
// Use domain-warped 1D noise to create river paths
float river_noise = abs(warped_noise(p * river_frequency));
float river = 1.0 - smoothstep(0, river_width, river_noise);
height -= river * river_depth;
```

- Creates meandering channels without flow simulation
- Place rivers in valleys (multiply by inverted height gradient)
- Width increases downstream (scale by accumulated flow proxy)

---

## 5. Coloring & Visual Polish

### Tri-Planar Texturing
For steep terrain where UV projection fails:

```c
// Blend 3 axis-aligned texture lookups based on normal
vec3 blend = abs(normal);
blend = normalize(max(blend, 0.001));
blend = pow(blend, vec3(sharpness));  // sharpness=4.0 for crisp transitions
blend /= (blend.x + blend.y + blend.z);

color = texture(tex, p.yz) * blend.x
      + texture(tex, p.xz) * blend.y
      + texture(tex, p.xy) * blend.z;
```

### Height-Based Color Gradients
Instead of flat biome colors, use gradients within each biome:

```c
// Within grassland biome:
vec3 grass_color = mix(
    vec3(0.2, 0.35, 0.05),  // Dark grass (valleys)
    vec3(0.4, 0.55, 0.1),   // Light grass (hilltops)
    saturate(local_height_normalized + color_noise * 0.2)
);
```

- Local height variation creates visual depth
- Noise perturbation prevents banding
- Each biome gets its own gradient (desert: tan→orange, tundra: gray→white)

### Ambient Occlusion from Voxels
For hex terrain at close range:

```c
// Sample 6 neighbors + above for occlusion
float ao = 0;
for each neighbor:
    if (solid) ao += neighbor_weight;
ao = 1.0 - ao * ao_strength;
color *= ao;
```

- Darkens concave corners and crevices
- Huge visual improvement for minimal cost
- Bake into vertex color during mesh generation

---

## 6. Notable Open-Source Projects & References

### Repositories

| Project | Language | Key Technique | Link |
|---------|----------|---------------|------|
| **libnoise** | C++ | Classic procedural noise library, module-based | [github.com/qknight/libnoise](https://github.com/qknight/libnoise) |
| **FastNoiseLite** | C/C++/C# | Single-header noise (already using) | [github.com/Auburn/FastNoiseLite](https://github.com/Auburn/FastNoiseLite) |
| **Accidental Noise Library** | C++ | Hierarchical noise expressions | [github.com/JTippetts/accern](https://github.com/JTippetts) |
| **Open-Simplex-2** | Multi | Improved simplex noise, no patent issues | [github.com/KdotJPG/OpenSimplex2](https://github.com/KdotJPG/OpenSimplex2) |
| **TerraForge3D** | C++ | Node-based terrain gen with GPU acceleration | [github.com/Jaysmito101/TerraForge3D](https://github.com/Jaysmito101/TerraForge3D) |
| **World Machine** | Commercial | Industry-standard terrain tool (study their approach) | [world-machine.com](https://www.world-machine.com/) |
| **Erosion** | C/C++ | Nick McDonald's particle erosion implementation | [github.com/weigert/SimpleHydrology](https://github.com/weigert/SimpleHydrology) |
| **SoilMachine** | C++ | Nick McDonald's full erosion + thermal + wind | [github.com/weigert/SoilMachine](https://github.com/weigert/SoilMachine) |
| **ProcGen Planet** | Various | Collection of procedural planet techniques | Search "procedural planet generation github" |

### Articles & Tutorials

| Title | Author | Key Takeaway |
|-------|--------|--------------|
| [Procedural Swiss Eroded Alps](https://www.decarpentier.nl/scape-procedural-extensions) | Giliam de Carpentier | Swiss turbulence, noise extensions for realistic mountains |
| [Domain Warping](https://iquilezles.org/articles/warp/) | Inigo Quilez | Cascaded warping for organic shapes |
| [Terrain from Noise](https://iquilezles.org/articles/morenoise/) | Inigo Quilez | Analytical derivatives for lighting/erosion |
| [Making Maps with Noise](https://www.redblobgames.com/maps/terrain-from-noise/) | Red Blob Games | Noise redistribution, island masks, biomes |
| [Simulating Hydraulic Erosion](https://jobtalle.com/simulating_hydraulic_erosion.html) | Job Talle | Particle erosion with clear code examples |
| [GPU Procedural Terrain (GDC)](https://www.gdcvault.com/) | Various | No Man's Sky & Houdini terrain talks |
| [Procedural Generation Wiki](http://pcg.wikidot.com/) | Community | Comprehensive collection of PCG techniques |
| [The Book of Shaders: FBM](https://thebookofshaders.com/13/) | Patricio Gonzalez Vivo | Interactive noise + FBM visualizations |
| [GPU Gems Ch.1: Procedural Terrain](https://developer.nvidia.com/gpugems/gpugems3/part-i-geometry/chapter-1-generating-complex-procedural-terrains-using-gpu) | NVIDIA | GPU-accelerated terrain with Perlin noise |

### Video Series

- **Sebastian Lague — Procedural Planets** (YouTube): Excellent walkthrough of noise layering, continent masks, biomes on a sphere
- **Sebastian Lague — Hydraulic Erosion** (YouTube): Clear visual explanation of particle-based erosion
- **Coding Adventure: Procedural Moons and Planets** (Sebastian Lague): Crater generation, atmosphere rendering

---

## 7. Recommended Improvements for Hex-Planets

Ordered by visual impact vs implementation effort:

### Tier 1: High Impact, Low Effort
1. **Slope-based material selection** — Use existing surface normals to pick rock on cliffs, soil on moderate slopes, biome ground on flat. Zero new noise samples needed.
2. **Multiplicative noise layering** — Multiply continental × mountain × detail instead of adding. Mountains naturally only on continents, detail only roughens existing features.
3. **Power redistribution** — `pow(height, 1.3-1.8)` flattens lowlands and steepens mountains. One line of code, dramatic visual change.
4. **Height-based color gradients** — Replace flat biome colors with per-biome gradients (dark valleys, light hilltops). Adds enormous depth.
5. **Vertex AO** — During hex mesh generation, sample neighbor voxels for ambient occlusion. Darkens crevices, adds 3D pop.

### Tier 2: High Impact, Moderate Effort
6. **Swiss turbulence** — Replace or blend with current ridged noise for interconnected mountain ridges. ~20 lines in noise sampling.
7. **Cascaded domain warping** — Increase warp strength to 0.5+ and add second cascade. More organic coastlines and mountain shapes.
8. **Whittaker biomes** — Temperature (latitude + altitude) × moisture (noise + ocean distance) for scientifically plausible biomes. Replace height-only bands.
9. **Analytical erosion channels** — Noise-based fake river channels on slopes. Looks like erosion, runs at noise speed.
10. **Continent masking** — Sharp land/ocean threshold with smooth coastline. All terrain features masked to land areas.

### Tier 3: Moderate Impact, Higher Effort
11. **Tri-planar texturing** — Proper texturing on vertical cliff faces (current UV projection stretches on steep surfaces).
12. **River noise carving** — Domain-warped 1D noise creates meandering river channels. Place in valleys using height gradient.
13. **Crater generation** — For moon/rocky planet variants. Poisson-distributed bowls with rims and central peaks.
14. **Voronoi plate boundaries** — Fake tectonic boundaries with mountain ridges at convergent zones, rifts at divergent zones.

### Current System Strengths (Keep These)
- Continental + ridged mountain + warp + detail noise stack is solid
- FastNoiseLite CPU sampling is fast and proven
- Height quantization to voxel layers works well with hex grid
- Domain warping already creates organic shapes

### What to Avoid
- **Full hydraulic erosion simulation** — Already tried, too expensive for real-time generation. Use analytical erosion patterns instead.
- **Thermal erosion** — Same problem as hydraulic. Slope-based material selection achieves similar visual results.
- **Compute shaders for terrain** — CPU-first approach is working well. GPU terrain adds complexity without proportional gain at current scale.
- **Marching cubes / dual contouring** — Hex voxels are the project's identity. These are alternatives, not improvements.

# Terrain Improvement Roadmap — Progress Tracker

## Phase 1: Improved Noise Terrain ✅ COMPLETE

### 1A. Continental-Scale Structure ✅
- Replaced single base FBM (freq=2.0, 4oct) with `create_continental_noise()` (freq=0.6, FBM 3oct)
- Low-frequency noise creates distinct landmass vs ocean basin shapes
- Domain warping strength increased 0.3 → 0.5 for more organic coastlines
- **File:** `src/lod.c` — `create_continental_noise()`, rewritten `sample_terrain_noise()`

### 1B. Ridged Multifractal Mountains ✅
- Added `create_mountain_noise()` (freq=1.5, RIDGED 5oct)
- Mountain values remapped [−1,1]→[0,1] then squared for sharp peaks
- Mountain mask via `smoothstepf(-0.05, 0.35, continent)` — mountains only on land interior
- Layered: `continent * 0.55 + mountain * 0.45 + detail * weight`
- **File:** `src/lod.c` — `create_mountain_noise()`, updated `sample_terrain_noise()`

### 1C. Biome Blending ✅
- Replaced hard if/else biome boundaries with smooth `vec3_lerp` + `smoothstepf` blending
- ±100m blend zones at each biome threshold (200m, 1500m, 3000m, 4500m above sea level)
- Ocean uses smooth deep→shallow gradient via `vec3_lerp`
- **File:** `src/lod.c` — rewritten `terrain_color_m()`

### Current Noise Stack (after Phase 1)
| Layer | Type | Frequency | Octaves | Purpose |
|-------|------|-----------|---------|---------|
| Continental | OpenSimplex2 FBM | 0.6 | 3 | Landmass/ocean shapes |
| Mountain | OpenSimplex2 RIDGED | 1.5 | 5 | Ridged mountain chains |
| Warp | OpenSimplex2 FBM | 4.0 | 3 | Domain warping (organic coastlines) |
| Detail | OpenSimplex2 RIDGED | 16.0 | 3 | Local surface variation |
| Color | OpenSimplex2 FBM | 40.0 | 2 | Visual color perturbation ±12% |

### Tunable Parameters
- Continental freq (0.6) — lower = fewer/bigger continents
- Mountain freq (1.5) — lower = broader ranges, higher = more jagged
- Warp strength (0.5) — higher = more distorted coastlines
- Layer weights: continent×0.55, mountain×0.45 — relative contribution
- Mountain mask ramp: smoothstep(−0.05, 0.35) — how far inland mountains start
- Biome blend width: 100m half-width at each boundary

---

## Phase 2: Dual Contouring at Depths 10–12 ✅ COMPLETE (core)

**Goal:** Replace smooth tessellation at depths 10–12 with isosurface extraction (voxelized terrain at distance).

### Completed
- [x] Density function: `density = altitude - max(terrain_height, sea_level)`
- [x] QEF solver (3×3 Cramer's rule with mass-point fallback, cell-clamped)
- [x] Voxel grid per patch with adaptive altitude range (pre-sampled min/max + margin)
- [x] Mesh extraction: sign-change edges → quads → triangles (interior edges only)
- [x] Integration: decision gate in `mesh_gen_worker()` — depth 10–12 → DC
- [x] DC meshes output LodVertex (pos+normal+color) — same pipeline as coarse meshes
- [x] Normal computation: central-difference gradients on density grid

### Not Yet Done
- [ ] Transvoxel seam stitching at DC↔coarse boundaries (using boundary skirts for now)
- [ ] Performance tuning (grid resolution, noise caching)

### Files
- `src/dual_contour.h` — DCGrid, DCVertex, DCMesh types + API
- `src/dual_contour.c` — QEF solver, grid gradient, mesh extraction algorithm
- `src/lod.c` — `generate_dc_mesh_for_triangle()`, mesh worker routing

### Grid Resolution (actual)
| Depth | Grid XZ | Grid Y | Notes |
|-------|---------|--------|-------|
| 10 | 32 | 12 | Coarsest DC |
| 11 | 48 | 16 | Mid-range |
| 12 | 64 | 20 | Finest DC (before hex at 13) |

### Reference
- See `memory/terrain-research.md` for DC algorithm details, QEF math, Transvoxel

---

## Phase 3: Hydraulic Erosion ⬜ NOT STARTED

**Goal:** Simulate water flow for realistic river valleys, alluvial fans, eroded ridges.

### Key Tasks
- [ ] Particle-based erosion algorithm (drop→flow→erode→deposit→evaporate)
- [ ] Run in tangent-plane space (reuse hex frame math)
- [ ] Per-patch erosion with boundary overlap to prevent seams
- [ ] Only for depth 10–13 (invisible at coarser depths)
- [ ] Runs on job system workers (already multithreaded)

### New Files Needed
- `src/erosion.h`, `src/erosion.c`

### Key Parameters
- Iterations: 50K, Inertia: 0.3, Capacity: 8.0
- Erosion rate: 0.7, Deposition: 0.02, Evaporation: 0.02, Brush radius: 3

---

## Phase 4: Temperature/Moisture Biome System ⬜ NOT STARTED

**Goal:** Replace height-only biomes with latitude + altitude + moisture driven biome map.

### Key Tasks
- [ ] Temperature = f(latitude, altitude): `base_temp * cos(latitude) - altitude_factor * elevation`
- [ ] Moisture = low-frequency noise + ocean proximity modifier
- [ ] Whittaker-style biome grid (4 temp × 3 moisture = 12 biomes)
- [ ] New `compute_biome()` function replacing `terrain_color_m()` and `lod_hex_atlas()`
- [ ] New texture atlas tiles for additional biome types
- [ ] Smooth transitions via Voronoi noise or smoothstep blending

### Biome Grid
| | Dry | Moderate | Wet |
|---|---|---|---|
| Hot | Desert | Savanna | Tropical Forest |
| Temperate | Steppe | Grassland | Temperate Forest |
| Cold | Tundra | Taiga | Snow Forest |
| Freezing | Ice Sheet | Snow | Glacier |

---

## Phase 5: Polish & Advanced Features ⬜ NOT STARTED

- [ ] Geomorphing (smooth LOD transitions by blending vertex positions)
- [ ] Overhangs/Caves (3D noise in density function — DC supports naturally)
- [ ] GPU compute terrain gen (port noise to HLSL compute shaders)
- [ ] Virtual texturing (decouple texture resolution from atlas)

---

## Implementation Order & Estimates
| Step | Phase | Status |
|------|-------|--------|
| 1 | 1A: Continental noise | ✅ Done |
| 2 | 1B: Ridged mountains | ✅ Done |
| 3 | 1C: Biome blending | ✅ Done |
| 4 | 2: Dual contouring core | ✅ Done |
| 5 | 2: Transvoxel seams | ⬜ (using skirts for now) |
| 6 | 3: Hydraulic erosion | ⬜ |
| 7 | 4: Biome system | ⬜ |
| 8 | 5: Polish features | ⬜ |

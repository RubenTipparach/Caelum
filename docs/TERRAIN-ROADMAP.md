# Terrain Improvement Roadmap — Progress Tracker

## Phase 1: Improved Noise Terrain ✅ COMPLETE

### 1A. Continental-Scale Structure ✅
- Replaced single base FBM (freq=2.0, 4oct) with `create_continental_noise()` (freq=0.6, FBM 3oct)
- Low-frequency noise creates distinct landmass vs ocean basin shapes
- Domain warping strength increased 0.3 → 0.5 for more organic coastlines
- **Files:** `src/lod.c`, `src/hex_terrain.c` — `ht_create_continental_noise()`, `ht_sample_terrain_noise()`

### 1B. Ridged Multifractal Mountains ✅
- Added `create_mountain_noise()` (freq=1.5, RIDGED 5oct)
- Mountain values remapped [−1,1]→[0,1] then squared for sharp peaks
- Mountain mask via `smoothstepf(-0.05, 0.35, continent)` — mountains only on land interior
- Layered: `continent * 0.55 + mountain * 0.45 + detail * weight`
- **Files:** `src/lod.c`, `src/hex_terrain.c` — `ht_create_mountain_noise()`

### 1C. Biome Blending & Dithering ✅
- Smooth `vec3_lerp` + `smoothstepf` blending for vertex colors (LOD rendering)
- Deterministic per-column dithering for voxel types at biome boundaries (`dither_hash()`)
- ±100m dither zones at each biome threshold (sand/grass@200m, grass/stone@1500m, stone/ice@4500m)
- **Files:** `src/lod.c` — `terrain_color_m()`, `src/hex_terrain.c` — `ht_voxel_type()` with dithering

### Current Noise Stack
| Layer | Type | Frequency | Octaves | Seed Offset | Purpose |
|-------|------|-----------|---------|-------------|---------|
| Continental | OpenSimplex2 FBM | 0.6 | 3 | +0 | Landmass/ocean shapes |
| Mountain | OpenSimplex2 RIDGED | 1.5 | 5 | +4000 | Ridged mountain chains |
| Warp | OpenSimplex2 FBM | 4.0 | 3 | +1000 | Domain warping (organic coastlines) |
| Detail | OpenSimplex2 RIDGED | 16.0 | 3 | +2000 | Local surface variation |
| Color | OpenSimplex2 FBM | 40.0 | 2 | — | Visual color perturbation ±12% (LOD only) |

### Tunable Parameters
- Continental freq (0.6) — lower = fewer/bigger continents
- Mountain freq (1.5) — lower = broader ranges, higher = more jagged
- Warp strength (0.5) — higher = more distorted coastlines
- Layer weights: continent×0.55, mountain×0.45 — relative contribution
- Mountain mask ramp: smoothstep(−0.05, 0.35) — how far inland mountains start
- Biome blend/dither width: 100m half-width at each boundary

---

## Phase 2: Voxel Hex Terrain System ✅ COMPLETE

**Goal:** Full 3D voxel terrain at close range with hex-prism geometry, block interactions, lighting, and physics.

*(Dual Contouring was prototyped but removed — replaced entirely by the hex voxel system.)*

### Terrain Generation
- [x] 4-layer noise stack (continental + mountain + warp + detail) shared with LOD system
- [x] `ht_sample_height_m()` — single source of truth for terrain height at any unit-sphere position
- [x] 3D voxel fill: surface type by height + depth layering (surface→dirt→stone)
- [x] Deterministic biome dithering at type transitions (`dither_hash()`)
- [x] Stone arch demo structure at (0,1,0) on unit sphere
- [x] Bedrock floor (3 unbreakable layers at bottom of each chunk)

### Chunk System
- [x] 32×32 hex grid per chunk, 256 vertical layers (HEX_CHUNK_LAYERS)
- [x] 128 layers of excavation padding below terrain (HEX_BASE_PADDING)
- [x] Up to 512 active chunks (HEX_MAX_CHUNKS), range-based activation/deactivation
- [x] Per-column min/max solid cache for fast culling
- [x] Threaded mesh generation via job system (up to 16 jobs/frame)
- [x] Second-pass job submission for dirty chunks missed by activation loop
- [x] Watchdog: force-dirty stuck chunks with no mesh

### Mesh Generation
- [x] Full hex-prism geometry: top face, bottom face, 6 side faces per exposed voxel
- [x] Textured via 7-tile atlas (water, sand, dirt, grass, stone, ice, snow) from 16×16 PNGs
- [x] Per-vertex sky light baked during mesh gen (BFS-computed, 0–32 range)
- [x] Slope normals from noise gradients (matches LOD shading)
- [x] Physics wireframe overlay (P key toggle, green hex-prism outlines)
- [x] Transition mode disabled — all chunks render full voxel prisms

### Sky Light (Minecraft-style BFS)
- [x] Pass 1: sky column fill (top-down scan, air above topmost solid → SKY_MAX=32)
- [x] Pass 2: BFS seeding from boundary blocks + cross-chunk height estimation
- [x] Pass 3: flood-fill with sky column rule (downward@MAX stays MAX, else −1 per step)
- [x] Direct sky_map reads for face lighting (no smoothed averaging)
- [x] Queue bounded at HEX_CHUNK_SIZE³ with overflow protection

### Block Interactions
- [x] Raycast: 0.1m step along camera forward, tangent-plane hex coord lookup
- [x] Face detection: top(0), bottom(1), side(2+) from air→solid transition direction
- [x] Break: set voxel to AIR, update column cache, mark chunk + neighbors dirty
- [x] Place: set voxel at placement position, validate chunk bounds
- [x] Bedrock protection: `VOXEL_BEDROCK` cannot be broken
- [x] Ctrl mode: `sides_only` raycast skips top/bottom for side-only placement

### Selection Highlight
- [x] White wireframe prism outline (36 line vertices, slightly oversized)
- [x] Green filled placement face (triangles, alpha 0.35)
- [x] Face types: hex fan for top/bottom (12 verts), quad for sides (6 verts)

### Configuration
| Define | Value | Description |
|--------|-------|-------------|
| HEX_RADIUS | 1.0m | Hex circumradius |
| HEX_HEIGHT | 1.0m | Vertical layer height |
| HEX_RANGE | 500m | Activation range from camera |
| HEX_CHUNK_SIZE | 32 | Hex columns per chunk side |
| HEX_CHUNK_LAYERS | 256 | Vertical layers per chunk |
| HEX_BASE_PADDING | 128 | Excavation room below terrain |
| HEX_BEDROCK_LAYERS | 3 | Unbreakable bottom layers |
| HEX_MAX_CHUNKS | 512 | Maximum active chunks |
| HEX_MAX_DRAW_ALT | 2000m | Hex terrain disabled above this |
| SKY_MAX | 32 | Maximum sky light level |

### Files
- `src/hex_terrain.h` — Config defines, HexChunk/HexTerrain/HexHitResult types, public API
- `src/hex_terrain.c` — Full implementation (~2300 lines): terrain gen, mesh building, BFS sky light, raycast, break/place
- `src/hex_vertex.h` — HexVertex format (pos, normal, uv, color, sky_light = 44 bytes)
- `shaders/hex_terrain.glsl` — Textured shader with sky light, sun angle, fog
- `shaders/highlight.glsl` — Wireframe/face highlight shader
- `assets/textures/` — 16×16 terrain tile PNGs (7 types)

---

## Phase 3: Hydraulic Erosion ⬜ NOT STARTED

**Goal:** Simulate water flow for realistic river valleys, alluvial fans, eroded ridges. Cache results on disk so erosion only runs once per region.

### Key Tasks
- [ ] Particle-based erosion algorithm (drop→flow→erode→deposit→evaporate)
- [ ] Heightmap-based: erode a 2D height grid, then use eroded heights during voxel fill
- [ ] Run in tangent-plane space (reuse hex terrain tangent frame math)
- [ ] Per-region erosion with boundary overlap to prevent seams
- [ ] Runs on job system workers (already multithreaded)
- [ ] Disk cache: store eroded heightmaps in SQLite after first generation
- [ ] Cache key: (region_coords, seed, erosion_params_hash) — invalidate on param change
- [ ] Applies to both LOD system (color/normals) and hex terrain (voxel heights)

### Disk Cache Strategy
- SQLite database (`world.db` or `terrain_cache.db`)
- Table: `erosion_cache(region_x, region_z, seed, params_hash, heightmap BLOB, timestamp)`
- Heightmap stored as compressed float array (zlib or LZ4)
- On chunk gen: check cache → hit = load heightmap, miss = run erosion + store
- Cache warming: background thread pre-erodes visible regions during idle frames

### New Files Needed
- `src/erosion.h`, `src/erosion.c` — Particle erosion algorithm + cache API

### Key Parameters
- Iterations: 50K per region, Inertia: 0.3, Capacity: 8.0
- Erosion rate: 0.7, Deposition: 0.02, Evaporation: 0.02, Brush radius: 3
- Region size: ~512×512 samples (covers multiple hex chunks with overlap)

---

## Phase 4: Temperature/Moisture Biome System ⬜ NOT STARTED

**Goal:** Replace height-only biomes with latitude + altitude + moisture driven biome map.

### Key Tasks
- [ ] Temperature = f(latitude, altitude): `base_temp * cos(latitude) - altitude_factor * elevation`
- [ ] Moisture = low-frequency noise + ocean proximity modifier
- [ ] Whittaker-style biome grid (4 temp × 3 moisture = 12 biomes)
- [ ] New `compute_biome()` function replacing `terrain_color_m()` and `ht_voxel_type()`
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
- [ ] Overhangs/Caves (3D noise carving in voxel fill phase)
- [ ] GPU compute terrain gen (port noise to HLSL compute shaders)
- [ ] Virtual texturing (decouple texture resolution from atlas)
- [ ] World persistence (SQLite: player edits, chunk modifications)

---

## Implementation Order & Estimates
| Step | Phase | Status |
|------|-------|--------|
| 1 | 1A: Continental noise | ✅ Done |
| 2 | 1B: Ridged mountains | ✅ Done |
| 3 | 1C: Biome blending + dithering | ✅ Done |
| 4 | 2: Voxel hex terrain system | ✅ Done |
| 5 | 2: Sky light BFS | ✅ Done |
| 6 | 2: Block interactions + highlight | ✅ Done |
| 7 | 3: Hydraulic erosion + disk cache | ⬜ |
| 8 | 4: Biome system | ⬜ |
| 9 | 5: Polish features | ⬜ |

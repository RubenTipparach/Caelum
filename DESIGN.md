# Hex-Voxel Planetary Engine -- Design Document

> A hexagon-based voxel game engine capable of generating massive planets in a
> full solar system, with sophisticated LOD, dynamic chunk streaming, voxel
> lighting, orbital mechanics, and realistic day/night cycles.

---

## Table of Contents

1. [Tech Stack Decision](#1-tech-stack-decision)
2. [Hex Sphere Generation](#2-hex-sphere-generation)
3. [Distortion Minimization](#3-distortion-minimization)
4. [LOD System](#4-lod-system)
5. [Chunking & Streaming](#5-chunking--streaming)
6. [Voxel Lighting](#6-voxel-lighting)
7. [Collision & Movement](#7-collision--movement)
8. [Orbital Mechanics & Day/Night](#8-orbital-mechanics--daynight)
9. [Reference Material](#9-reference-material)
10. [Multiplayer & Networking](#10-multiplayer--networking)
11. [Persistence & Data Storage](#11-persistence--data-storage)

---

## 1. Tech Stack Decision

### Hard Requirements

| Requirement | Rationale |
|---|---|
| **Multithreading** (non-negotiable) | Chunk meshing, noise generation, lighting propagation -- all must run on worker threads. A single-threaded engine cannot generate planet-scale hex terrain in real time. |
| **Mature graphics library** (raylib-like) | We want a batteries-included framework with simple APIs, not raw Vulkan/D3D12 boilerplate. Must handle windowing, input, audio, basic 3D rendering, and shaders. |
| **Web browser support** (bonus) | Ship a playable build via WebAssembly + WebGPU/WebGL. Reaches the widest audience with zero install friction. |
| **Custom mesh generation** | We build every triangle ourselves (hex prisms, LOD transitions). Need low-level vertex/index buffer control. |
| **No Vulkan** | We explicitly avoid Vulkan backends. The maintenance overhead and complexity are not worth it for this project. D3D11/Metal/WebGPU cover all target platforms. |

### Detailed Comparison: Sokol vs raylib vs bgfx

These are the three viable C-friendly rendering libraries for this project.
Below is a feature-by-feature comparison to make the decision concrete.

#### Backend Support

| Backend | Sokol | raylib | bgfx |
|---------|-------|--------|------|
| OpenGL / GLES3 | ✅ | ✅ (GL 3.3 is the *only* native backend) | ✅ |
| D3D11 | ✅ | ❌ (only via ANGLE, experimental) | ✅ |
| D3D12 | ❌ | ❌ | ✅ |
| Metal | ✅ | ❌ (only via ANGLE, experimental) | ✅ |
| WebGL2 | ✅ | ✅ (via Emscripten) | ✅ |
| WebGPU | ✅ (first-class) | ❌ | ✅ (Dawn-based, maturing) |

**Note:** Vulkan backends exist in bgfx and experimentally in Sokol, but we
explicitly exclude Vulkan from consideration. D3D11/Metal/WebGPU cover all
target platforms without the maintenance overhead.

#### Rendering Capabilities

| Feature | Sokol | raylib | bgfx |
|---------|-------|--------|------|
| Custom vertex/index buffers | ✅ Full low-level control (`sg_make_buffer`) | ⚠️ Via `Mesh` struct or `rlgl.h` low-level calls | ✅ Full low-level control (static/dynamic/transient buffers) |
| Instancing | ✅ (`sg_draw` with instance count) | ✅ (`DrawMeshInstanced`) | ✅ (dedicated instance data buffers) |
| Multi-Draw Indirect | ❌ | ❌ | ✅ (GPU-driven rendering pipelines) |
| Vertex pulling (storage buffers) | ✅ (explicitly supported with examples) | ❌ | ⚠️ Via compute-writable vertex buffers |
| Draw call batching | Manual (you control draw order) | Internal 2D batch system only | Automatic sort-based bucketing (minimizes state changes) |
| Shader toolchain | `sokol-shdc` (GLSL → cross-platform, offline) | Raw GLSL, loaded at runtime | `shaderc` (GLSL/HLSL → all backends, offline) |

#### Compute Shader Support

| Aspect | Sokol | raylib | bgfx |
|--------|-------|--------|------|
| Compute shaders? | ✅ (added 2025) | ✅ (since v4.0, OpenGL 4.3 flag) | ✅ |
| Available on web? | ✅ (WebGPU backend) | ❌ (GL 4.3 not available on web) | ✅ (WebGPU backend) |
| Available on mobile? | ❌ (GLES3 lacks compute) | ❌ | ✅ (Metal on iOS, GLES 3.1 on Android) |
| Available on desktop? | ✅ (Metal, D3D11, desktop GL) | ✅ (desktop GL 4.3 only) | ✅ (all desktop backends) |

See "Do We Need Compute Shaders?" below for whether this even matters.

#### Multithreading

| Feature | Sokol | raylib | bgfx |
|---------|-------|--------|------|
| Thread safety | ❌ Single-threaded. All `sg_*` calls on one thread. | ❌ Single-threaded. OpenGL context is main-thread only. | ✅ Built-in multithreaded architecture. |
| Render thread separation | ❌ (community wrapper exists) | ❌ | ✅ Dedicated render thread via `bgfx::renderFrame()`. |
| Multi-thread draw submission | ❌ | ❌ | ✅ Encoder API allows up to 8 threads to submit draws simultaneously. |
| Async resource loading | ✅ `sokol_fetch.h` for async file I/O. GPU resource creation on main thread. | ❌ No built-in async I/O. GPU ops on main thread. | ✅ Resource API calls are mutex-guarded and thread-safe from any thread. |

**Key insight:** All three libraries require our own job system for background
chunk meshing/noise/lighting. The difference is draw submission:
bgfx can submit from multiple threads, Sokol and raylib cannot. For our use case
(hundreds of chunk meshes per frame), this matters less than it sounds — the
bottleneck is mesh generation, not draw submission.

#### Built-In Features

| Feature | Sokol | raylib | bgfx |
|---------|-------|--------|------|
| Windowing | ✅ (`sokol_app.h`) | ✅ (built-in) | ❌ (needs SDL/GLFW) |
| Input | ✅ (keyboard, mouse, touch, gamepad) | ✅ (keyboard, mouse, touch, gamepad) | ❌ |
| Audio | ✅ (`sokol_audio.h`) | ✅ (WAV, OGG, MP3, FLAC, QOA) | ❌ |
| Camera | ❌ | ✅ (built-in 2D/3D camera) | ❌ |
| Font/Text | ⚠️ Via `sokol_fontstash.h` | ✅ (TTF, OTF, BDF) | ❌ |
| UI | ⚠️ Via `sokol_imgui.h` / `sokol_nuklear.h` | ✅ `raygui` (separate lib) | ❌ |
| Model loading | ❌ | ✅ (OBJ, GLTF, IQM, VOX) | ❌ |
| Image loading | ⚠️ Via stb_image internally | ✅ (PNG, BMP, JPG, GIF, DDS, HDR, etc.) | ❌ |

#### Performance & Scale

| Metric | Sokol | raylib | bgfx |
|--------|-------|--------|------|
| Max draw calls @ 60fps | ~10K-30K (D3D11/Metal efficient) | ~1K-5K (single-thread GL overhead) | ~64K default (configurable) |
| Triangles per frame | Millions (GPU-bound) | Millions possible but CPU-bound sooner | Millions (designed for throughput) |
| Voxel engine suitability | Good — vertex pulling + instancing for chunks | Workable but CPU bottleneck with many chunks | Best — MDI + sort batching + multi-thread submit |
| WASM binary size | ~33 KB (basic demo) | Larger (full library linked in) | Moderate |
| Desktop binary footprint | Very small (only link what you use) | ~1-2 MB static library | ~1-3 MB (varies by backends) |
| Compile time | Fast (pure C, single headers) | Fast (pure C99) | Slower (C++, multiple backends) |

#### Community & Ecosystem

| Metric | Sokol | raylib | bgfx |
|--------|-------|--------|------|
| GitHub stars | ~9.6K | ~25K+ | ~15.8K |
| Language | C99 (STB-style headers) | C99 | C++ (with auto-generated C API) |
| Bindings | Zig, Odin, Rust, Nim, D, V, Python | 60+ language bindings | C, C#, Rust, Zig, Java, Go, Swift |
| Examples | ~80 samples | 140+ with interactive web demos | 49 examples |
| Documentation | Excellent blog posts, inline header docs | Website, cheatsheet, wiki, large community | Formal docs site, API reference |
| License | zlib/libpng | zlib/libpng | BSD 2-Clause |
| Shipped games | Solar Storm, Medivia Online, various indie | 900+ community projects, popular for jams | Minecraft Bedrock, Rotwood (Klei), Talking Tom franchise |

#### Notable Engines & Shipped Games

**Sokol — Engines & Projects**

| Project | Description |
|---------|-------------|
| Medivia Online | Free-to-play MMORPG (15K+ active players), sponsored sokol_gp development |
| Solar Storm | Turn-based sci-fi artillery game (Odin + Sokol), custom shadowmapping pipeline |
| Soluna | Lua-based 2D game framework using sokol_gfx |
| z_impact | Zig-based 2D action game engine with sokol as platform backend |
| Doom (sokol port) | Classic Doom shareware running in-browser via WASM, by floooh |
| Tiny 8-bit Emulators | ZX Spectrum, C64, CPC emulators by floooh — cross-platform via sokol |
| Outline 2024 demo | Demo by Aras Pranckevičius (former Unity CTO), 4th place Newskool Demo |

**raylib — Engines & Projects**

| Project | Description |
|---------|-------------|
| CAT & ONION | Whimsical cat adventure (Odin + raylib), first Odin game on Steam |
| Lurking I: Immortui | Party-based retro fantasy RPG, rated 4.8/5 |
| MasterPlan | Project management / visual idea board (Go + raylib), available on Steam |
| Waste of Space | Spaceship builder with deterministic rollback netcode co-op |
| rayfork | Single-header C99 game library forked from raylib internals |
| Galaxy Engine | Gravity physics simulator, 100K+ particle rendering |
| Gloamvault | Roguelike monster collector + first-person dungeon crawler |

**bgfx — Engines & Projects**

| Project | Description |
|---------|-------------|
| Minecraft: Bedrock Edition | Best-selling game of all time — RenderDragon is a closed-source bgfx fork |
| Starlite Engine (Outfit7) | Powers entire Talking Tom franchise (billions of downloads) |
| Rotwood (Klei) | Beat-em-up roguelike by the makers of Don't Starve |
| Football Manager 2018 | Sports sim by Sports Interactive / Sega |
| Braid: Anniversary Edition | Jonathan Blow's remastered puzzle-platformer (Android) |
| Unity Project Tiny | Unity's experimental lightweight runtime used bgfx rendering backend |
| Babylon Native (Microsoft) | Cross-platform framework for running Babylon.js natively via bgfx |
| MAME | Multiple Arcade Machine Emulator — bgfx is the primary shader/rendering backend |
| FFNx | Next-gen modding platform for Final Fantasy VII & VIII |
| Smith and Winston | Destructible voxel twin-stick shooter (PC, PS4, Xbox One) |

**Summary**

| Criteria | Sokol | raylib | bgfx |
|----------|-------|--------|------|
| Biggest shipped game | Medivia Online (MMORPG) | CAT & ONION (indie) | Minecraft: Bedrock Edition |
| AAA / major studio use | ❌ | ❌ | ✅ (Mojang, Sega, Klei, Outfit7) |
| Console shipping confirmed | ❌ | ❌ | ✅ (PS4, Xbox, Switch) |
| Mobile at scale | ❌ | ❌ | ✅ (Talking Tom — billions of downloads) |
| Engine adoption | Small frameworks | rayfork, CopyCat | Starlite, Unity Tiny, Babylon Native, Crown, NeoAxis |
| Major tool adoption | Emulators, demos | MasterPlan | MAME, FFNx |

### Ranking

**1st: Sokol** — Best balance of simplicity, modern backends, and web support.
First-class WebGPU gives us the best web target. The API is clean enough for
custom hex mesh generation without fighting abstractions. Vertex pulling via
storage buffers is well-suited for voxel chunk rendering. The 33 KB WASM demos
are unmatched. Modular headers mean we only pull in what we need.

**2nd: bgfx** — Most capable renderer with the best multithreading story and
MDI support. Strongest track record in shipped commercial games. However: it's
C++ (not C), requires a separate windowing library (SDL3), has a more complex
build system, and its API is more idiosyncratic. The extra power is nice but
may not be needed — our bottleneck is chunk generation, not draw submission.

**3rd: raylib** — Best for prototyping and learning, but weakest fit for
production. OpenGL-only backend limits us to WebGL on web (no WebGPU path), the
single-threaded rendering ceiling is lower, and draw call throughput is
significantly behind the other two. Excellent for quickly testing hex grid math
or noise algorithms, but not the final renderer.

### Do We Need Compute Shaders?

Minecraft has no compute shaders. It uses CPU-side multithreaded chunk meshing
and is the best-selling game of all time. Let's evaluate whether we actually
need them.

**Task-by-task analysis:**

| Task | CPU-Only Viable? | Compute Speedup | Complexity Cost | Recommendation |
|------|-----------------|-----------------|-----------------|----------------|
| **Mesh generation** (greedy meshing) | Yes — 50-200 μs per 32³ chunk. 8 threads = 40K-160K chunks/sec. | Modest (sequential merge step is hard to parallelize) | High (indirect draw pipeline) | **CPU** |
| **Noise evaluation** (terrain gen) | Yes — FastNoiseLite SIMD, ~0.5-2 ms per chunk. 8 threads = 50-200 chunks/sec. | 10-100x (embarrassingly parallel) | Moderate | **CPU first, GPU later if needed** |
| **Light propagation** (BFS flood fill) | Yes — Minecraft-proven. Per-chunk BFS is independent. | Moderate (tight data dependencies) | High (jump flooding, sync) | **CPU** |

**Verdict: Start CPU-only with multithreading.** Binary greedy meshing +
FastNoiseLite + BFS lighting on a thread pool. This is the Minecraft-proven
path and handles planetary scale when combined with LOD streaming.

**Design for compute later:** Keep chunk data in GPU-friendly formats. Use
vertex pulling (which Sokol supports) so vertex layout is already storage-buffer
compatible. If we ever hit the CPU bottleneck (fast travel, mass LOD transitions),
GPU noise generation is the first thing to offload — it's embarrassingly parallel
and gives 10-100x speedup.

**This means compute shader support is a nice-to-have, not a blocker.** All three
libraries can do CPU-only voxel rendering. Sokol and bgfx give us the option to
add compute later if needed.

### Recommendation

```
┌─────────────────────────────────────────────────────┐
│              RECOMMENDED STACK                       │
├──────────┬──────────────────────────────────────────┤
│ Language │  C                                       │
│ Renderer │  Sokol (sokol_gfx + sokol_app)           │
│          │  → WebGPU backend for web                │
│          │  → Metal / D3D11 for native              │
│ Threading│  C11 threads + custom job system          │
│ Audio    │  sokol_audio (or miniaudio)              │
│ Math     │  HandmadeMath.h                          │
│ Noise    │  FastNoiseLite (CPU, multithreaded)      │
│ Web      │  Emscripten + WebGPU (native in Sokol)   │
│ Compute  │  CPU-only initially; GPU noise via        │
│          │  Sokol compute if needed later            │
└──────────┴──────────────────────────────────────────┘
```

**Why Sokol:** Cleanest API for custom mesh generation, first-class WebGPU for
web builds, modular headers (only use what you need), tiny footprint, and
compute shaders available when we need them. The draw call ceiling (~10K-30K)
is sufficient — with greedy meshing, 1000 visible chunks × ~10 quads/chunk
average = well within budget.

**Why not bgfx:** More power than we need, C++ dependency, requires SDL3 for
windowing, more complex build. If we hit Sokol's draw call ceiling, bgfx is
the upgrade path.

**Why not raylib for production:** OpenGL-only limits web to WebGL (no WebGPU),
lower draw call throughput, no path to compute shaders on web. Good for
prototyping individual systems though.

### Multithreading Architecture

```
Main Thread          ─── Input, game logic, render submission
Render Thread        ─── sokol_gfx command encoding (automatic in Sokol)
Worker Pool (N cores)─── Chunk meshing, noise generation, lighting BFS,
                         LOD tree traversal
Async I/O Thread     ─── Chunk save/load, streaming
```

**Job System Design:**
- Lock-free work-stealing queue (Chase-Lev deque per worker)
- Jobs are small, self-contained units: `MeshChunk(chunk_id)`,
  `PropagateLight(chunk_id)`, `GenerateNoise(chunk_id, lod_level)`
- Dependencies expressed as job chains: `GenerateNoise → MeshChunk → UploadGPU`
- Double-buffered chunk data: workers write to back buffer, render reads front buffer

---

## 2. Hex Sphere Generation

### The Goldberg Polyhedron Approach

Our planets use **Goldberg polyhedra** -- the mathematical dual of geodesic
spheres. Start with an icosahedron (20 triangular faces), subdivide, and take
the dual to get a sphere tiled with hexagons and exactly **12 pentagons**
(a topological invariant from Euler's formula: V - E + F = 2).

**Construction: GP(m, n)**

The notation GP(m,n) describes the subdivision. Starting from one pentagon,
take `m` steps in one direction on the hex grid, turn 60 degrees left, then
take `n` steps to reach the next pentagon.

- **T-number** = m² + mn + n² (triangulation number)
- **Hex count** = 10T - 10
- **Pentagon count** = always 12
- **Total cells** = 10T + 2

| GP(m,n) | T | Total Cells | Approx Cell Width (Earth-size) | Use Case |
|---------|---|-------------|-------------------------------|----------|
| GP(4,0) | 16 | 162 | ~1,800 km | Solar system overview |
| GP(16,0) | 256 | 2,562 | ~450 km | Continent scale |
| GP(64,0) | 4,096 | 40,962 | ~120 km | Strategy game |
| GP(256,0) | 65,536 | 655,362 | ~30 km | Detailed terrain |
| GP(1024,0) | 1,048,576 | 10,485,762 | ~7 km | Fine terrain |
| GP(4096,0) | 16,777,216 | 167,772,162 | ~1.8 km | Near-surface |
| Finest LOD | — | ~billions | ~1-10 m | Walking scale |

### Construction Algorithm

**Phase 1 -- Build the planar hex grid template:**
1. Lay out a regular hexagonal tiling in 2D (axial coordinates q, r).
2. Define the "Goldberg triangle" from the (m, n) parameters -- an equilateral
   triangle containing T hexagonal cells.
3. This triangle is our fundamental patch.

**Phase 2 -- Map onto the icosahedron:**
1. Take the 20 equilateral triangular faces of an icosahedron.
2. Map the Goldberg triangle onto each face via barycentric coordinates.
3. Stitch adjacent face patches along shared edges.
4. Project all vertex positions outward onto the sphere (normalize to radius R).
5. The 12 original icosahedron vertices become the 12 pentagonal cell centers.

**Phase 3 -- Build the hex-prism voxels:**
1. At each hex cell center on the sphere, define a vertical column extending
   inward (toward planet core) and outward (atmosphere).
2. Each column is a stack of hexagonal prism voxels.
3. Terrain height determines which prisms are solid vs. air.

### Handling the 12 Pentagons

The pentagons are unavoidable. Strategies (in order of preference):

1. **Treat as 5-neighbor hex cells.** All game logic checks `cell.neighbor_count`
   instead of hardcoding 6. Pentagons participate in lighting, meshing, and
   physics identically to hexagons -- they just have one fewer neighbor.

2. **Place strategically.** The 12 pentagons sit at icosahedron vertices.
   Position them in oceans, at poles, or under terrain features where the
   visual difference is minimal.

3. **At high subdivision levels, they're invisible.** At GP(64,0)+, pentagons
   are 0.03% of cells and nearly indistinguishable from hexagons.

### Coordinate System

Use **hierarchical face-local coordinates:**

```
CellID = {
    u8  face;        // 0-19 (icosahedral face)
    u8  lod_level;   // subdivision depth
    i32 q, r;        // axial hex coordinates within the face
    i16 y;           // vertical layer (voxel stack index)
}
```

For algorithms (pathfinding, lighting, distance), use **cube coordinates**
(q, r, s where q + r + s = 0) -- arithmetic is uniform and intuitive.
See Red Blob Games' canonical hex grid reference for the full treatment.

For storage, use **axial coordinates** (q, r) since s = -q - r is redundant.

For indexing, pack into a **64-bit cell ID** (similar to Uber H3):
```
Bits 63-60:  LOD level (0-15)
Bits 59-55:  Icosahedral face (0-19)
Bits 54-32:  q coordinate (signed 23-bit)
Bits 31-16:  r coordinate (signed 16-bit)
Bits 15-0:   y layer (signed 16-bit)
```

---

## 3. Distortion Minimization

### The Problem

A sphere cannot be tiled with perfectly regular hexagons. Any mapping from a
flat hex grid to a sphere introduces distortion in area, shape, or both. We
need to choose a projection that minimizes the distortion most harmful to
gameplay and rendering.

### Projection Comparison

| Projection | Equal Area? | Area Variation | Shape Distortion | Implementation Complexity |
|---|---|---|---|---|
| **ISEA (Snyder Equal-Area)** | Yes | 1.0:1 | Moderate (~11° max) | Medium |
| **Gnomonic (per face)** | No | ~1.3:1 | Low near center | Low |
| **Fuller / Dymaxion** | No | ~1.2:1 | Low (uniform) | High (hard to invert) |
| **S2 (Google)** | No | ~1.4:1 | Low | Medium (but uses quads) |
| **H3 (Uber)** | No | ~1.3:1 | Low-medium | Low (library exists) |

### Recommendation: ISEA with Aperture 4

**Icosahedral Snyder Equal-Area (ISEA)** projection with **aperture 4**
subdivision is the best fit:

- **Equal area:** Every hex cell covers the same surface area. This means
  terrain generation, biome distribution, and resource placement are uniform.
  No region of the planet is more "stretched" than another.

- **Aperture 4:** Each hex subdivides into 4 children (axis-aligned). This
  gives us a clean quadtree-like LOD hierarchy -- easy to reason about,
  efficient to traverse, simple parent-child relationships.

  Alternative: Aperture 7 (each hex → 7 children, rotated ~19.1°) preserves
  hexagonal symmetry better but creates a more complex tree and slight
  rotation artifacts between LOD levels.

- **Max angular distortion ~11°:** Hex cells near icosahedral face edges are
  slightly elongated. At gameplay scales this is imperceptible.

**How ISEA works:**
1. Project the sphere onto the 20 faces of an inscribed icosahedron.
2. Within each face, apply Snyder's modified Lambert azimuthal equal-area
   transform -- this maps the spherical triangle to a flat equilateral
   triangle while preserving area ratios.
3. Tile the flat triangle with a regular hex grid.
4. Inverse-project hex centers back to the sphere.

Reference implementations: DGGRID by Kevin Sahr, geogrid library, and the
PROJ library's `isea` projection.

---

## 4. LOD System

### The Scale Problem

Earth-radius planet (~6,371 km) viewed from surface to orbit:
- Detail range: **~1:10,000,000** (1m hexes to full-planet view)
- At 1m resolution: **~510 trillion** surface hex cells
- With depth layers: multiply by underground depth
- Obviously cannot fit in memory → **LOD is mandatory**

### Architecture: Hybrid Static/Dynamic LOD

The LOD system has two modes that blend based on camera distance:

```
Distance from surface:   Rendering Mode:
─────────────────────   ──────────────
> 100 km                STATIC: Pre-baked impostor spheres,
                        low-poly icosphere with baked textures
                        (like KSP's "Scaled Space")

10 km - 100 km          STATIC: Pre-generated chunked meshes at
                        fixed LOD levels, streamed from disk/memory

0.5 km - 10 km          DYNAMIC: GPU-generated hex meshes at
                        variable LOD, geomorphing transitions

0 - 0.5 km              DYNAMIC: Full-resolution hex prism voxels,
                        breakable/buildable, collision-active
```

### LOD Tree Structure

The LOD hierarchy follows the icosahedral subdivision:

```
Level 0:   20 icosahedral faces (top-level chunks)
Level 1:   20 × 4 = 80 chunks (aperture 4 subdivision)
Level 2:   80 × 4 = 320 chunks
Level 3:   320 × 4 = 1,280 chunks
...
Level N:   20 × 4^N chunks
```

Each level doubles the linear resolution of hex cells. The tree is traversed
top-down each frame:

```
function selectLOD(node, camera):
    screenSize = node.worldSize / distance(node.center, camera)
    if screenSize < LOD_THRESHOLD or node.isLeaf():
        render(node)
    else:
        for child in node.children:
            selectLOD(child, camera)
```

**Target budget:** ~500,000 visible hex cells at any time, regardless of
viewing distance. The LOD system distributes this budget: many fine cells
nearby, fewer coarse cells far away.

### Geomorphing (Smooth LOD Transitions)

To prevent visible "popping" when LOD levels switch:

1. Each hex vertex stores two positions: current LOD position and parent LOD
   position.
2. A **morph factor** (0.0 to 1.0) based on camera distance smoothly blends
   between them in the vertex shader.
3. When morph reaches 1.0, the coarser LOD substitutes in without any visual
   discontinuity.

```glsl
// Vertex shader geomorphing
vec3 morphed = mix(position_current_lod, position_parent_lod, morph_factor);
float morph_factor = smoothstep(lod_near, lod_far, camera_distance);
```

### Hex Growing/Shrinking for Dynamic LOD

As the camera approaches a chunk, its hex cells **subdivide** (grow in count,
shrink in size). As the camera retreats, cells **merge** (shrink in count,
grow in size).

With aperture 4:
- **Subdivide:** 1 parent hex → 4 child hexes (each ~half the linear size)
- **Merge:** 4 child hexes → 1 parent hex (double the linear size)

The transition uses geomorphing: child hexes start at the parent's size/position
and smoothly interpolate to their true size/position over a distance range.

### Seam Handling Between LOD Levels

Where two adjacent chunks have different LOD levels, a visible crack can appear
(T-junction). Solutions (use both):

1. **Skirts:** Each chunk extends a thin vertical "skirt" mesh downward along
   its edges. Cracks between LODs are hidden beneath the skirt. Simple and
   cheap.

2. **Transition strips:** At the boundary, insert a special strip that matches
   the fine-LOD edge on one side and the coarse-LOD edge on the other. Similar
   to the Transvoxel algorithm for cubic voxels, adapted for hex grids.

### Floating Point Precision

At planetary scale, 32-bit floats lose precision. At 1 AU (~150 billion meters),
float32 resolution is ~10 km -- useless for surface rendering.

**Solution: Camera-Relative Rendering**

```
CPU (double precision):
    planet_pos_f64 = evaluate_orbit(time)
    vertex_world_f64 = chunk_origin_f64 + local_vertex_f32
    vertex_camera_relative_f64 = vertex_world_f64 - camera_pos_f64
    vertex_for_gpu = (float)vertex_camera_relative_f64

GPU (float32):
    // All positions are small offsets from camera -- high precision
    gl_Position = projection * view * vec4(vertex_camera_relative, 1.0);
```

Combined with a **logarithmic depth buffer** to handle the z-fighting problem
across the enormous depth range (1m to millions of km):

```glsl
// Logarithmic depth buffer (fragment shader)
float Fcoef = 2.0 / log2(far_plane + 1.0);
gl_FragDepth = log2(frag_depth) * Fcoef * 0.5;
```

---

## 5. Chunking & Streaming

### Chunk Definition

A **chunk** is the fundamental unit of loading, meshing, and rendering:

```
Chunk = {
    CellID     origin;          // bottom-corner cell of this chunk
    u8         lod_level;       // resolution level
    HexVoxel   voxels[W×W×H];  // hex prism voxel data (e.g., 32×32×64)
    MeshData   mesh;            // generated triangle data
    LightData  light;           // per-voxel light values
    u8         state;           // UNLOADED | GENERATING | MESHED | ACTIVE
}
```

**Chunk size tradeoffs:**
- Smaller chunks (16³): More granular streaming, less wasted work, but more
  draw calls and more boundary stitching.
- Larger chunks (64³): Fewer draw calls, less boundary overhead, but more
  wasted work when regenerating and larger memory footprint per chunk.
- **Recommended: 32×32×64** (32 hex cells wide, 64 layers tall). Balances
  meshing speed, draw call count, and memory.

### Streaming Architecture

```
                    ┌─────────────────────┐
                    │    LOD Tree          │
                    │    Traversal         │
                    │  (main thread)       │
                    └─────────┬───────────┘
                              │ chunk requests
                              ▼
                    ┌─────────────────────┐
                    │  Priority Queue      │
                    │  (distance-sorted)   │
                    └─────────┬───────────┘
                              │
              ┌───────────────┼───────────────┐
              ▼               ▼               ▼
    ┌──────────────┐ ┌──────────────┐ ┌──────────────┐
    │  Worker 1    │ │  Worker 2    │ │  Worker N    │
    │  Generate    │ │  Generate    │ │  Generate    │
    │  Noise →     │ │  Noise →     │ │  Noise →     │
    │  Voxelize →  │ │  Voxelize →  │ │  Voxelize →  │
    │  Mesh →      │ │  Mesh →      │ │  Mesh →      │
    │  Light       │ │  Light       │ │  Light       │
    └──────┬───────┘ └──────┬───────┘ └──────┬───────┘
           │                │                │
           ▼                ▼                ▼
    ┌─────────────────────────────────────────────┐
    │  Upload Queue (main thread picks up &       │
    │  submits to GPU each frame, N chunks/frame) │
    └─────────────────────────────────────────────┘
```

### Loading/Unloading Policy

- **Load radius:** Chunks within a spherical shell around the camera (distance
  depends on LOD level -- finer LOD = smaller radius).
- **Unload radius:** Load radius × 1.5 (hysteresis prevents thrashing when
  the camera is near a boundary).
- **Priority:** Chunks directly in front of the camera (within the view
  frustum) are highest priority. Chunks behind the camera are lowest.
- **Budget:** Maximum N chunk uploads per frame (e.g., 4-8) to avoid stalling
  the GPU. Spread heavy work across frames.

### Memory Management

**Vertex Pool:** A single large GPU buffer (e.g., 1.5 GB) pre-allocated at
startup. Chunk meshes are sub-allocated from this pool using a free-list
allocator. No per-chunk buffer creation/destruction (which causes stutters).

```
┌──────────────────────────────────────────────────┐
│                  GPU Vertex Pool                  │
│  [Chunk A mesh][Chunk B mesh][  free  ][Chunk C] │
│  offset: 0     offset: 4200  offset: 8000       │
└──────────────────────────────────────────────────┘
```

**Multi-Draw Indirect:** All visible chunks are drawn in a single draw call
using an array of indirect draw commands. Each command references its
sub-allocation in the vertex pool via a base-vertex offset.

**Procedural Generation as Compression:** The planet is defined by a seed +
noise parameters. Voxel data is generated on-the-fly, never fully stored.
Only player modifications (block placed/removed) are persisted as a sparse
diff overlay on top of the procedural base.

---

## 6. Voxel Lighting

### Dual-Channel Light Model (Minecraft-inspired)

Each hex voxel stores two 4-bit light values (packed into 1 byte):

- **Block light** (0-15): Emitted by light sources (torches, lava, glowing ore).
- **Sky light** (0-15): Propagated from the sky downward.

The final light level at a voxel = `max(block_light, sky_light * sun_intensity)`.

### Flood-Fill BFS Algorithm (Adapted for Hex Grid)

Light propagation uses breadth-first search across the hex voxel neighbor graph.

**Hex prism neighbors (8 total):**
- 6 lateral (in-plane hex neighbors via cube coordinates)
- 2 vertical (above and below in the column stack)

```
function propagate_light(source_cell, light_value):
    queue = [(source_cell, light_value)]
    while queue is not empty:
        (cell, value) = queue.pop_front()
        for neighbor in cell.get_8_neighbors():  // 6 lateral + 2 vertical
            attenuation = 1 + neighbor.opacity
            new_value = value - attenuation
            if new_value > neighbor.light_level:
                neighbor.light_level = new_value
                queue.push_back((neighbor, new_value))
```

**Light removal** (when a light source is destroyed):
1. **Removal BFS:** Starting from the removed source, flood outward, zeroing
   light values and collecting "boundary" cells (cells that still have light
   from other sources).
2. **Re-propagation BFS:** From the boundary cells, re-propagate light
   normally. This correctly restores the lighting from remaining sources.

### Sky Light on a Spherical Planet

On a sphere, "up" varies by location, and sunlight comes from a direction
that depends on the planet's rotation and orbital position.

**Implementation:**
1. For each surface chunk, compute the local sun direction:
   `sun_dir = normalize(star_position - planet_position)` rotated into the
   planet's body frame accounting for axial rotation.
2. For each hex column, if the top is exposed to the sky AND `dot(sun_dir,
   local_up) > 0` (sun is above horizon), set sky light to 15.
3. Sky light propagates downward without attenuation (straight down the column)
   and laterally with -1 attenuation per step (same BFS as block light).
4. At the terminator (day/night boundary), sky light = 0 for columns where
   the sun is below the horizon.

**Amortized updates:** The sun moves slowly. Update sky light for a few chunks
per frame in round-robin order, prioritizing chunks near the terminator.

### Ambient Occlusion

Per-vertex ambient occlusion for hex prisms, adapted from the cubic AO
technique by Mikola Lysenko:

- For each vertex of a hex face, sample the 3 adjacent hex cells that share
  that corner.
- Count how many are solid (occluding).
- AO factor = 1.0 (no occlusion) down to 0.0 (fully occluded).
- Interpolate AO across the face for smooth contact shadows.

### Performance Budget

| Operation | Budget | Strategy |
|---|---|---|
| Light propagation (block placed/removed) | < 1 ms | Incremental BFS, typically < 2000 cells |
| Chunk initial lighting | < 2 ms | Full BFS on worker thread during chunk generation |
| Sky light update (per chunk) | < 0.5 ms | Amortized, 4-8 chunks per frame |
| AO computation | Baked into mesh | Computed during mesh generation, zero runtime cost |

### Future: Advanced Lighting

If we want to go beyond Minecraft-style discrete lighting later:

- **Light Propagation Volumes (LPV):** GPU-based, one-bounce diffuse GI using
  spherical harmonics on a 3D grid. Medium quality, medium cost.
- **Voxel Cone Tracing:** Voxelize the scene into a 3D texture mipmap, trace
  cones for diffuse + specular GI. High quality, high cost.
- **Screen-Space AO (SSAO/GTAO):** Layer on top of baked AO for dynamic
  objects. Low cost, good visual improvement.

Start with flood-fill BFS + baked AO. It's proven, fast, and fits the hex
aesthetic. Add SSAO on top for visual polish.

---

## 7. Collision & Movement

### No Physics Engine

This is a voxel game like Minecraft. Minecraft doesn't use a physics engine
and neither do we. All collision is direct voxel lookup — no Jolt, no Rapier,
no Bullet. The world *is* the collision data.

**Why no physics library:**
- Voxel terrain is axis-aligned prisms — collision is a grid lookup, not a
  geometry intersection problem.
- Player movement is simple: walk, jump, fall, swim. No ragdolls, no rigid
  body stacking, no joints.
- A physics engine adds a large dependency, threading complexity, and a whole
  category of bugs (tunneling, jitter, solver instability) for zero benefit.
- Minecraft handles millions of players with pure voxel collision. It works.

### AABB vs Hex Prism Collision

Hex prisms are trickier than cubes for collision because their cross-section
isn't axis-aligned. Two approaches:

**Option A: Approximate with AABB (recommended to start)**

Treat each hex prism as its axis-aligned bounding box. The AABB of a regular
hexagon with circumradius `r` is `2r × r√3` (width × height in the hex plane).

```
       ___
      /   \        AABB approximation:
     /     \       ┌─────────┐
    /       \      │         │
    \       /      │   hex   │
     \     /       │         │
      \___/        └─────────┘
```

- Simple, fast, identical to Minecraft's collision model
- Slight inaccuracy at hex corners (AABB extends ~15% beyond the hex edge)
- At voxel scale (1-2 m cells), the error is imperceptible
- All standard swept-AABB algorithms apply directly

**Option B: True hex prism collision (upgrade later if needed)**

For precise collision, test against the hex polygon using a point-in-hexagon
test (6 half-plane checks). Only needed if gameplay demands pixel-perfect hex
boundaries (unlikely for a Minecraft-style game).

### Movement & Collision Algorithm

Standard Minecraft-style swept collision against the voxel grid:

```
function move_player(player, velocity, dt):
    desired = velocity * dt

    // Resolve each axis independently (separating axis)
    // Order: Y first (gravity), then X, then Z
    for axis in [Y, X, Z]:
        new_pos = player.pos + desired[axis]
        player_aabb = aabb_at(new_pos, player.size)

        // Check all voxels the AABB overlaps
        for cell in get_overlapping_cells(player_aabb):
            if cell.is_solid():
                // Push player out along this axis
                desired[axis] = resolve_overlap(player_aabb, cell_aabb, axis)

        player.pos[axis] += desired[axis]

    // Ground check: is the voxel directly below the player solid?
    player.on_ground = is_solid(player.pos - (0, 0.01, 0))
```

**Key details:**
- **Y-first resolution** ensures gravity works before lateral movement
  (prevents sliding off ledges, sticking to walls)
- **Separating axis** resolution handles corners correctly
- **Sub-cell precision:** Player position is floating-point, collision grid
  is integer cell coordinates. `get_overlapping_cells` converts the AABB to
  a cell range and iterates.

### Collision on a Spherical Planet

On a sphere, "down" varies by location. The movement algorithm needs
adaptation:

**Local tangent frame:** At the player's current position, define a local
coordinate frame where:
- **Up** = surface normal (away from planet center)
- **Forward/Right** = tangent to the sphere surface

All movement and collision math runs in this local frame. The frame rotates
as the player moves across the surface.

```
local_up = normalize(player_pos - planet_center)
local_forward = normalize(cross(world_up_hint, local_up))
local_right = cross(local_up, local_forward)

// Gravity always points toward planet center
gravity = -local_up * g  // g ≈ 9.81 m/s²
```

**Voxel lookup on a sphere:** Given a world-space position, convert to the
hex grid coordinate system (face + axial q,r + layer y) to look up whether
that cell is solid. This is the same coordinate conversion used for rendering
and chunk streaming — no additional data structure needed.

### Performance

Collision is essentially free compared to mesh generation and rendering:

| Operation | Cost | Notes |
|---|---|---|
| Player collision (per frame) | < 0.01 ms | ~20-50 cell lookups |
| Entity collision (per entity) | < 0.01 ms | Same algorithm |
| 100 entities | < 1 ms total | Trivially parallelizable |

---

## 8. Orbital Mechanics & Day/Night

### Keplerian Orbits (On-Rails)

All celestial bodies (planets, moons, asteroids) follow pre-computed
**Keplerian orbits**. Given the 6 orbital elements, the position at any time
T is computed analytically -- no numerical integration, O(1) per body.

**Orbital elements:**
```
a   = semi-major axis (orbit size)
e   = eccentricity (0 = circle, 0-1 = ellipse)
i   = inclination (tilt relative to reference plane)
Ω   = longitude of ascending node (orientation in reference plane)
ω   = argument of periapsis (orientation of ellipse within orbital plane)
M₀  = mean anomaly at epoch (starting position)
```

**Position at time T:**
1. Compute mean anomaly: `M = M₀ + n * (T - T₀)` where `n = sqrt(μ/a³)`
2. Solve Kepler's equation iteratively: `M = E - e * sin(E)` (Newton-Raphson,
   converges in 3-5 iterations)
3. Compute true anomaly ν from eccentric anomaly E
4. Compute radius: `r = a * (1 - e * cos(E))`
5. Convert (r, ν) to 3D position using the orbital element rotation matrices

**Advantage:** Orbits are perfectly stable forever. No energy drift. Time warp
at any speed is free (just evaluate at T + dt). This is how KSP handles
"on-rails" bodies.

### Patched Conics for Spacecraft

When the player's ship is actively thrusting or transitioning between
gravitational spheres of influence:

1. **Sphere of Influence (SOI):** Each body has an SOI radius:
   `r_SOI = R * (m_body / m_parent)^(2/5)`
2. While inside body A's SOI, solve a two-body Keplerian orbit around A.
3. When the trajectory crosses the SOI boundary, transform the state vector
   to the parent body's frame and compute a new orbit.
4. An interplanetary transfer = 3 patched conics: departure hyperbola →
   transfer ellipse → arrival hyperbola.

### Reference Frame Hierarchy

```
Solar Frame (Heliocentric Inertial)
 ├── Planet A Frame (Planet-Centered Inertial)
 │    ├── Planet A Surface Frame (Rotating with planet)
 │    │    └── Local Tangent Frame (physics, rendering)
 │    └── Moon A1 Frame
 │         └── Moon A1 Surface Frame
 └── Planet B Frame
      └── ...
```

**Switching logic:**
- On a planet's surface → Surface Frame (rotates with planet)
- In low orbit → Planet-Centered Inertial Frame
- On interplanetary trajectory → Solar Frame
- Transitions are seamless: convert state vector between frames

### Day/Night Cycle

**Computing sun direction at any surface point:**

```
1. sun_world = normalize(star_pos - planet_pos)     // in solar frame
2. rotation_angle = 2π × time / rotation_period     // planet spin
3. axial_tilt = rotation_matrix(obliquity_axis, obliquity_angle)
4. spin = rotation_matrix(pole_axis, rotation_angle)
5. sun_local = inverse(axial_tilt × spin) × sun_world  // planet body frame
6. For surface point at (lat, lon):
     local_up = surface_normal(lat, lon)
     sun_elevation = asin(dot(sun_local, local_up))
     → Daytime if sun_elevation > 0
     → Night if sun_elevation < 0
     → Dawn/dusk in the transition zone
```

**Visual implementation:**
- Sun elevation drives sky light intensity for voxel lighting
- Sun elevation drives voxel face shading (see below)
- Sun elevation drives simple sky color gradient

### Voxel Face Shading

Shading is simple — we're rendering voxels, not PBR materials. Each hex prism
face gets a brightness based on two factors:

1. **Face normal vs sun direction** (directional light)
2. **Face normal vs local up** (ambient directional bias)

```glsl
// Fragment shader: simple voxel face shading
float sun_factor = max(dot(face_normal, sun_direction), 0.0);
float ambient = 0.3;  // base ambient so shadowed faces aren't black

// Voxel light from BFS (Section 6) provides the other component
float voxel_light = texelFetch(light_data, cell_id).r / 15.0;

float brightness = max(ambient + sun_factor * 0.7, voxel_light);
vec3 color = block_color * brightness;
```

**Per-face shading constants** (Minecraft-style, baked into vertex color):

| Face Direction | Brightness Multiplier | Rationale |
|---|---|---|
| Top | 1.0 | Fully lit by sky |
| Bottom | 0.5 | Darkest (rarely seen) |
| Side (sun-facing) | 0.8 - 1.0 | Depends on `dot(normal, sun_dir)` |
| Side (shadow-facing) | 0.6 | Ambient only |

This is the same approach Minecraft uses. Combined with the BFS voxel lighting
from Section 6, it produces the characteristic voxel look with minimal shader
complexity.

### Sky Rendering

A simple gradient sky based on sun elevation — no atmospheric scattering
simulation needed:

```glsl
// Sky color: blend between day/sunset/night based on sun elevation
vec3 day_color   = vec3(0.4, 0.7, 1.0);   // blue sky
vec3 dusk_color  = vec3(1.0, 0.5, 0.2);   // orange horizon
vec3 night_color = vec3(0.02, 0.02, 0.05); // dark blue

float sun_elev = dot(sun_direction, local_up);  // -1 to +1

vec3 sky;
if (sun_elev > 0.1)
    sky = day_color;
else if (sun_elev > -0.1)
    sky = mix(dusk_color, day_color, (sun_elev + 0.1) / 0.2);
else
    sky = mix(night_color, dusk_color, max(sun_elev + 0.3, 0.0) / 0.2);
```

**Planet surface angle:** The key insight for shading on a spherical planet is
that `local_up` (surface normal) varies across the planet. A hex face on the
equator facing the sun gets full light. The same face near the terminator
(day/night boundary) gets grazing light. This happens naturally from
`dot(face_normal, sun_direction)` because face normals are in world space,
already accounting for the planet's curvature.

### Future: Visual Upgrades

If we want to improve visuals later, these are additive — they don't change
the core shading model:

- **Simple shadow mapping:** Single shadow map from the sun for nearby terrain.
  No cascades needed initially.
- **SSAO:** Screen-space ambient occlusion layered on top of baked per-vertex AO.
  Low cost, good visual improvement.
- **Simple atmospheric scattering:** O'Neil single-scattering model (GPU Gems 2)
  for sunrise/sunset color. Only worth adding if the gradient sky looks too flat.

---

## 10. Multiplayer & Networking

### Architecture: Host-as-Server with WebRTC

One player acts as the **authoritative host** (dedicated or player-hosted). All
other players connect to the host via **WebRTC data channels**. This gives us
P2P-like latency without requiring a dedicated game server for every session.

```
┌──────────────────────────────────────────────────────────────┐
│                    NETWORK TOPOLOGY                           │
│                                                              │
│              ┌─────────────────────┐                         │
│              │   Vercel Signaling  │ ◄── Serverless funcs    │
│              │   (matchmaking +    │     (zero idle cost)    │
│              │    SDP relay)       │                         │
│              └──────────┬──────────┘                         │
│                 ▲   ▲   │   ▲                                │
│     SDP offer/  │   │   │   │  SDP offer/answer             │
│     answer      │   │   │   │                                │
│              ┌──┴───┴───┴───┴──┐                             │
│              │   HOST (Player)  │ ◄── Authoritative state    │
│              │   - World state  │     - Validates edits      │
│              │   - Chunk gen    │     - Runs collision         │
│              │   - Persistence  │     - Broadcasts updates   │
│              └──┬────┬────┬────┘                             │
│     WebRTC      │    │    │   WebRTC                         │
│     DataChannel │    │    │   DataChannel                    │
│              ┌──┴──┐ │ ┌──┴──┐                               │
│              │ P2  │ │ │ P3  │  ... up to ~8-16 peers        │
│              └─────┘ │ └─────┘                               │
│                   ┌──┴──┐                                    │
│                   │ P4  │                                    │
│                   └─────┘                                    │
└──────────────────────────────────────────────────────────────┘
```

**Why host-as-server over full mesh:**
- Authoritative state prevents desync (one source of truth for voxel edits)
- Bandwidth scales linearly (N connections) vs quadratically (N² for full mesh)
- Simpler conflict resolution for simultaneous block edits
- Host can run headless for dedicated server scenarios

**Why WebRTC over WebSocket:**
- UDP-like unreliable channels for position updates (lower latency)
- Reliable ordered channels for voxel edits (TCP-like guarantees)
- P2P after signaling — traffic doesn't route through a relay server
- Works in browsers (web build support)

### Signaling & Matchmaking (Vercel)

The signaling server runs as **Vercel serverless functions** + a lightweight
real-time layer. It handles:

1. **Lobby management:** Create/list/join game sessions
2. **SDP relay:** Exchange WebRTC Session Description Protocol offers/answers
3. **ICE candidate relay:** Forward STUN/TURN ICE candidates between peers
4. **NAT traversal:** Provide STUN server list; fall back to TURN relay for
   symmetric NAT (worst case)

```
Signaling Flow:
─────────────────────────────────────────────────────────
Host                  Vercel                  Client
  │                     │                       │
  │── POST /lobby ─────►│                       │
  │◄── lobby_id ────────│                       │
  │                     │                       │
  │                     │◄── GET /lobbies ──────│
  │                     │── lobby list ────────►│
  │                     │                       │
  │                     │◄── POST /join ────────│
  │◄── client SDP offer─│                       │
  │── SDP answer ──────►│                       │
  │                     │── SDP answer ────────►│
  │                     │                       │
  │◄════ WebRTC DataChannel established ═══════►│
  │         (Vercel no longer in the path)       │
  │                     │                       │
```

**Vercel stack:**
- **API Routes** (serverless functions): `/api/lobby/create`, `/api/lobby/list`,
  `/api/lobby/join`, `/api/signal`
- **Vercel KV** (Redis): Transient lobby state, active session registry
- **STUN:** Google's public STUN servers (`stun.l.google.com:19302`)
- **TURN fallback:** Xirsys, Twilio, or Cloudflare Calls (for symmetric NAT)

**Cost model:** Signaling is lightweight (~1-5 KB per connection setup). Vercel's
free tier handles hundreds of concurrent matchmaking requests. Game traffic is
pure P2P after signaling — zero ongoing server cost.

### WebRTC Data Channels

Each client-to-host connection uses **two data channels** with different
reliability modes:

| Channel | Mode | Use Case | Rate |
|---|---|---|---|
| `state` | Unreliable, unordered | Position, rotation, animation state | 20-30 Hz |
| `world` | Reliable, ordered | Voxel edits, chat, events, chunk requests | On demand |

**Why two channels:** Position updates must be low-latency and can tolerate
loss (old positions are superseded). Voxel edits must be reliable and ordered
(a place then break sequence must not arrive reversed).

### Synchronization Model

#### World State (Host-Authoritative)

```
┌─────────────────────────────────────────────────────────┐
│                HOST GAME LOOP                            │
│                                                         │
│  1. Receive client inputs (movement, edits)             │
│  2. Validate edits (anti-cheat, permission checks)      │
│  3. Apply valid edits to world state                    │
│  4. Step collision & movement                           │
│  5. Broadcast state updates to all clients              │
│  6. Persist edits to local database (async)             │
└─────────────────────────────────────────────────────────┘
```

#### Chunk Streaming to Clients

Clients **do not** independently generate terrain from the seed alone, because
the host may have persisted player edits that modify the procedural base. Instead:

1. Client enters a region → sends `ChunkRequest(chunk_id, lod_level)` to host
2. Host checks if the chunk has modifications:
   - **No modifications:** Send only the seed + noise parameters. Client
     generates locally (saves bandwidth).
   - **Has modifications:** Send seed + noise parameters + sparse edit diff.
     Client generates base terrain, then applies the diff.
3. Client generates mesh locally from the received data.

```
Message: ChunkResponse
{
    chunk_id:    u64,
    lod_level:   u8,
    has_edits:   bool,
    // If no edits, client generates from seed (zero voxel data sent)
    // If edits exist, send only the diff:
    edit_count:  u16,
    edits:       [{cell_id: u64, old_voxel: u16, new_voxel: u16}],
}
```

**Bandwidth estimate:** An unmodified chunk = ~20 bytes (chunk ID + LOD +
generation params). A modified chunk = 20 + (8 bytes × edit count). For a
typical build of 500 edits, that's ~4 KB — trivial.

#### Player Position Sync

```
Message: PlayerState (unreliable channel, 20 Hz)
{
    player_id:   u16,
    position:    [f64; 3],    // camera-relative on receiver side
    velocity:    [f32; 3],
    orientation: [f32; 4],    // quaternion
    animation:   u8,          // current animation state
    timestamp:   u32,         // for interpolation ordering
}
```

**Client-side prediction:** Each client simulates its own movement locally for
instant responsiveness. The host validates and sends corrections. On correction,
the client replays inputs from the corrected state (reconciliation).

**Entity interpolation:** Remote players are rendered with a ~100 ms interpolation
buffer. Two recent states are interpolated between for smooth motion, even with
jitter or packet loss.

#### Voxel Edit Sync

```
Edit Flow:
────────────────────────────────────────
Client                  Host              Other Clients
  │                       │                    │
  │── EditRequest ──────►│                    │
  │   {cell, new_voxel}  │                    │
  │                      │── validate ──►     │
  │                      │                    │
  │◄── EditConfirm ──────│── EditBroadcast ──►│
  │   {cell, new_voxel,  │   {cell, new_voxel,│
  │    sequence_num}     │    player_id}      │
  │                      │                    │
```

- **Optimistic application:** Client applies the edit locally immediately for
  responsiveness, but marks it as unconfirmed.
- **Rollback on rejection:** If the host rejects the edit (e.g., protected area,
  invalid position, rate limit), the client reverts.
- **Sequence numbers:** Ensure edits are applied in the correct order. The host
  assigns a monotonically increasing sequence number to each confirmed edit.

### Bandwidth & Scalability

| Players | Position Traffic (20 Hz) | Typical Edit Traffic | Total Estimate |
|---|---|---|---|
| 2 | ~2 KB/s | ~0.5 KB/s | ~2.5 KB/s |
| 4 | ~6 KB/s | ~1.5 KB/s | ~7.5 KB/s |
| 8 | ~14 KB/s | ~3 KB/s | ~17 KB/s |
| 16 | ~30 KB/s | ~6 KB/s | ~36 KB/s |

These are well within typical home upload bandwidth (~5-50 Mbps). The bottleneck
is chunk serving when many players explore different regions simultaneously.

**Mitigation for chunk serving:** Clients that have already generated a chunk
can serve it to newly connecting peers (BitTorrent-style chunk sharing), reducing
host bandwidth for chunk data.

### Session Management

- **Host migration:** If the host disconnects, the player with the lowest latency
  to the most peers is promoted to new host. The departing host's persisted world
  state must be transferred (or the new host regenerates from seed + available
  edit history).
- **Save on disconnect:** The host periodically checkpoints the world state.
  On graceful disconnect, a final save is triggered.
- **Reconnection:** Clients reconnect via the signaling server and receive a
  full state sync (current chunk edits + player positions) on rejoin.

### Anti-Cheat (Lightweight)

Since the host is authoritative:
- **Edit rate limiting:** Max N edits per second per player
- **Range checks:** Player can only edit cells within interaction range
- **Material validation:** Player must possess the block type they're placing
- **Movement validation:** Host checks for impossible speeds/teleportation

This is not a competitive game, so heavy anti-cheat (kernel drivers, etc.) is
unnecessary. Trust the host, validate the basics.

---

## 11. Persistence & Data Storage

### Design Principles

1. **Procedural generation is the base layer.** The planet is fully defined by
   its seed + noise parameters. We never store unmodified terrain.
2. **Only store the diff.** Player modifications (blocks placed, removed, or
   changed) are stored as a sparse overlay on top of the procedural base.
3. **Host-local storage.** The host player's machine stores the world database.
   No cloud storage required (but could be added later).

### Database: SQLite

**SQLite** is the storage engine:

| Requirement | SQLite Fit |
|---|---|
| Embedded (no server process) | Single file, zero config |
| Concurrent reads + single writer | WAL mode supports this |
| Structured queries (find edits near X) | Full SQL with spatial indexing |
| Cross-platform | C library, works everywhere including WASM |
| Battle-tested | Used in every smartphone, browser, aircraft |
| Backup/transfer | Copy one file |

**Why not LevelDB/RocksDB:** They're key-value stores — good for simple lookups
but lack the relational queries we need (e.g., "all edits within chunk X",
"all builds by player Y"). SQLite gives us both fast point lookups and flexible
queries with a single dependency.

**Why not a custom binary format:** We'd end up reimplementing half of SQLite
poorly. The overhead of SQLite over raw file I/O is negligible for our data
volumes, and we get transactions, crash recovery, and indexing for free.

### Schema

```sql
-- World metadata
CREATE TABLE world (
    id          INTEGER PRIMARY KEY,
    seed        INTEGER NOT NULL,          -- procedural generation seed
    name        TEXT NOT NULL,
    created_at  INTEGER NOT NULL,          -- unix timestamp
    play_time   INTEGER DEFAULT 0,         -- total seconds played
    version     INTEGER DEFAULT 1          -- schema migration version
);

-- Noise parameters per planet (procedural generation config)
CREATE TABLE planet (
    id              INTEGER PRIMARY KEY,
    world_id        INTEGER NOT NULL REFERENCES world(id),
    name            TEXT,
    seed            INTEGER NOT NULL,
    radius          REAL NOT NULL,              -- meters
    orbital_elements TEXT NOT NULL,             -- JSON: {a, e, i, omega, Omega, M0}
    noise_config    BLOB NOT NULL,             -- serialized noise layer stack
    biome_config    BLOB NOT NULL              -- serialized biome parameters
);

-- Sparse voxel edit overlay (the core persistence table)
CREATE TABLE voxel_edit (
    id          INTEGER PRIMARY KEY,
    planet_id   INTEGER NOT NULL REFERENCES planet(id),
    cell_id     INTEGER NOT NULL,              -- packed 64-bit cell ID
    chunk_id    INTEGER NOT NULL,              -- for efficient chunk-level queries
    old_voxel   INTEGER NOT NULL,              -- original procedural voxel type
    new_voxel   INTEGER NOT NULL,              -- player-modified voxel type
    player_id   INTEGER,                       -- who made the edit (NULL = system)
    timestamp   INTEGER NOT NULL,              -- unix timestamp
    UNIQUE(planet_id, cell_id)                 -- one edit per cell (latest wins)
);
CREATE INDEX idx_voxel_edit_chunk ON voxel_edit(planet_id, chunk_id);

-- Player data
CREATE TABLE player (
    id              INTEGER PRIMARY KEY,
    name            TEXT NOT NULL,
    position        BLOB,                      -- serialized: planet_id + lat/lon/alt
    inventory       BLOB,                      -- serialized inventory state
    last_seen       INTEGER                    -- unix timestamp
);

-- Named locations / waypoints
CREATE TABLE waypoint (
    id          INTEGER PRIMARY KEY,
    planet_id   INTEGER NOT NULL REFERENCES planet(id),
    player_id   INTEGER NOT NULL REFERENCES player(id),
    name        TEXT NOT NULL,
    cell_id     INTEGER NOT NULL,
    created_at  INTEGER NOT NULL
);

-- Structures / builds (grouped edits for efficient loading)
CREATE TABLE structure (
    id          INTEGER PRIMARY KEY,
    planet_id   INTEGER NOT NULL REFERENCES planet(id),
    player_id   INTEGER REFERENCES player(id),
    name        TEXT,
    origin_cell INTEGER NOT NULL,              -- anchor cell ID
    bbox_min    INTEGER NOT NULL,              -- bounding box for spatial queries
    bbox_max    INTEGER NOT NULL,
    created_at  INTEGER NOT NULL
);

-- Link edits to structures for group operations (copy, undo, protect)
CREATE TABLE structure_edit (
    structure_id INTEGER NOT NULL REFERENCES structure(id),
    edit_id      INTEGER NOT NULL REFERENCES voxel_edit(id),
    PRIMARY KEY (structure_id, edit_id)
);
```

### Read/Write Patterns

**Loading a chunk (host):**
```
1. Client requests chunk (planet_id, chunk_id, lod_level)
2. SELECT cell_id, new_voxel FROM voxel_edit
   WHERE planet_id = ? AND chunk_id = ?
3. If no rows → chunk is unmodified, client generates from seed
4. If rows exist → send seed params + edit list to client
```

**Saving an edit (host):**
```
1. Host validates the edit
2. INSERT OR REPLACE INTO voxel_edit (planet_id, cell_id, chunk_id,
   old_voxel, new_voxel, player_id, timestamp)
   VALUES (?, ?, ?, ?, ?, ?, ?)
3. Broadcast edit to connected clients
```

**Write batching:** Edits are buffered in memory and flushed to SQLite in
batches (every 5 seconds or every 100 edits, whichever comes first) inside a
single transaction. This avoids the overhead of one transaction per edit while
keeping data loss window small.

```c
// Pseudocode: batched write
if (edit_buffer.count >= 100 || time_since_flush >= 5.0) {
    sqlite3_exec(db, "BEGIN");
    for (edit in edit_buffer) {
        sqlite3_bind_and_step(insert_stmt, edit);
    }
    sqlite3_exec(db, "COMMIT");
    edit_buffer.clear();
}
```

### World Save Files

A world save is a single `.hexworld` directory:

```
saves/
└── my_planet/
    ├── world.db          ← SQLite database (all tables above)
    ├── world.db-wal      ← WAL file (auto-managed by SQLite)
    └── thumbnails/
        └── planet_1.png  ← Preview image for save selection screen
```

**Backup:** Copy the directory. SQLite in WAL mode requires copying both
`world.db` and `world.db-wal` for a consistent backup (or use SQLite's
online backup API).

**World transfer (multiplayer):** When a new host takes over, the departing
host can send the `world.db` file to the new host. At ~8 bytes per edit,
even 1 million edits = ~8 MB — easily transferable.

### Data Lifecycle

```
┌──────────────────────────────────────────────────────────┐
│                  DATA LIFECYCLE                            │
│                                                          │
│  Procedural Base (seed)                                  │
│      │                                                   │
│      ▼                                                   │
│  Generate chunk voxels on-the-fly                        │
│      │                                                   │
│      ▼                                                   │
│  Player makes edit ──► Buffer in memory                  │
│      │                      │                            │
│      ▼                      ▼                            │
│  Broadcast to clients   Flush to SQLite (batched)        │
│      │                      │                            │
│      ▼                      ▼                            │
│  Clients apply edit     On next chunk load, edits are    │
│  to local mesh          merged with procedural base      │
│                                                          │
│  On world close ──► Final flush + WAL checkpoint         │
│  On world open  ──► Load world.db, resume                │
└──────────────────────────────────────────────────────────┘
```

### Undo System

The `old_voxel` field in `voxel_edit` enables undo:

- **Single undo:** Restore `old_voxel` at `cell_id`, delete the edit row.
- **Structure undo:** Delete all edits linked via `structure_edit`, restoring
  original terrain.
- **History limit:** Keep last N edits per player for undo. Older edits become
  permanent (delete `old_voxel` to save space, or archive).

### Storage Estimates

| Scenario | Edit Count | DB Size | Notes |
|---|---|---|---|
| Light play (1 hour) | ~1,000 | ~50 KB | Small builds |
| Medium play (10 hours) | ~50,000 | ~2 MB | Several structures |
| Heavy building (100 hours) | ~1,000,000 | ~40 MB | Large base + tunnels |
| Extreme (1000+ hours) | ~10,000,000 | ~400 MB | Massive megastructure |

These are tiny by modern standards. SQLite handles databases up to 281 TB.
Performance remains excellent up to millions of rows with proper indexing.

### Future: Cloud Sync

If cloud persistence is desired later:
- **Turso** (libSQL): SQLite-compatible, edge-distributed, Vercel-friendly
- Sync the local `world.db` to Turso for cross-device play
- Or use Vercel Blob Storage for simple file backup of the `.hexworld` directory

---

## 9. Reference Material

### Essential Papers & Articles

| Topic | Reference |
|---|---|
| Goldberg polyhedra construction | Bapat et al. 2022, "Extending Goldberg's method" (Royal Society) |
| Hex grid math (definitive) | Red Blob Games: "Hexagonal Grids" |
| ISEA projection | Snyder 1992, "An Equal-Area Map Projection for Polyhedral Globes" |
| CDLOD terrain | Strugar 2009, "Continuous Distance-Dependent LOD for Heightmaps" |
| Transvoxel (LOD seams) | Lengyel 2010, "Voxel-Based Terrain for Real-Time Virtual Simulations" |
| Voxel lighting | 0fps blog: "Voxel Lighting" by Mikola Lysenko |
| Flood-fill light propagation | Seeds of Andromeda: "Fast Flood Fill Lighting" |
| Ambient occlusion for voxels | 0fps: "AO for Minecraft-like Worlds" |
| Greedy meshing | 0fps: "Meshing in a Minecraft Game" |
| Voxel collision (Minecraft-style) | Minecraft Wiki: "Entity physics", 0fps: "Collision Detection for Voxels" |
| Vertex pooling for voxels | Nick McDonald: "High Performance Voxel Engine" |
| Logarithmic depth buffer | Outerra blog, Brano Kemen |
| Patched conics | Wikipedia: "Patched Conic Approximation" |
| Camera-relative rendering | Godot: "Emulating Double Precision on GPU" |
| DGGS theory | Sahr, White, Kimerling 2003: "Geodesic Discrete Global Grid Systems" |
| WebRTC networking | MDN: WebRTC API, "High Performance Browser Networking" by Ilya Grigorik |
| Client-side prediction | Gabriel Gambetta: "Fast-Paced Multiplayer" (4-part series) |
| Netcode for voxel games | Gaffer on Games: "Networked Physics" by Glenn Fiedler |
| SQLite for games | sqlite.org: "Appropriate Uses for SQLite" |

### Key Codebases to Study

| Project | What to Learn |
|---|---|
| Uber H3 (github.com/uber/h3) | Hex sphere indexing, 64-bit cell IDs, neighbor traversal |
| DGGRID (discreteglobalgrids.org) | ISEA projection implementation, aperture variants |
| sp4cerat/Planet-LOD | Adaptive spherical LOD with triangle subdivision |
| fstrugar/CDLOD | Reference CDLOD implementation |
| SebLague/Procedural-Planets | Clean Unity planet generation (noise layering) |
| cgerikj/binary-greedy-meshing | Ultra-fast bitmask voxel meshing |
| google/s2geometry | Hierarchical sphere indexing (quad-based, but useful patterns) |
| floooh/sokol-samples | Sokol rendering examples and patterns |

### Comparable Engines

| Engine | Lesson |
|---|---|
| Outerra | Seamless space-to-ground, logarithmic depth, cube-sphere LOD |
| KSP (PQS system) | Layered procedural terrain mods, scaled-space dual representation |
| No Man's Sky | 3D voxel density field, GPU marching cubes, seed-based determinism |
| Hytale (The Forge) | Voxel vertex packing (77% memory reduction), visibility buffers |
| Minecraft | Flood-fill lighting, greedy meshing, chunk streaming, AABB voxel collision |

---

## Next Steps

### Phase 1: Foundation (Weeks 1-4)
- [ ] Set up project with Sokol (or chosen stack) + build system
- [ ] Implement basic Goldberg polyhedron generation (GP(4,0) → GP(64,0))
- [ ] Render a wireframe hex sphere with camera controls
- [ ] Implement axial coordinate system and 64-bit cell ID

### Phase 2: Terrain (Weeks 5-8)
- [ ] Implement ISEA projection with aperture 4 subdivision
- [ ] Add procedural noise-based terrain height (FastNoiseLite)
- [ ] Generate hex prism voxel meshes for terrain chunks
- [ ] Implement basic chunked LOD (2-3 levels)

### Phase 3: Core Systems (Weeks 9-12)
- [ ] Multithreaded chunk generation pipeline
- [ ] Flood-fill voxel lighting (block light + sky light)
- [ ] Minecraft-style AABB voxel collision (swept, per-axis resolution)
- [ ] First-person camera walking on planet surface

### Phase 4: Scale (Weeks 13-16)
- [ ] Full LOD tree (8+ levels, space to surface)
- [ ] Camera-relative rendering + logarithmic depth buffer
- [ ] Chunk streaming with priority queue
- [ ] Vertex pool GPU memory management

### Phase 5: Solar System (Weeks 17-20)
- [ ] Keplerian orbit evaluation
- [ ] Multiple planets with LOD
- [ ] Day/night cycle with sky light propagation
- [ ] Sun-angle voxel face shading + gradient sky
- [ ] Reference frame switching (surface ↔ orbit ↔ solar)

### Phase 6: Polish (Weeks 21-24)
- [ ] Geomorphing LOD transitions
- [ ] Per-vertex ambient occlusion (baked into mesh)
- [ ] Web build (Emscripten + WebGPU)
- [ ] Player terrain modification (place/break hex voxels)

### Phase 7: Persistence (Weeks 25-28)
- [ ] Integrate SQLite (single-file world.db)
- [ ] Implement voxel edit table + sparse diff overlay
- [ ] Batched write pipeline (buffer → flush every 5s / 100 edits)
- [ ] World save/load (.hexworld directory)
- [ ] Undo system (per-player edit history)
- [ ] Planet table with seed + noise config storage

### Phase 8: Multiplayer (Weeks 29-34)
- [ ] WebRTC data channel abstraction (reliable + unreliable channels)
- [ ] Vercel signaling server (lobby create/list/join, SDP relay)
- [ ] Host-authoritative game loop (validate edits, broadcast state)
- [ ] Chunk streaming to clients (seed-only vs seed+diff)
- [ ] Player position sync (client prediction + interpolation)
- [ ] Voxel edit sync (optimistic apply + rollback on rejection)
- [ ] Host migration on disconnect
- [ ] NAT traversal (STUN + TURN fallback)

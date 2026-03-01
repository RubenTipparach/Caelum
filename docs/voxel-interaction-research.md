# Voxel Interaction & Physics Research

## Sources: Tenebris (TypeScript/Three.js cubic voxels) + Hex-Planets (C/Sokol hex voxels)

---

## Current State (Hex-Planets)

### What Works
- **Raycast:** Step-based ray marching (0.1m steps), tangent-plane projection to hex grid, radius comparison for hit detection
- **Break/Place:** LMB = decrement column height, RMB = increment + set to STONE. Chunk marked dirty for remesh
- **Highlight:** 6-edge hex wireframe outline via `highlight.glsl` (LINES primitive)
- **Ground collision:** `lod_tree_terrain_height()` returns terrain + arch height. Player snaps to ground with smoothed interpolation (instant snap down, 3Hz lowpass up)
- **Jetpack:** Double-tap SPACE, speed scales with altitude

### What's Missing (the full Minecraft experience)
1. **3D voxel storage** — currently per-column `int16_t heights[]` + `uint8_t types[]` (only stores top surface, not full column)
2. **Wall collision** — player teleports to top of any voxel column instead of being blocked by walls
3. **Headroom check** — can walk into 1-block-tall spaces
4. **Multi-block structures** — no way to have floating blocks, tunnels, arches as actual voxel data
5. **Per-voxel break/place** — currently only raises/lowers column height, can't carve a tunnel or place a block on the side

---

## 1. 3D Voxel Storage (Required Foundation)

### Current: Heightfield (1 value per column)
```c
// hex_terrain.h
int16_t heights[32][32];  // terrain height per hex column
uint8_t types[32][32];    // voxel type at surface only
```
This can't represent caves, tunnels, floating slabs, or any structure with air gaps.

### Target: Full 3D Voxel Array (1 value per cell per layer)
```c
// Like Tenebris:
// blocks[x + z*CHUNK_SIZE + y*CHUNK_SIZE*CHUNK_SIZE]
// Each byte = BlockType (AIR, STONE, DIRT, GRASS, WATER, etc.)
```

### Proposed Hex Approach
Each hex column stores a **sparse vertical array** of voxel types:

```c
#define HEX_CHUNK_SIZE    32
#define HEX_CHUNK_LAYERS  128   // 128m vertical range per chunk

typedef struct HexChunk {
    int cx, cz;
    // 3D voxel data: [col][row][layer]
    uint8_t voxels[HEX_CHUNK_SIZE][HEX_CHUNK_SIZE][HEX_CHUNK_LAYERS];
    int base_layer;             // world layer offset (for chunks at different altitudes)
    bool dirty;
    // ... mesh data
} HexChunk;
```

**Memory per chunk:** 32 x 32 x 128 = 131,072 bytes (128 KB). Reasonable.

### Alternative: Run-Length Encoding (RLE) per column
For mostly-solid terrain, RLE compresses well:
```c
typedef struct {
    uint8_t type;     // AIR, STONE, etc.
    uint16_t count;   // number of consecutive layers
} VoxelRun;
```
A typical terrain column (50m stone + 78m air) = 2 runs = 6 bytes vs 128 bytes raw.

**Recommendation:** Start with raw 3D array (simple, fast random access). Optimize to RLE later if memory is a concern.

---

## 2. Raycast System (Tenebris DDA vs Current Step-March)

### Current (Hex-Planets): Step-Based
- March along ray in 0.1m steps
- At each step: project to tangent plane, convert to hex grid, check radius vs terrain height
- **Problem:** Only checks column top. Can't detect side faces or individual voxel layers.

### Tenebris: DDA (Digital Differential Analyzer)
```
1. Compute step direction per axis: step = sign(ray_dir)
2. Compute deltaDist = |1 / ray_dir| per axis (distance to cross one voxel)
3. Compute initial sideDist to first boundary per axis
4. Loop:
   a. Advance along axis with smallest sideDist
   b. Check if voxel at current position is solid
   c. If hit: return position + normal (derived from which axis crossed)
   d. If exceeded max distance: return miss
```

### Proposed: Hex DDA
Hex grids don't have orthogonal axes, so standard DDA needs adaptation:

**Option A: Tangent-plane 2D DDA + vertical stepping**
1. Project ray onto tangent plane (lx, lz) + radial component (ly)
2. March through hex cells using **hex line traversal** (cube coordinates)
3. At each hex cell, check all voxel layers for intersection with the vertical component

**Option B: Keep step-march, but check 3D voxels**
1. Same 0.1m step approach (already works for curved surface)
2. At each step: project to hex grid → get (col, row)
3. Compute altitude layer from radius
4. Check `voxels[col][row][layer]` for solidity
5. Compute face normal from which boundary was crossed (top, side, etc.)

**Recommendation:** Option B is simpler and already proven. Just extend the existing step-march to check layer + 6 neighbor faces.

### Face Normal Detection
When a hit is detected, determine which face was hit:
- **Top/bottom face:** Previous step was in a different layer (above/below)
- **Side face:** Previous step was in a different hex column
- This determines where a placed block goes (on top vs. on the side)

---

## 3. Block Breaking & Placing (Per-Voxel)

### Current: Column-Only
```c
// Break: decrement column height
chunk->heights[col][row]--;
// Place: increment column height
chunk->heights[col][row]++;
```

### Target: True 3D Editing
```c
// Break: set specific voxel to AIR
chunk->voxels[col][row][layer] = VOXEL_AIR;

// Place: set adjacent voxel to block type
// placement_pos = hit_pos + hit_normal (offset to empty cell)
chunk->voxels[place_col][place_row][place_layer] = selected_type;
```

### Placement Logic (from Tenebris)
1. Raycast → get hit position + surface normal
2. Placement position = `hit_cell + normal_offset`
   - Hit top face → place one layer above
   - Hit side face → place in adjacent hex column at same layer
   - Hit bottom face → place one layer below
3. **Player collision check:** Don't place block if it would intersect player AABB
4. Mark chunk dirty for remesh

### Mesh Rebuild
With 3D voxels, mesh generation changes:
- **Top face:** Only emit if voxel above is AIR
- **Bottom face:** Only emit if voxel below is AIR (for ceilings!)
- **Side faces:** Only emit if neighbor hex at same layer is AIR
- This naturally creates tunnels, caves, overhangs, floating structures

---

## 4. Collision System (The Big Upgrade)

### Current Problem
Player position snaps to `lod_tree_terrain_height()` which returns the highest solid point. Walking into a 100m cliff teleports you to the top. No wall blocking.

### Tenebris Collision (3-tier system)

**Tier 1: Step Height Check**
- Scan destination column for walkable floors (air-to-solid transitions)
- Only allow stepping UP by `AUTO_STEP_HEIGHT` (0.3m)
- Allow unlimited drops (gravity handles falling)

**Tier 2: Wall Collision**
- Check destination + all neighbor cells within PLAYER_RADIUS
- For each cell, scan entire voxel column
- Block is a "wall" if: `block_top > ground_height + AUTO_STEP_HEIGHT`
- If player AABB overlaps wall horizontally → reject movement

**Tier 3: Headroom Check**
- At destination, check 2m above feet for solid blocks
- Prevents walking into 1-block tunnels

### Proposed Hex Collision System

```
Parameters:
  PLAYER_RADIUS  = 0.3m
  PLAYER_HEIGHT  = 2.0m (need 2m clearance)
  STEP_HEIGHT    = 0.3m (auto-step up small ledges)
  EYE_HEIGHT     = 1.7m

Movement check(old_pos, new_pos):
  1. Find ground at new_pos:
     - Get hex column at new_pos
     - Scan voxels top-down to find highest solid with AIR above
     - This is "ground_layer"

  2. Step check:
     - If ground_layer > old_ground + STEP_HEIGHT → BLOCKED (wall too tall)
     - User should see the wall, not teleport on top

  3. Headroom check:
     - From ground_layer, check layers ground+1 through ground+PLAYER_HEIGHT
     - If ANY are solid → BLOCKED (ceiling too low)

  4. Neighbor wall check (within PLAYER_RADIUS):
     - For each of 6 hex neighbors at new_pos
     - Check if any solid voxel overlaps player's vertical extent
     - If horizontal distance < PLAYER_RADIUS → BLOCKED

  5. If all pass → allow movement, set feet to ground_layer
```

### The "2m Wall Block" Rule
> "any hex voxel that is 2m or higher should block my character from passing through"

Implementation:
```c
int height_diff = neighbor_ground - my_ground;
if (height_diff > STEP_HEIGHT_LAYERS) {
    // This is a wall, not a step. Block movement.
    // STEP_HEIGHT_LAYERS = 1 (since 1 layer ≈ 1m, and step = 0.3m)
    // Any 1+ layer difference blocks. For 2m rule: use threshold of 2
    return COLLISION_BLOCKED;
}
```

### Slide Along Walls
When blocked, don't just stop — project velocity along the wall:
1. Try full movement → blocked
2. Try X-only component → if clear, slide along X
3. Try Z-only component → if clear, slide along Z
4. Both blocked → stop completely

This gives smooth wall-sliding behavior (Tenebris does this).

---

## 5. Mesh Generation Changes for 3D Voxels

### Current: Heightfield Mesh
- Top cap at column height
- Side walls only where column is taller than neighbor
- No bottom faces (terrain is always solid below)

### Target: Full 3D Face Culling
For each solid voxel at (col, row, layer):

| Face | Condition to Emit |
|------|------------------|
| **Top** | `voxels[col][row][layer+1] == AIR` |
| **Bottom** | `voxels[col][row][layer-1] == AIR` |
| **Side 0-5** | `neighbor_hex_voxels[layer] == AIR` (6 hex neighbors) |

This automatically handles:
- **Tunnels:** Air voxels carved through mountain → ceiling faces + wall faces
- **Overhangs:** Solid blocks with air below → bottom face rendered
- **Floating structures:** All 8 faces (top + bottom + 6 sides) rendered
- **Caves:** Enclosed air pockets with visible interior surfaces

### Vertex Count Impact
Worst case (checkerboard pattern): every voxel exposes all faces.
Typical terrain: mostly interior voxels are culled (only surface matters).
A 32x32 chunk with 128 layers but only ~2-5 surface layers exposed = similar vertex count to current system.

---

## 6. LOD Integration for Structures

### The Arch Test Proved
- Coarse mesh (depths 0-12): can raise terrain height for structures
- Hex mesh (depth 13): can render hex prism structures with walls/caps
- Collision: `lod_tree_terrain_height()` can account for structures

### For True 3D Voxels at Distance
At depths 0-12 (far away), individual voxels aren't visible. Options:

**Option A: Heightfield approximation at distance**
- At coarse LOD, just use max height per column (current approach)
- Tunnels/caves invisible at distance (acceptable — Minecraft does this)
- Only render full 3D voxels at depth 13 (close range)

**Option B: Silhouette-preserving LOD**
- At medium distance, merge adjacent voxels into larger blocks
- Maintain overall structure shape but reduce triangle count
- More complex but better visual quality for large structures

**Recommendation:** Option A first. Minecraft proves you only need full voxel detail up close.

---

## 7. Implementation Roadmap

### Step 1: 3D Voxel Storage
- Replace `heights[]`/`types[]` with `voxels[col][row][layer]` in HexChunk
- Convert terrain generation to fill 3D array (solid below surface, air above)
- Update break/place to modify individual voxels

### Step 2: 3D Mesh Generation
- Rewrite hex mesh builder to emit faces based on neighbor air checks
- Top, bottom, and 6 side faces per voxel
- Only emit faces adjacent to AIR

### Step 3: 3D Raycast
- Extend step-march raycast to check voxel layers
- Return hit face normal (top/bottom/side)
- Place blocks adjacent to hit face

### Step 4: Wall Collision
- Replace ground-snap with proper collision:
  - Step height check (only step up ≤ 0.3m)
  - Wall blocking (≥ 2m height diff = wall)
  - Headroom check (need 2m clearance)
  - Wall sliding (project velocity along blocked axis)

### Step 5: LOD Integration
- At depth 13: full 3D voxel rendering
- At depths 0-12: heightfield approximation (max column height)
- Structures visible as raised terrain at distance

---

## Appendix: Tenebris Code Reference

| Feature | File | Lines |
|---------|------|-------|
| DDA Raycast | `src/world/World.ts` | 114-185 |
| Block Break/Place | `src/player/BlockInteraction.ts` | 116-156 |
| Collision (walls) | `src/player/PlanetPlayer.ts` | 1681-1897 |
| Collision (step) | `src/player/PlanetPlayer.ts` | 1637-1679 |
| Collision (headroom) | `src/player/PlanetPlayer.ts` | 1899+ |
| Chunk voxel data | `src/world/Chunk.ts` | 1-40 |
| Mesh generation | `src/world/Chunk.ts` | 118-200 |
| Wireframe highlight | `src/player/BlockInteraction.ts` | 27-34 |
| Hotbar / block select | `src/player/BlockInteraction.ts` | 50-80 |

| Feature | File (hex-planets) | Lines |
|---------|-------------------|-------|
| Hex raycast | `src/hex_terrain.c` | 1005-1073 |
| Hex break/place | `src/hex_terrain.c` | 1075-1108 |
| Hex highlight | `src/hex_terrain.c` | 1110-1152 |
| Ground collision | `src/camera.c` | 136-405 |
| LOD terrain height | `src/lod.c` | `lod_tree_terrain_height()` |
| Hex mesh gen | `src/lod.c` | `generate_hex_mesh_for_triangle()` |

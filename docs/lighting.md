# Voxel Lighting System — Research & Design

## Minecraft Lighting Reference

### Data Structure

Every block stores **two separate 4-bit values** (0-15), packed into one byte:
- **Sky light** (upper nibble): Light from the sky/sun
- **Block light** (lower nibble): Light from emitting blocks (torches, lava, etc.)

Final brightness at render time:
```
internal_light = max(sky_light_adjusted_for_time_of_day, block_light)
```

Minecraft uses a **16x16 lightmap texture** where U=block light, V=sky light. The GPU samples this to convert the two 0-15 values into final brightness+color. Block light is warm/orange, sky light is neutral white/blue.

### Sky Light Propagation

**Critical special rule**: Sky light at level 15 propagates downward through transparent blocks with **NO decay** (stays at 15). This creates infinite columns of full brightness from sky to first opaque block.

**Lateral and upward propagation**: Decays by -1 per block, same as block light.

**Once sky light drops below 15** (e.g. after one lateral step), it propagates identically to block light in all directions (-1 per block).

**Implementation:**
1. For each (X,Z) column, compute heightmap (Y of highest opaque block)
2. All transparent blocks above heightmap get sky light = 15
3. Seed BFS queue from level-15 blocks at heightmap boundary
4. BFS propagates: downward from 15 stays 15, all other directions -1

**Light-filtering blocks:**
- Air/glass: -1
- Water: -2 to -3
- Leaves: -1 to -2
- Tinted glass: blocks all light

### Block Light Propagation

Light-emitting blocks have fixed levels:
- Glowstone, Sea Lantern, Beacon: 15
- Torch, Fire: 14
- Lit Furnace: 13
- Soul Campfire: 10
- Redstone Torch: 7
- Magma Block: 3

Block light decays -1 per block in **all 6 directions** (no special vertical rule). Creates diamond-shaped illumination. A torch (14) reaches 13 blocks before going dark.

### BFS Flood-Fill Algorithm

**Propagation (adding light):**
```c
queue<LightNode> bfs_queue;

// Seed: set source block's light level, push to queue
set_light(x, y, z, emission_level);
bfs_queue.push({x, y, z});

while (!bfs_queue.empty()) {
    node = bfs_queue.pop();
    light_level = get_light(node.x, node.y, node.z);

    for each neighbor in 6 directions {
        if (neighbor.is_opaque) continue;

        neighbor_light = get_light(neighbor);
        new_level = light_level - 1;

        // Special: sky light going DOWN from 15 => no decay
        if (is_sky && light_level == 15 && dir == DOWN)
            new_level = 15;

        if (neighbor_light < new_level) {
            set_light(neighbor, new_level);
            bfs_queue.push(neighbor);
        }
    }
}
```

Each voxel visited exactly **once** — O(n) where n = illuminated blocks.

**Removal (two-queue algorithm):**
```c
queue<RemovalNode> removal_queue;   // {position, old_light_level}
queue<LightNode>   propagation_queue;

// Seed: record old level, zero it out
old_level = get_light(x, y, z);
set_light(x, y, z, 0);
removal_queue.push({x, y, z, old_level});

// Phase 1: Removal BFS
while (!removal_queue.empty()) {
    node = removal_queue.pop();
    for each neighbor in 6 directions {
        neighbor_light = get_light(neighbor);
        if (neighbor_light != 0 && neighbor_light < node.old_level) {
            // Lit by removed source => darken
            set_light(neighbor, 0);
            removal_queue.push({neighbor, neighbor_light});
        } else if (neighbor_light >= node.old_level) {
            // Lit by DIFFERENT source => re-propagate
            propagation_queue.push(neighbor);
        }
    }
}

// Phase 2: Re-propagation (fills gaps left by removal)
// Run standard propagation BFS using propagation_queue
```

Key insight: neighbors with lower light were dependent on removed source (darken them). Neighbors with equal/higher light are from other sources (become re-propagation seeds).

### Per-Face Shading

Fixed brightness multipliers by face normal (no actual sun direction):

| Face     | Multiplier |
|----------|-----------|
| Top (Y+) | 1.0       |
| Bottom (Y-) | 0.5    |
| North/South (Z) | 0.8 |
| East/West (X) | 0.6  |

### Ambient Occlusion (Per-Vertex)

For each vertex of an exposed face, check 3 neighboring blocks:
- **side1**: block adjacent along one edge
- **side2**: block adjacent along the other edge
- **corner**: block diagonally adjacent (touching vertex)

```c
int vertex_ao(bool side1, bool side2, bool corner) {
    if (side1 && side2) return 0;  // fully occluded
    return 3 - (side1 + side2 + corner);  // 0=dark, 3=bright
}
```

AO level → brightness: 3=1.0, 2=0.75, 1=0.5, 0=0.2

**Diagonal flip fix**: When a quad's 4 AO values differ, flip the triangle diagonal to align with AO gradient:
```c
if (a00 + a11 > a01 + a10)
    // flip diagonal
```

### Overhangs, Caves, Arches

This is where the sky light system shines:

1. Sky light at 15 propagates straight down to overhang top surface (blocked)
2. At overhang edges, sky light 15 continues down past the overhang
3. From those level-15 edge columns, **lateral BFS** carries light under the overhang, -1 per block
4. 5-block-deep overhang → back wall gets sky light 10
5. Cave entrance → light fades from 14, 13, 12... to 0 about 15 blocks in

This creates a **natural gradient** from bright at openings to dark deep inside.

### Combined Rendering Formula

```
vertex_brightness = texture_color
                  * face_direction_multiplier   // 0.5-1.0
                  * lightmap(block_light, sky_light)  // from 16x16 texture
                  * ao_multiplier               // 0.2-1.0
```

---

## Adaptation for Hex-Planets

### Differences from Minecraft

1. **Hex neighbors**: 6 lateral + 2 vertical = 8 neighbors (vs 6 cubic)
2. **Curved surface**: "down" = toward planet center, not -Y
3. **Per-face shading**: Can't use fixed axis multipliers. Options:
   - Use `dot(face_normal, sun_direction)` (current approach via slope normal)
   - Or use fixed multipliers mapped to face type: top=1.0, bottom=0.5, side=0.7
4. **AO**: Same `vertex_ao(side1, side2, corner)` formula, just identify which hex neighbors are side1/side2/corner for each vertex

### Current Implementation (depth-based, no flood-fill)

Sky light approximated by depth below natural terrain surface:
```c
sky_light = pow(0.8, depth_below_surface)
```
- Simple, no propagation needed
- Doesn't handle overhangs (blocks at surface under arch get full light)
- Doesn't handle lateral light entry into caves

### Future: Full BFS Flood-Fill

To properly light arches and caves:
1. Store per-voxel sky light (4 bits) in chunk data
2. Compute heightmap per column (outermost solid voxel)
3. Seed sky light = 15 above heightmap
4. BFS flood-fill: down from 15 = no decay, lateral = -1 per block
5. Re-run on chunk dirty (break/place)
6. Pass light level as vertex attribute (already have `sky_light` field in HexVertex)

**Performance**: BFS is O(n) per chunk update. With 32x32x128 chunks, worst case ~130K voxels but typical case much less (only surface + cave voxels). Can run on worker threads alongside mesh gen.

**Chunk boundaries**: Need to propagate light across chunk borders. When chunk A's light changes at boundary, mark chunk B dirty for re-lighting.

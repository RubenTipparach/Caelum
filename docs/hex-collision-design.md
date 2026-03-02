# Hex Terrain Collision Design

## Player Model
- Height: 2m (2 hex layers)
- Eye height: 1.7m above feet
- Radius: 0.3m (for future horizontal AABB)
- Collision is per-hex-column: check destination column's voxels before allowing movement

## Constants
- `AUTO_STEP_HEIGHT` = 1.0m — max height the player auto-steps onto without jumping
- `PLAYER_HEIGHT` = 2.0m — player capsule height (2 layers clearance needed)

## Collision Cases

Each diagram shows two adjacent hex columns (left = current, right = destination).
`P` = player body, `X` = solid voxel, `A` = air, `_` = ground below.

### Case 1: Full Wall — BLOCKED
```
  cur  dest
  [A]  [X]  ← head level: solid wall
  [P]  [X]  ← feet level: solid wall
  [_]  [_]
```
Height diff >= 2m > AUTO_STEP_HEIGHT → blocked.

### Case 2: Overhang — BLOCKED
```
  cur  dest
  [A]  [X]  ← solid at player head height
  [P]  [A]  ← air at feet level
  [_]  [_]
```
Ground at dest is same height, but headroom check fails:
solid voxel within 2 layers above dest ground → blocked.

### Case 3: Step Up (≤ 1m) — ALLOWED
```
  cur  dest
  [A]  [A]  ← 2m+ air above dest ground
  [A]  [A]
  [P]  [X]  ← dest ground is 1 layer higher
  [_]  [_]
```
Height diff = 1m ≤ AUTO_STEP_HEIGHT, 2+ layers air above → allowed.
Player auto-steps up smoothly.

### Case 4: High Step (> 1m) — BLOCKED
```
  cur  dest
  [A]  [A]
  [P]  [X]  ← dest ground is 2 layers higher
  [P]  [X]
  [_]  [_]
```
Height diff = 2m > AUTO_STEP_HEIGHT → blocked. Player must jump.

### Case 5: Walking Under Arch/Cave — ALLOWED
```
  cur  dest
  [X]  [X]  ← arch/cave ceiling (solid)
  [P]  [A]  ← air: inside arch
  [P]  [A]  ← air: inside arch
  [_]  [_]  ← ground (solid)
```
Dest ground = ground below arch (scan down from player altitude).
Height diff = 0, headroom = 2 layers of air → allowed.

**Key implementation detail:** ground height must be found relative to the
player's current altitude, not the absolute topmost surface. Otherwise
`ground_height` returns the arch roof and the player can't enter.

### Case 6: Low Ceiling — BLOCKED
```
  cur  dest
  [A]  [X]  ← ceiling 1 layer above ground
  [P]  [A]  ← only 1m clearance
  [_]  [_]
```
Headroom check: only 1 layer of air (< PLAYER_HEIGHT = 2) → blocked.

### Case 7: Falling Off Edge — ALLOWED
```
  cur  dest
  [P]  [A]
  [P]  [A]
  [_]  [A]  ← no ground nearby = cliff edge
       [_]
```
Height diff < 0 (lower ground) → always allowed. Gravity handles the fall.

### Case 8: Step Down (any height) — ALLOWED
```
  cur  dest
  [A]  [A]  ← plenty of air
  [P]  [A]
  [_]  [A]
       [_]  ← ground 1-2 layers lower
```
Height diff < 0 → allowed. Gravity snaps player to new ground.

## Algorithm

```
For each movement attempt (full, forward-only, side-only):
  1. Compute destination hex column (gcol, grow)
  2. Find dest ground: scan DOWN from player's current feet layer + step_height
     to find topmost solid-with-air-above at or below that search ceiling
  3. Step height check:
     - height_diff = dest_ground - cur_ground
     - if height_diff > AUTO_STEP_HEIGHT (1m): BLOCKED
     - if height_diff < 0: always allowed (falling)
  4. Headroom check:
     - Check PLAYER_HEIGHT (2) layers of air above dest ground
     - Any solid in that range: BLOCKED
  5. If all checks pass: MOVE

Wall sliding:
  - Try full movement → blocked?
  - Try forward component only → blocked?
  - Try side component only → blocked?
  - All blocked → stop
```

## Ground Height Search (Arch-Aware)

The naive approach (scan from chunk top) returns the arch roof when the player
is inside the arch. Fix: limit the scan to start from the player's current
altitude + step margin.

```c
// Find ground at or below max_world_layer
float hex_terrain_ground_height_below(ht, gcol, grow, max_world_layer);
```

This scans down from `max_world_layer` and returns the first solid with air
above, which correctly finds the floor inside arches/caves.

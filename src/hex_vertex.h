#ifndef HEX_VERTEX_H
#define HEX_VERTEX_H

// Vertex format for textured hex terrain meshes.
// Separate from LodVertex (which uses vertex colors) so the LOD system stays unchanged.
typedef struct HexVertex {
    float pos[3];       // World position (relative to floating origin)
    float normal[3];    // Surface normal
    float uv[2];        // Texture atlas UV
    float color[3];     // Terrain color (for distance fade to LOD vertex color)
    float sky_light;    // Depth-based sky light [0..1] (tenebris: 0.8^depth, 1.0 = surface)
} HexVertex;

#endif

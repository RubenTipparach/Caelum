#ifndef HEX_VERTEX_H
#define HEX_VERTEX_H

// Vertex format for textured hex terrain meshes.
// Separate from LodVertex (which uses vertex colors) so the LOD system stays unchanged.
typedef struct HexVertex {
    float pos[3];       // World position (relative to floating origin)
    float normal[3];    // Surface normal
    float uv[2];        // Texture atlas UV
} HexVertex;

#endif

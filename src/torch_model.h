// torch_model.h -- Procedural torch geometry (header-only, data generation)
// Generates a torch model consisting of:
//   1. Brown cylinder handle (h=0.4, r=0.03, 8 segments)
//   2. Dark gray head (h=0.1, r=0.06, 8 segments)
//   3. Orange cone flame (h=0.15, r=0.05, 8 segments)
//
// Vertex format: float3 pos + float3 normal + float3 color + float anim_weight = 10 floats
// Total: 3 parts * 8 segments * 2 triangles * 3 verts = 144 tris (side faces)
//        + 2 caps * 8 triangles * 3 verts = 48 tris (top caps)
//        = ~192 triangles total (simplified: only side faces + top caps)

#ifndef TORCH_MODEL_H
#define TORCH_MODEL_H

#include <math.h>

#define TORCH_VERT_FLOATS  10   // pos(3) + normal(3) + color(3) + anim_weight(1)
#define TORCH_SEGMENTS      8

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Cylinder section: generates side faces + top cap for one section of the torch.
// origin_y = bottom of section, height = section height, radius = section radius
// color_r/g/b = vertex color, anim = animation weight (0 or 1)
// out_verts: pointer to float array, out_count: pointer to vertex count (in vertices)
static inline void torch_gen_section(float* out_verts, int* out_count,
                                     float origin_y, float height, float radius,
                                     float top_radius,
                                     float color_r, float color_g, float color_b,
                                     float anim) {
    int vi = *out_count * TORCH_VERT_FLOATS;
    float bot_y = origin_y;
    float top_y = origin_y + height;

    // Precompute circle positions
    float cx[TORCH_SEGMENTS], cz[TORCH_SEGMENTS];
    for (int i = 0; i < TORCH_SEGMENTS; i++) {
        float angle = (float)(2.0 * M_PI * i / TORCH_SEGMENTS);
        cx[i] = cosf(angle);
        cz[i] = sinf(angle);
    }

    // Side faces (2 triangles per segment)
    for (int i = 0; i < TORCH_SEGMENTS; i++) {
        int j = (i + 1) % TORCH_SEGMENTS;

        float bx0 = cx[i] * radius,    bz0 = cz[i] * radius;
        float bx1 = cx[j] * radius,    bz1 = cz[j] * radius;
        float tx0 = cx[i] * top_radius, tz0 = cz[i] * top_radius;
        float tx1 = cx[j] * top_radius, tz1 = cz[j] * top_radius;

        // Normal: average of the two corner normals (good enough for 8 segments)
        float nx = (cx[i] + cx[j]) * 0.5f;
        float nz = (cz[i] + cz[j]) * 0.5f;
        float nl = sqrtf(nx * nx + nz * nz);
        if (nl > 0.001f) { nx /= nl; nz /= nl; }

        // Triangle 1: bot0, top1, bot1 (CCW when viewed from outside)
        #define EMIT_VERT(px, py, pz, nnx, nny, nnz, aw) do { \
            out_verts[vi++] = (px); out_verts[vi++] = (py); out_verts[vi++] = (pz); \
            out_verts[vi++] = (nnx); out_verts[vi++] = (nny); out_verts[vi++] = (nnz); \
            out_verts[vi++] = color_r; out_verts[vi++] = color_g; out_verts[vi++] = color_b; \
            out_verts[vi++] = (aw); \
            (*out_count)++; \
        } while(0)

        EMIT_VERT(bx0, bot_y, bz0,  nx, 0.0f, nz, anim);
        EMIT_VERT(tx1, top_y, tz1,  nx, 0.0f, nz, anim);
        EMIT_VERT(bx1, bot_y, bz1,  nx, 0.0f, nz, anim);

        // Triangle 2: bot0, top0, top1
        EMIT_VERT(bx0, bot_y, bz0,  nx, 0.0f, nz, anim);
        EMIT_VERT(tx0, top_y, tz0,  nx, 0.0f, nz, anim);
        EMIT_VERT(tx1, top_y, tz1,  nx, 0.0f, nz, anim);
    }

    // Top cap (fan from center)
    for (int i = 0; i < TORCH_SEGMENTS; i++) {
        int j = (i + 1) % TORCH_SEGMENTS;
        EMIT_VERT(0.0f, top_y, 0.0f,  0.0f, 1.0f, 0.0f, anim);
        EMIT_VERT(cx[i] * top_radius, top_y, cz[i] * top_radius,  0.0f, 1.0f, 0.0f, anim);
        EMIT_VERT(cx[j] * top_radius, top_y, cz[j] * top_radius,  0.0f, 1.0f, 0.0f, anim);
    }

    #undef EMIT_VERT
}

// Returns the total number of vertices in the torch model
static inline int torch_model_vertex_count(void) {
    // Each section: 8 segments * 2 side tris * 3 verts + 8 cap tris * 3 verts = 72 verts
    // 3 sections = 216 verts
    return 3 * (TORCH_SEGMENTS * 2 * 3 + TORCH_SEGMENTS * 3);
}

// Generate the complete torch model into out_verts buffer.
// Buffer must hold at least torch_model_vertex_count() * TORCH_VERT_FLOATS floats.
// Model is centered at (0,0,0) with Y-up, handle starts at y=0.
static inline void torch_model_generate(float* out_verts) {
    int count = 0;

    // 1. Handle: brown cylinder (y=0 to y=0.4, r=0.03)
    torch_gen_section(out_verts, &count,
                      0.0f, 0.4f, 0.03f, 0.03f,
                      0.45f, 0.28f, 0.12f,  // brown
                      0.0f);                  // no animation

    // 2. Head: dark gray block (y=0.4 to y=0.5, r=0.06)
    torch_gen_section(out_verts, &count,
                      0.4f, 0.1f, 0.06f, 0.06f,
                      0.3f, 0.3f, 0.3f,     // dark gray
                      0.0f);                  // no animation

    // 3. Flame: orange cone (y=0.5 to y=0.65, r=0.05 → 0.0 at tip)
    torch_gen_section(out_verts, &count,
                      0.5f, 0.15f, 0.05f, 0.0f,
                      1.0f, 0.6f, 0.1f,     // bright orange
                      1.0f);                  // full animation
}

#endif // TORCH_MODEL_H

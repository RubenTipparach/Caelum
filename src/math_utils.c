#include "math_utils.h"
#include <math.h>

#define PHI 1.6180339887498949f

// 12 icosahedron vertices (will be normalized at runtime)
const HMM_Vec3 ICO_VERTICES[ICO_VERTEX_COUNT] = {
    {{ 0.0f,  1.0f,  PHI}},
    {{ 0.0f, -1.0f,  PHI}},
    {{ PHI,   0.0f,  1.0f}},
    {{ 1.0f,  PHI,   0.0f}},
    {{-1.0f,  PHI,   0.0f}},
    {{-PHI,   0.0f,  1.0f}},
    {{ 0.0f,  1.0f, -PHI}},
    {{ PHI,   0.0f, -1.0f}},
    {{ 1.0f, -PHI,   0.0f}},
    {{-1.0f, -PHI,   0.0f}},
    {{-PHI,   0.0f, -1.0f}},
    {{ 0.0f, -1.0f, -PHI}},
};

// 20 triangular faces (vertex indices, counter-clockwise winding from outside)
const int ICO_FACES[ICO_FACE_COUNT][3] = {
    { 0,  1,  2}, { 0,  2,  3}, { 0,  3,  4}, { 0,  4,  5}, { 0,  5,  1},
    { 1,  8,  2}, { 2,  7,  3}, { 3,  6,  4}, { 4, 10,  5}, { 5,  9,  1},
    { 1,  9,  8}, { 2,  8,  7}, { 3,  7,  6}, { 4,  6, 10}, { 5, 10,  9},
    {11,  8,  9}, {11,  7,  8}, {11,  6,  7}, {11, 10,  6}, {11,  9, 10},
};

HMM_Vec3 vec3_normalize(HMM_Vec3 v) {
    float len = sqrtf(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
    if (len < 1e-10f) return (HMM_Vec3){{0, 0, 0}};
    float inv = 1.0f / len;
    return (HMM_Vec3){{v.X * inv, v.Y * inv, v.Z * inv}};
}

float vec3_dot(HMM_Vec3 a, HMM_Vec3 b) {
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}

HMM_Vec3 vec3_cross(HMM_Vec3 a, HMM_Vec3 b) {
    return (HMM_Vec3){{
        a.Y * b.Z - a.Z * b.Y,
        a.Z * b.X - a.X * b.Z,
        a.X * b.Y - a.Y * b.X
    }};
}

HMM_Vec3 vec3_lerp(HMM_Vec3 a, HMM_Vec3 b, float t) {
    return (HMM_Vec3){{
        a.X + (b.X - a.X) * t,
        a.Y + (b.Y - a.Y) * t,
        a.Z + (b.Z - a.Z) * t
    }};
}

HMM_Vec3 vec3_scale(HMM_Vec3 v, float s) {
    return (HMM_Vec3){{v.X * s, v.Y * s, v.Z * s}};
}

HMM_Vec3 vec3_add(HMM_Vec3 a, HMM_Vec3 b) {
    return (HMM_Vec3){{a.X + b.X, a.Y + b.Y, a.Z + b.Z}};
}

HMM_Vec3 vec3_sub(HMM_Vec3 a, HMM_Vec3 b) {
    return (HMM_Vec3){{a.X - b.X, a.Y - b.Y, a.Z - b.Z}};
}

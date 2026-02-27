#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "HandmadeMath.h"

// Icosahedron geometry
#define ICO_VERTEX_COUNT 12
#define ICO_FACE_COUNT 20

extern const HMM_Vec3 ICO_VERTICES[ICO_VERTEX_COUNT];
extern const int ICO_FACES[ICO_FACE_COUNT][3];

HMM_Vec3 vec3_normalize(HMM_Vec3 v);
float vec3_dot(HMM_Vec3 a, HMM_Vec3 b);
HMM_Vec3 vec3_cross(HMM_Vec3 a, HMM_Vec3 b);
HMM_Vec3 vec3_lerp(HMM_Vec3 a, HMM_Vec3 b, float t);
HMM_Vec3 vec3_scale(HMM_Vec3 v, float s);
HMM_Vec3 vec3_add(HMM_Vec3 a, HMM_Vec3 b);
HMM_Vec3 vec3_sub(HMM_Vec3 a, HMM_Vec3 b);

#endif

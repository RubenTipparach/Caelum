// sokol_impl.c — Single compilation unit for all Sokol implementations.
// Platform backend is set via CMake define (SOKOL_D3D11, SOKOL_METAL, etc.)

#define SOKOL_IMPL
#include "sokol_app.h"
#include "sokol_gfx.h"
#include "sokol_glue.h"
#include "sokol_log.h"
#include "sokol_time.h"

#define SOKOL_DEBUGTEXT_IMPL
#include "util/sokol_debugtext.h"

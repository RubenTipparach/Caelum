// torch.c -- Torch rendering system implementation
#include "torch.h"
#include "torch_model.h"
#include "torch.glsl.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void torch_init(TorchSystem* ts) {
    memset(ts, 0, sizeof(*ts));

    // Generate procedural torch geometry
    ts->model_vertex_count = torch_model_vertex_count();
    int float_count = ts->model_vertex_count * TORCH_VERT_FLOATS;
    float* verts = (float*)malloc(float_count * sizeof(float));
    torch_model_generate(verts);

    ts->model_buf = sg_make_buffer(&(sg_buffer_desc){
        .usage.vertex_buffer = true,
        .data = { verts, float_count * sizeof(float) },
        .label = "torch-model-vbuf",
    });
    free(verts);

    // Create pipeline
    sg_shader shd = sg_make_shader(torch_shader_desc(sg_query_backend()));
    ts->pip = sg_make_pipeline(&(sg_pipeline_desc){
        .shader = shd,
        .layout = {
            .attrs = {
                [ATTR_torch_a_position]    = { .format = SG_VERTEXFORMAT_FLOAT3 },
                [ATTR_torch_a_normal]      = { .format = SG_VERTEXFORMAT_FLOAT3 },
                [ATTR_torch_a_color]       = { .format = SG_VERTEXFORMAT_FLOAT3 },
                [ATTR_torch_a_anim_weight] = { .format = SG_VERTEXFORMAT_FLOAT },
            }
        },
        .depth = {
            .compare = SG_COMPAREFUNC_LESS_EQUAL,
            .write_enabled = true,
        },
        .cull_mode = SG_CULLMODE_BACK,
        .face_winding = SG_FACEWINDING_CCW,
        .label = "torch-pipeline",
    });

    printf("[TORCH] Initialized: %d model verts\n", ts->model_vertex_count);
    fflush(stdout);
}

void torch_destroy(TorchSystem* ts) {
    if (ts->model_buf.id != SG_INVALID_ID) {
        sg_destroy_buffer(ts->model_buf);
    }
    if (ts->pip.id != SG_INVALID_ID) {
        sg_destroy_pipeline(ts->pip);
    }
    memset(ts, 0, sizeof(*ts));
}

void torch_update(TorchSystem* ts, float dt) {
    for (int i = 0; i < TORCH_MAX_INSTANCES; i++) {
        if (ts->instances[i].active) {
            ts->instances[i].time += dt;
        }
    }
}

int torch_add(TorchSystem* ts, float world_x, float world_y, float world_z,
              float up_x, float up_y, float up_z,
              int chunk_cx, int chunk_cz, int gcol, int grow, int layer) {
    // Find free slot
    for (int i = 0; i < TORCH_MAX_INSTANCES; i++) {
        if (!ts->instances[i].active) {
            TorchInstance* inst = &ts->instances[i];
            inst->active = true;
            inst->world_pos[0] = world_x;
            inst->world_pos[1] = world_y;
            inst->world_pos[2] = world_z;
            inst->up_dir[0] = up_x;
            inst->up_dir[1] = up_y;
            inst->up_dir[2] = up_z;
            inst->chunk_cx = chunk_cx;
            inst->chunk_cz = chunk_cz;
            inst->gcol = gcol;
            inst->grow = grow;
            inst->layer = layer;
            inst->time = 0.0f;
            ts->instance_count++;
            return i;
        }
    }
    return -1;  // full
}

bool torch_remove(TorchSystem* ts, int gcol, int grow, int layer) {
    for (int i = 0; i < TORCH_MAX_INSTANCES; i++) {
        if (ts->instances[i].active &&
            ts->instances[i].gcol == gcol &&
            ts->instances[i].grow == grow &&
            ts->instances[i].layer == layer) {
            ts->instances[i].active = false;
            ts->instance_count--;
            return true;
        }
    }
    return false;
}

void torch_remove_chunk(TorchSystem* ts, int chunk_cx, int chunk_cz) {
    for (int i = 0; i < TORCH_MAX_INSTANCES; i++) {
        if (ts->instances[i].active &&
            ts->instances[i].chunk_cx == chunk_cx &&
            ts->instances[i].chunk_cz == chunk_cz) {
            ts->instances[i].active = false;
            ts->instance_count--;
        }
    }
}

void torch_render(TorchSystem* ts, HMM_Mat4 vp,
                  const double camera_offset[3], const double camera_offset_low[3],
                  const double world_origin[3],
                  float fcoef, float far_plane, float z_bias,
                  HMM_Vec3 sun_dir) {
    if (ts->instance_count <= 0) return;

    sg_apply_pipeline(ts->pip);

    sg_bindings bind = {
        .vertex_buffers[0] = ts->model_buf,
    };
    sg_apply_bindings(&bind);

    // Fragment shader uniforms (same for all instances)
    torch_fs_params_t fs_p = {
        .sun_direction = {{ sun_dir.X, sun_dir.Y, sun_dir.Z, 0.0f }},
    };
    sg_apply_uniforms(UB_torch_fs_params, &SG_RANGE(fs_p));

    // Draw each torch instance
    for (int i = 0; i < TORCH_MAX_INSTANCES; i++) {
        if (!ts->instances[i].active) continue;

        TorchInstance* inst = &ts->instances[i];

        torch_vs_params_t vs_p = {
            .vp = vp,
            .camera_offset = {{
                (float)camera_offset[0],
                (float)camera_offset[1],
                (float)camera_offset[2],
                0.0f
            }},
            .camera_offset_low = {{
                (float)camera_offset_low[0],
                (float)camera_offset_low[1],
                (float)camera_offset_low[2],
                0.0f
            }},
            .log_depth = {{ fcoef, far_plane, z_bias, 0.0f }},
            .torch_pos = {{
                inst->world_pos[0] - (float)world_origin[0],
                inst->world_pos[1] - (float)world_origin[1],
                inst->world_pos[2] - (float)world_origin[2],
                inst->time
            }},
            .torch_up = {{
                inst->up_dir[0],
                inst->up_dir[1],
                inst->up_dir[2],
                0.3f  // flicker amount
            }},
        };
        sg_apply_uniforms(UB_torch_vs_params, &SG_RANGE(vs_p));
        sg_draw(0, ts->model_vertex_count, 1);
    }
}

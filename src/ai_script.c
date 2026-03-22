#include "ai_script.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Simple JSON string finder (reuse pattern from ai_npc.c)
static bool jfind_str(const char* json, const char* key, char* out, int out_size) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    if (*p != '"') return false;
    p++;
    int i = 0;
    while (*p && *p != '"' && i < out_size - 1) {
        if (*p == '\\' && *(p+1)) { p++; }
        out[i++] = *p++;
    }
    out[i] = '\0';
    return true;
}

static bool jfind_int(const char* json, const char* key, int* out) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    *out = atoi(p);
    return true;
}

static bool jfind_float(const char* json, const char* key, float* out) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    *out = (float)atof(p);
    return true;
}

static ScriptActionType parse_action_type(const char* name) {
    if (strcmp(name, "move_to") == 0)       return SCRIPT_ACT_MOVE_TO;
    if (strcmp(name, "place") == 0)          return SCRIPT_ACT_PLACE;
    if (strcmp(name, "break") == 0)          return SCRIPT_ACT_BREAK;
    if (strcmp(name, "say") == 0)            return SCRIPT_ACT_SAY;
    if (strcmp(name, "wait") == 0)           return SCRIPT_ACT_WAIT;
    if (strcmp(name, "scan") == 0)           return SCRIPT_ACT_SCAN;
    if (strcmp(name, "scan_and_report") == 0) return SCRIPT_ACT_SCAN_REPORT;
    if (strcmp(name, "look_at") == 0)        return SCRIPT_ACT_LOOK_AT;
    if (strcmp(name, "jump") == 0)           return SCRIPT_ACT_JUMP;
    if (strcmp(name, "loop") == 0)           return SCRIPT_ACT_LOOP;
    if (strcmp(name, "move_rel") == 0 || strcmp(name, "move_relative") == 0) return SCRIPT_ACT_MOVE_REL;
    if (strcmp(name, "place_rel") == 0 || strcmp(name, "place_relative") == 0) return SCRIPT_ACT_PLACE_REL;
    if (strcmp(name, "break_rel") == 0 || strcmp(name, "break_relative") == 0) return SCRIPT_ACT_BREAK_REL;
    if (strcmp(name, "set_anchor") == 0) return SCRIPT_ACT_SET_ANCHOR;
    return SCRIPT_ACT_NONE;
}

bool ai_script_parse(const char* json, AiScript* out) {
    memset(out, 0, sizeof(AiScript));

    jfind_str(json, "name", out->name, sizeof(out->name));
    jfind_str(json, "description", out->description, sizeof(out->description));

    // Find "steps" array
    const char* steps = strstr(json, "\"steps\"");
    if (!steps) return false;
    const char* bracket = strchr(steps, '[');
    if (!bracket) return false;

    // Parse each step object
    const char* cursor = bracket + 1;
    while (out->step_count < AI_SCRIPT_MAX_STEPS) {
        const char* obj_start = strchr(cursor, '{');
        if (!obj_start) break;

        // Find matching closing brace (simple: no nested objects)
        const char* obj_end = strchr(obj_start, '}');
        if (!obj_end) break;

        int obj_len = (int)(obj_end - obj_start + 1);
        char obj[512];
        if (obj_len >= (int)sizeof(obj)) obj_len = (int)sizeof(obj) - 1;
        memcpy(obj, obj_start, obj_len);
        obj[obj_len] = '\0';

        AiScriptStep* step = &out->steps[out->step_count];
        memset(step, 0, sizeof(AiScriptStep));

        char action_name[32];
        if (jfind_str(obj, "action", action_name, sizeof(action_name))) {
            step->action = parse_action_type(action_name);
            jfind_int(obj, "q", &step->q);
            jfind_int(obj, "r", &step->r);
            jfind_int(obj, "layer", &step->layer);
            jfind_int(obj, "dq", &step->dq);
            jfind_int(obj, "dr", &step->dr);
            jfind_int(obj, "dlayer", &step->dlayer);
            jfind_str(obj, "block", step->block, sizeof(step->block));
            jfind_str(obj, "text", step->text, sizeof(step->text));
            jfind_float(obj, "seconds", &step->seconds);
            jfind_int(obj, "radius", &step->radius);
            jfind_int(obj, "from", &step->from_step);
            jfind_int(obj, "count", &step->loop_count);
            step->_loop_remaining = step->loop_count;
            out->step_count++;
        }

        cursor = obj_end + 1;
    }

    printf("[SCRIPT] Parsed '%s': %d steps\n", out->name, out->step_count);
    fflush(stdout);
    return out->step_count > 0;
}

bool ai_script_load(const char* path, AiScript* out) {
    FILE* f = fopen(path, "rb");
    if (!f) return false;

    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);

    if (size <= 0 || size > 65536) { fclose(f); return false; }

    char* buf = (char*)malloc(size + 1);
    fread(buf, 1, size, f);
    buf[size] = '\0';
    fclose(f);

    bool ok = ai_script_parse(buf, out);
    free(buf);

    if (ok) {
        printf("[SCRIPT] Loaded '%s' from %s\n", out->name, path);
    } else {
        printf("[SCRIPT] Failed to load from %s\n", path);
    }
    fflush(stdout);
    return ok;
}

bool ai_script_save(const AiScript* script, const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) return false;

    fprintf(f, "{\n  \"name\": \"%s\",\n  \"description\": \"%s\",\n  \"steps\": [\n",
            script->name, script->description);

    for (int i = 0; i < script->step_count; i++) {
        const AiScriptStep* s = &script->steps[i];
        fprintf(f, "    {\"action\": \"");
        switch (s->action) {
            case SCRIPT_ACT_MOVE_TO: fprintf(f, "move_to\", \"q\": %d, \"r\": %d", s->q, s->r); break;
            case SCRIPT_ACT_PLACE: fprintf(f, "place\", \"q\": %d, \"r\": %d, \"layer\": %d, \"block\": \"%s\"", s->q, s->r, s->layer, s->block); break;
            case SCRIPT_ACT_BREAK: fprintf(f, "break\", \"q\": %d, \"r\": %d, \"layer\": %d", s->q, s->r, s->layer); break;
            case SCRIPT_ACT_SAY: fprintf(f, "say\", \"text\": \"%s\"", s->text); break;
            case SCRIPT_ACT_WAIT: fprintf(f, "wait\", \"seconds\": %.1f", s->seconds); break;
            case SCRIPT_ACT_SCAN: fprintf(f, "scan\", \"radius\": %d", s->radius); break;
            case SCRIPT_ACT_SCAN_REPORT: fprintf(f, "scan_and_report\""); break;
            case SCRIPT_ACT_LOOK_AT: fprintf(f, "look_at\", \"q\": %d, \"r\": %d", s->q, s->r); break;
            case SCRIPT_ACT_JUMP: fprintf(f, "jump\""); break;
            case SCRIPT_ACT_LOOP: fprintf(f, "loop\", \"from\": %d, \"count\": %d", s->from_step, s->loop_count); break;
            case SCRIPT_ACT_MOVE_REL: fprintf(f, "move_rel\", \"dq\": %d, \"dr\": %d", s->dq, s->dr); break;
            case SCRIPT_ACT_PLACE_REL: fprintf(f, "place_rel\", \"dq\": %d, \"dr\": %d, \"dlayer\": %d, \"block\": \"%s\"", s->dq, s->dr, s->dlayer, s->block); break;
            case SCRIPT_ACT_BREAK_REL: fprintf(f, "break_rel\", \"dq\": %d, \"dr\": %d, \"dlayer\": %d", s->dq, s->dr, s->dlayer); break;
            case SCRIPT_ACT_SET_ANCHOR: fprintf(f, "set_anchor\", \"q\": %d, \"r\": %d, \"layer\": %d", s->q, s->r, s->layer); break;
            default: fprintf(f, "none\""); break;
        }
        fprintf(f, "}%s\n", i < script->step_count - 1 ? "," : "");
    }

    fprintf(f, "  ]\n}\n");
    fclose(f);

    printf("[SCRIPT] Saved '%s' to %s (%d steps)\n", script->name, path, script->step_count);
    fflush(stdout);
    return true;
}

void ai_script_runner_init(AiScriptRunner* runner) {
    memset(runner, 0, sizeof(AiScriptRunner));
    runner->state = SCRIPT_IDLE;
}

void ai_script_run(AiScriptRunner* runner, const AiScript* script) {
    runner->script = *script;
    runner->state = SCRIPT_RUNNING;
    runner->current_step = 0;
    runner->wait_timer = 0.0f;
    runner->error[0] = '\0';
    // Anchor will be set by agent update (current hex position) or by set_anchor step
    runner->anchor_set = false;

    // Reset loop counters
    for (int i = 0; i < runner->script.step_count; i++) {
        runner->script.steps[i]._loop_remaining = runner->script.steps[i].loop_count;
    }

    printf("[SCRIPT] Running '%s' (%d steps)\n", script->name, script->step_count);
    fflush(stdout);
}

void ai_script_stop(AiScriptRunner* runner) {
    runner->state = SCRIPT_IDLE;
    runner->current_step = 0;
    printf("[SCRIPT] Stopped\n");
    fflush(stdout);
}

ScriptActionType ai_script_update(AiScriptRunner* runner, float dt) {
    if (runner->state != SCRIPT_RUNNING && runner->state != SCRIPT_WAITING)
        return SCRIPT_ACT_NONE;

    // Waiting on a timer
    if (runner->state == SCRIPT_WAITING) {
        runner->wait_timer -= dt;
        if (runner->wait_timer > 0.0f) return SCRIPT_ACT_NONE;
        runner->state = SCRIPT_RUNNING;
        runner->current_step++;
    }

    // Check if done
    if (runner->current_step >= runner->script.step_count) {
        runner->state = SCRIPT_DONE;
        printf("[SCRIPT] '%s' completed\n", runner->script.name);
        fflush(stdout);
        return SCRIPT_ACT_NONE;
    }

    AiScriptStep* step = &runner->script.steps[runner->current_step];

    switch (step->action) {
        case SCRIPT_ACT_WAIT:
            runner->state = SCRIPT_WAITING;
            runner->wait_timer = step->seconds;
            return SCRIPT_ACT_NONE;

        case SCRIPT_ACT_LOOP:
            if (step->_loop_remaining > 0) {
                step->_loop_remaining--;
                if (step->from_step >= 0 && step->from_step < runner->script.step_count) {
                    runner->current_step = step->from_step;
                    return SCRIPT_ACT_NONE;
                }
            }
            // Loop exhausted, continue
            runner->current_step++;
            return SCRIPT_ACT_NONE;

        case SCRIPT_ACT_SAY:
            printf("[SCRIPT] Say: \"%s\"\n", step->text);
            fflush(stdout);
            // Fall through — caller handles the chat log push
            return step->action;

        case SCRIPT_ACT_SET_ANCHOR:
            // Explicit anchor — use step coords if provided, else agent sets it
            if (step->q != 0 || step->r != 0) {
                runner->anchor_q = step->q;
                runner->anchor_r = step->r;
                runner->anchor_layer = step->layer;
            }
            runner->anchor_set = true;
            printf("[SCRIPT] Anchor set: q=%d r=%d layer=%d\n",
                   runner->anchor_q, runner->anchor_r, runner->anchor_layer);
            fflush(stdout);
            runner->current_step++;
            return SCRIPT_ACT_NONE;

        case SCRIPT_ACT_MOVE_REL:
        case SCRIPT_ACT_PLACE_REL:
        case SCRIPT_ACT_BREAK_REL:
            // Resolve relative to anchor, then return as absolute for agent to handle
            // Agent update will read dq/dr/dlayer and add to anchor
            return step->action;

        case SCRIPT_ACT_MOVE_TO:
        case SCRIPT_ACT_PLACE:
        case SCRIPT_ACT_BREAK:
        case SCRIPT_ACT_SCAN:
        case SCRIPT_ACT_SCAN_REPORT:
        case SCRIPT_ACT_LOOK_AT:
        case SCRIPT_ACT_JUMP:
            // These need external handling by the agent
            return step->action;

        default:
            runner->current_step++;
            return SCRIPT_ACT_NONE;
    }
}

const AiScriptStep* ai_script_current_step(const AiScriptRunner* runner) {
    if (runner->state != SCRIPT_RUNNING && runner->state != SCRIPT_MOVING &&
        runner->state != SCRIPT_BUILDING)
        return NULL;
    if (runner->current_step >= runner->script.step_count) return NULL;
    return &runner->script.steps[runner->current_step];
}

void ai_script_advance(AiScriptRunner* runner) {
    if (runner->state == SCRIPT_RUNNING || runner->state == SCRIPT_MOVING) {
        runner->current_step++;
        runner->state = SCRIPT_RUNNING;
        if (runner->current_step >= runner->script.step_count) {
            runner->state = SCRIPT_DONE;
            printf("[SCRIPT] '%s' completed\n", runner->script.name);
            fflush(stdout);
        }
    }
}

bool ai_script_idle(const AiScriptRunner* runner) {
    return runner->state == SCRIPT_IDLE || runner->state == SCRIPT_DONE;
}

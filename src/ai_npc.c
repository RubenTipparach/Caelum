#include "ai_npc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ---- Minimal JSON helpers (same style as lobby.c) ----

static bool json_find_string(const char* json, const char* key, char* out, int out_size) {
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
        if (*p == '\\' && *(p+1)) { p++; } // skip escape
        out[i++] = *p++;
    }
    out[i] = '\0';
    return true;
}

static bool json_find_int(const char* json, const char* key, int* out) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    *out = atoi(p);
    return true;
}

// ---- JSON string escaping ----

static int json_escape(const char* src, char* dst, int dst_size) {
    int j = 0;
    for (int i = 0; src[i] && j < dst_size - 2; i++) {
        char c = src[i];
        if (c == '"' || c == '\\') {
            if (j + 2 >= dst_size) break;
            dst[j++] = '\\';
            dst[j++] = c;
        } else if (c == '\n') {
            if (j + 2 >= dst_size) break;
            dst[j++] = '\\';
            dst[j++] = 'n';
        } else if (c == '\r') {
            // skip
        } else if (c == '\t') {
            if (j + 2 >= dst_size) break;
            dst[j++] = '\\';
            dst[j++] = 't';
        } else {
            dst[j++] = c;
        }
    }
    dst[j] = '\0';
    return j;
}

// ---- Chat logging (JSONL) ----

static void ai_log_message(const char* agent_name, const char* type, const char* text) {
    // Write to cache/agents/{name}/chat_log.jsonl
    char dir[256], path[512];
    snprintf(dir, sizeof(dir), "cache/agents/%s", agent_name);
    snprintf(path, sizeof(path), "%s/chat_log.jsonl", dir);

    FILE* f = fopen(path, "a");
    if (!f) return;

    // Escape text for JSON
    char escaped[1024];
    json_escape(text, escaped, sizeof(escaped));

    fprintf(f, "{\"type\":\"%s\",\"agent\":\"%s\",\"text\":\"%s\"}\n",
            type, agent_name, escaped);
    fclose(f);
}

// ---- Config parser for ai_config.yaml ----

void ai_load_config(AiNpc* ai) {
    ai->provider = AI_PROVIDER_LOCAL;
    ai->claude_api_key[0] = '\0';
    snprintf(ai->claude_model, sizeof(ai->claude_model), "claude-sonnet-4-6");

    FILE* f = fopen("config/ai_config.yaml", "r");
    if (!f) {
        printf("[AI] config/ai_config.yaml not found, using local provider.\n");
        return;
    }

    char line[512];
    while (fgets(line, sizeof(line), f)) {
        // Skip comments and empty lines
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;

        char key[64], val[256];
        if (sscanf(p, "%63[^:]: %255[^\n]", key, val) == 2) {
            // Strip quotes from value
            char* v = val;
            while (*v == ' ' || *v == '"') v++;
            int vlen = (int)strlen(v);
            while (vlen > 0 && (v[vlen-1] == '"' || v[vlen-1] == '\n' || v[vlen-1] == '\r')) vlen--;
            v[vlen] = '\0';

            if (strcmp(key, "provider") == 0) {
                if (strcmp(v, "claude") == 0) ai->provider = AI_PROVIDER_CLAUDE;
            } else if (strcmp(key, "api_key") == 0) {
                snprintf(ai->claude_api_key, sizeof(ai->claude_api_key), "%s", v);
            } else if (strcmp(key, "model") == 0 && ai->provider == AI_PROVIDER_CLAUDE) {
                snprintf(ai->claude_model, sizeof(ai->claude_model), "%s", v);
            }
        }
    }
    fclose(f);

    printf("[AI] Provider: %s\n", ai->provider == AI_PROVIDER_CLAUDE ? "Claude API" : "Local");
    if (ai->provider == AI_PROVIDER_CLAUDE) {
        printf("[AI] Model: %s\n", ai->claude_model);
        printf("[AI] API key: %s\n", ai->claude_api_key[0] ? "(set)" : "(empty!)");
    }
    fflush(stdout);
}

// ---- Build Claude API request body ----

static int build_claude_request_body(const AiNpc* ai, char* buf, int buf_size) {
    int pos = 0;

    char escaped_system[2048];
    json_escape(ai->system_prompt, escaped_system, sizeof(escaped_system));

    int max_tok = ai->script_mode ? 4096 : 100;
    pos += snprintf(buf + pos, buf_size - pos,
        "{\"model\":\"%s\",\"max_tokens\":%d,\"temperature\":1.0,"
        "\"system\":\"%s\",\"messages\":[",
        ai->claude_model, max_tok, escaped_system);

    for (int i = 0; i < ai->history_count; i++) {
        bool is_last_user = (strcmp(ai->history[i].role, "user") == 0 &&
                            i == ai->history_count - 1);
        char content_buf[3072];
        if (is_last_user && ai->context_prefix[0]) {
            snprintf(content_buf, sizeof(content_buf), "%s\n%s",
                     ai->context_prefix, ai->history[i].content);
        } else {
            snprintf(content_buf, sizeof(content_buf), "%s", ai->history[i].content);
        }
        char escaped[3072];
        json_escape(content_buf, escaped, sizeof(escaped));
        pos += snprintf(buf + pos, buf_size - pos,
            "{\"role\":\"%s\",\"content\":\"%s\"}%s",
            ai->history[i].role, escaped,
            (i < ai->history_count - 1) ? "," : "");
    }

    pos += snprintf(buf + pos, buf_size - pos, "]}");
    return pos;
}

// ---- Build chat completions request body ----

// Builds the JSON payload for /v1/chat/completions with conversation history.
static int build_request_body(const AiNpc* ai, char* buf, int buf_size) {
    int pos = 0;

    pos += snprintf(buf + pos, buf_size - pos,
        "{\"model\":\"any\",\"messages\":[");

    // System message
    if (ai->system_prompt[0]) {
        char escaped[2048];
        json_escape(ai->system_prompt, escaped, sizeof(escaped));
        pos += snprintf(buf + pos, buf_size - pos,
            "{\"role\":\"system\",\"content\":\"%s\"},", escaped);
    }

    // Conversation history
    for (int i = 0; i < ai->history_count; i++) {
        bool is_last_user = (strcmp(ai->history[i].role, "user") == 0 &&
                            i == ai->history_count - 1);
        // For last user message: prepend context_prefix (mood, surroundings)
        char content_buf[3072];
        if (is_last_user && ai->context_prefix[0]) {
            snprintf(content_buf, sizeof(content_buf), "%s\n%s",
                     ai->context_prefix, ai->history[i].content);
        } else {
            snprintf(content_buf, sizeof(content_buf), "%s", ai->history[i].content);
        }
        char escaped[3072];
        json_escape(content_buf, escaped, sizeof(escaped));
        const char* suffix = is_last_user ? " /no_think" : "";
        pos += snprintf(buf + pos, buf_size - pos,
            "{\"role\":\"%s\",\"content\":\"%s%s\"}%s",
            ai->history[i].role, escaped, suffix,
            (i < ai->history_count - 1) ? "," : "");
    }

    int max_tok = ai->script_mode ? 4096 : 100;
    pos += snprintf(buf + pos, buf_size - pos,
        "],\"temperature\":1.05,\"max_tokens\":%d,\"stream\":false,"
        "\"cache_prompt\":false,\"seed\":-1,"
        "\"repeat_penalty\":1.3,\"repeat_last_n\":128}", max_tok);

    return pos;
}

// ---- Parse the LLM response ----

// Extract the "content" field from the chat completions response,
// then parse the inner JSON (dialogue + actions).
static void parse_response(AiNpc* ai, const char* body) {
    // Clear all output buffers first
    ai->dialogue[0] = '\0';
    ai->script_buf[0] = '\0';
    ai->action_count = 0;

    // Find "content" in the response (nested in choices[0].message.content)
    const char* content_key = "\"content\"";
    const char* p = strstr(body, content_key);
    if (!p) {
        snprintf(ai->error, sizeof(ai->error), "No 'content' in response");
        ai->state = AI_ERROR;
        return;
    }
    p += strlen(content_key);
    while (*p == ' ' || *p == ':' || *p == '\t' || *p == '\n' || *p == '\r') p++;

    // The content value is a JSON string containing our grammar-constrained JSON.
    // It may be a raw string "..." or contain escaped JSON.
    if (*p != '"') {
        snprintf(ai->error, sizeof(ai->error), "Content is not a string");
        ai->state = AI_ERROR;
        return;
    }
    p++; // skip opening quote

    // Unescape the content string into a buffer
    char content[4096];
    int ci = 0;
    while (*p && ci < (int)sizeof(content) - 1) {
        if (*p == '"' && (ci == 0 || *(p-1) != '\\')) break;
        if (*p == '\\' && *(p+1)) {
            p++;
            if (*p == 'n') content[ci++] = '\n';
            else if (*p == 't') content[ci++] = '\t';
            else if (*p == '"') content[ci++] = '"';
            else if (*p == '\\') content[ci++] = '\\';
            else content[ci++] = *p;
        } else {
            content[ci++] = *p;
        }
        p++;
    }
    content[ci] = '\0';

    // Strip Qwen3 <think>...</think> blocks if present
    {
        char* ts = strstr(content, "<think>");
        while (ts) {
            char* te = strstr(ts, "</think>");
            if (te) {
                memmove(ts, te + 8, strlen(te + 8) + 1);
            } else {
                *ts = '\0'; // unclosed think tag, just truncate
            }
            ts = strstr(content, "<think>");
        }
        // Trim leading whitespace after stripping
        char* trimmed = content;
        while (*trimmed == ' ' || *trimmed == '\n' || *trimmed == '\r' || *trimmed == '\t') trimmed++;
        if (trimmed != content) memmove(content, trimmed, strlen(trimmed) + 1);
    }

    // Now parse the inner JSON (our grammar-constrained output)
    printf("[AI] Raw content (%.300s%s)\n", content, strlen(content) > 300 ? "..." : "");
    fflush(stdout);
    // Extract dialogue
    json_find_string(content, "dialogue", ai->dialogue, sizeof(ai->dialogue));

    // Extract actions array
    ai->action_count = 0;
    const char* actions_start = strstr(content, "\"actions\"");
    const char* bracket = actions_start ? strchr(actions_start, '[') : NULL;
    if (!bracket) {
        // No actions array — still a valid response (dialogue only)
        // If no dialogue either, treat the whole content as dialogue
        if (!ai->dialogue[0] && content[0]) {
            snprintf(ai->dialogue, sizeof(ai->dialogue), "%s", content);
        }
        goto finish_parse;
    }

    // Parse each action object
    const char* cursor = bracket + 1;
    while (ai->action_count < AI_MAX_ACTIONS) {
        // Find next action object
        const char* obj_start = strchr(cursor, '{');
        if (!obj_start) break;

        // Find matching closing brace
        const char* obj_end = strchr(obj_start, '}');
        if (!obj_end) break;

        // Extract into a temporary buffer
        int obj_len = (int)(obj_end - obj_start + 1);
        char obj[512];
        if (obj_len >= (int)sizeof(obj)) obj_len = (int)sizeof(obj) - 1;
        memcpy(obj, obj_start, obj_len);
        obj[obj_len] = '\0';

        AiAction* act = &ai->actions[ai->action_count];
        memset(act, 0, sizeof(AiAction));

        char type_str[32];
        if (json_find_string(obj, "type", type_str, sizeof(type_str))) {
            if (strcmp(type_str, "place_block") == 0) {
                act->type = AI_ACTION_PLACE_BLOCK;
                json_find_int(obj, "q", &act->q);
                json_find_int(obj, "r", &act->r);
                json_find_int(obj, "layer", &act->layer);
                json_find_string(obj, "block", act->block, sizeof(act->block));
            } else if (strcmp(type_str, "break_block") == 0) {
                act->type = AI_ACTION_BREAK_BLOCK;
                json_find_int(obj, "q", &act->q);
                json_find_int(obj, "r", &act->r);
                json_find_int(obj, "layer", &act->layer);
            } else if (strcmp(type_str, "say") == 0) {
                act->type = AI_ACTION_SAY;
                json_find_string(obj, "text", act->text, sizeof(act->text));
            } else if (strcmp(type_str, "move_to") == 0) {
                act->type = AI_ACTION_MOVE_TO;
                json_find_int(obj, "q", &act->q);
                json_find_int(obj, "r", &act->r);
            } else if (strcmp(type_str, "look_at") == 0) {
                act->type = AI_ACTION_LOOK_AT;
                json_find_int(obj, "q", &act->q);
                json_find_int(obj, "r", &act->r);
            }
            ai->action_count++;
        }

        cursor = obj_end + 1;
    }

finish_parse:
    // Add assistant response to history
    if (ai->history_count < AI_MAX_HISTORY) {
        AiMessage* msg = &ai->history[ai->history_count++];
        snprintf(msg->role, sizeof(msg->role), "assistant");
        snprintf(msg->content, sizeof(msg->content), "%s", ai->dialogue);
    }

    ai->state = AI_SUCCESS;
    ai_log_message("agent", "agent_msg", ai->dialogue);
}

// ---- Platform-specific implementation ----

#ifdef __EMSCRIPTEN__

// Web build: AI is disabled, stub everything
void ai_npc_init(AiNpc* ai, const char* model_path, const char* grammar_path, int port) {
    (void)model_path; (void)grammar_path; (void)port;
    memset(ai, 0, sizeof(AiNpc));
    printf("[AI] Disabled in web build.\n");
}
void ai_npc_shutdown(AiNpc* ai) { (void)ai; }
void ai_npc_set_system_prompt(AiNpc* ai, const char* prompt) { (void)ai; (void)prompt; }
void ai_npc_send(AiNpc* ai, const char* player_message) {
    (void)player_message;
    ai->state = AI_ERROR;
    snprintf(ai->error, sizeof(ai->error), "AI not available in web build");
}
void ai_npc_poll(AiNpc* ai) { (void)ai; }
bool ai_npc_server_healthy(AiNpc* ai) { (void)ai; return false; }

#elif defined(_WIN32)

#include <windows.h>
#include <winhttp.h>

// ---- Server process management ----

static HANDLE g_server_process = NULL;
static char   g_model_path[512] = "";
static char   g_grammar_path[512] = "";

static void launch_server(AiNpc* ai) {
    if (g_server_process) return; // already running

    // Build command line
    char cmdline[2048];
    snprintf(cmdline, sizeof(cmdline),
        "tools\\ai\\llama-server.exe"
        " -m \"%s\""
        " --port %d"
        " -ngl 20"
        " --grammar-file \"%s\""
        " --ctx-size 4096"
        " --log-disable",
        g_model_path, ai->server_port, g_grammar_path);

    printf("[AI] Launching: %s\n", cmdline);
    fflush(stdout);

    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    memset(&si, 0, sizeof(si));
    si.cb = sizeof(si);
    // Hide the server console window
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = SW_HIDE;
    memset(&pi, 0, sizeof(pi));

    if (CreateProcessA(NULL, cmdline, NULL, NULL, FALSE,
                       CREATE_NO_WINDOW, NULL, NULL, &si, &pi)) {
        g_server_process = pi.hProcess;
        CloseHandle(pi.hThread);
        ai->server_running = true;
        printf("[AI] llama-server started (PID %lu, port %d)\n",
               pi.dwProcessId, ai->server_port);
    } else {
        printf("[AI] Failed to start llama-server (error %lu)\n", GetLastError());
        ai->server_running = false;
    }
    fflush(stdout);
}

static void kill_server(AiNpc* ai) {
    if (g_server_process) {
        TerminateProcess(g_server_process, 0);
        WaitForSingleObject(g_server_process, 3000);
        CloseHandle(g_server_process);
        g_server_process = NULL;
        ai->server_running = false;
        printf("[AI] llama-server stopped.\n");
        fflush(stdout);
    }
}

// ---- HTTP thread (same pattern as lobby.c) ----

typedef struct {
    AiNpc* ai;
    char   post_body[16384]; // large enough for conversation history
    HANDLE thread;
} AiHttpThread;

static AiHttpThread* g_active_ai_http = NULL;

static DWORD WINAPI ai_http_thread_func(LPVOID param) {
    AiHttpThread* data = (AiHttpThread*)param;
    AiNpc* ai = data->ai;

    wchar_t whost[] = L"127.0.0.1";
    INTERNET_PORT port = (INTERNET_PORT)ai->server_port;

    HINTERNET hSession = WinHttpOpen(L"HexPlanets-AI/1.0",
        WINHTTP_ACCESS_TYPE_DEFAULT_PROXY,
        WINHTTP_NO_PROXY_NAME, WINHTTP_NO_PROXY_BYPASS, 0);
    if (!hSession) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpOpen failed");
        ai->state = AI_ERROR;
        return 0;
    }

    HINTERNET hConnect = WinHttpConnect(hSession, whost, port, 0);
    if (!hConnect) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpConnect failed (is llama-server running?)");
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hSession);
        return 0;
    }

    HINTERNET hRequest = WinHttpOpenRequest(hConnect, L"POST",
        L"/v1/chat/completions",
        NULL, WINHTTP_NO_REFERER, WINHTTP_DEFAULT_ACCEPT_TYPES, 0); // no HTTPS for localhost
    if (!hRequest) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpOpenRequest failed");
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    WinHttpAddRequestHeaders(hRequest,
        L"Content-Type: application/json", -1,
        WINHTTP_ADDREQ_FLAG_ADD | WINHTTP_ADDREQ_FLAG_REPLACE);

    // Set longer timeout for LLM inference (120 seconds)
    WinHttpSetTimeouts(hRequest, 60000, 60000, 120000, 120000);

    DWORD body_len = (DWORD)strlen(data->post_body);
    BOOL ok = WinHttpSendRequest(hRequest,
        WINHTTP_NO_ADDITIONAL_HEADERS, 0,
        (LPVOID)data->post_body, body_len, body_len, 0);

    if (!ok) {
        DWORD err = GetLastError();
        snprintf(ai->error, sizeof(ai->error),
            "WinHttpSendRequest failed (%lu) — llama-server may still be loading", err);
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hRequest);
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    ok = WinHttpReceiveResponse(hRequest, NULL);
    if (!ok) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpReceiveResponse failed");
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hRequest);
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    // Get status code
    DWORD status_code = 0;
    DWORD size = sizeof(status_code);
    WinHttpQueryHeaders(hRequest,
        WINHTTP_QUERY_STATUS_CODE | WINHTTP_QUERY_FLAG_NUMBER,
        WINHTTP_HEADER_NAME_BY_INDEX, &status_code, &size, WINHTTP_NO_HEADER_INDEX);

    // Read response body
    char body[16384];
    int body_pos = 0;
    DWORD bytes_read = 0;
    while (body_pos < (int)sizeof(body) - 1) {
        DWORD avail = 0;
        WinHttpQueryDataAvailable(hRequest, &avail);
        if (avail == 0) break;
        if (avail > (DWORD)(sizeof(body) - 1 - body_pos))
            avail = (DWORD)(sizeof(body) - 1 - body_pos);
        WinHttpReadData(hRequest, body + body_pos, avail, &bytes_read);
        body_pos += (int)bytes_read;
    }
    body[body_pos] = '\0';

    WinHttpCloseHandle(hRequest);
    WinHttpCloseHandle(hConnect);
    WinHttpCloseHandle(hSession);

    if (status_code >= 200 && status_code < 300) {
        parse_response(ai, body);
        printf("[AI] Response (HTTP %d, %d actions): \"%s\"\n",
               (int)status_code, ai->action_count, ai->dialogue);
    } else {
        ai->state = AI_ERROR;
        snprintf(ai->error, sizeof(ai->error), "HTTP %d from llama-server", (int)status_code);
        printf("[AI] Error: HTTP %d: %.200s\n", (int)status_code, body);
    }
    fflush(stdout);
    return 0;
}

// ---- Claude API HTTP thread ----

static void parse_claude_response(AiNpc* ai, const char* body) {
    // Clear all output buffers first
    ai->dialogue[0] = '\0';
    ai->script_buf[0] = '\0';
    ai->action_count = 0;

    // Claude returns {"content":[{"type":"text","text":"..."}], ...}
    int blen = (int)strlen(body);
    printf("[AI] Claude raw (first 500): %.500s\n", body);
    if (blen > 500)
        printf("[AI] Claude raw (last 300): %s\n", body + blen - 300);
    fflush(stdout);

    // Find "text" key that's followed by ":" (skip "type":"text" which has "," after)
    const char* p = body;
    while ((p = strstr(p, "\"text\"")) != NULL) {
        p += 6; // skip "text"
        // Skip whitespace
        while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
        if (*p == ':') { p++; break; } // found "text": — this is the content field
    }
    if (!p) {
        snprintf(ai->error, sizeof(ai->error), "No 'text' field in Claude response");
        ai->state = AI_ERROR;
        return;
    }
    while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
    if (*p != '"') {
        snprintf(ai->error, sizeof(ai->error), "Claude 'text' not a string");
        ai->state = AI_ERROR;
        return;
    }
    p++;

    // Unescape full text into a large temp buffer first
    char full_text[AI_MAX_SCRIPT_BUF];
    int ci = 0;
    while (*p && ci < (int)sizeof(full_text) - 1) {
        if (*p == '\\' && *(p+1)) {
            p++;
            if (*p == 'n') full_text[ci++] = '\n';
            else if (*p == '"') full_text[ci++] = '"';
            else if (*p == '\\') full_text[ci++] = '\\';
            else if (*p == 't') full_text[ci++] = '\t';
            else full_text[ci++] = *p;
        } else if (*p == '"') {
            break;
        } else {
            full_text[ci++] = *p;
        }
        p++;
    }
    full_text[ci] = '\0';

    // Log full response
    ai_log_message("agent", "agent_msg", full_text);

    // Extract [SCRIPT]...[/SCRIPT] into script_buf, remove from dialogue
    ai->script_buf[0] = '\0';
    char* script_start = strstr(full_text, "[SCRIPT]");
    char* script_end = script_start ? strstr(script_start, "[/SCRIPT]") : NULL;
    if (script_start) {
        const char* content_start = script_start + 8;
        int slen;
        if (script_end) {
            slen = (int)(script_end - content_start);
        } else {
            // Truncated — take everything after [SCRIPT] and try to close the JSON
            slen = (int)strlen(content_start);
            printf("[AI] Script truncated (max_tokens hit), salvaging %d bytes\n", slen);
            fflush(stdout);
        }
        if (slen > 0 && slen < AI_MAX_SCRIPT_BUF - 4) {
            memcpy(ai->script_buf, content_start, slen);
            ai->script_buf[slen] = '\0';
            // If truncated, try to close the JSON by finding last complete step
            if (!script_end) {
                // Find last complete '}' in a step
                char* last_brace = strrchr(ai->script_buf, '}');
                if (last_brace) {
                    // Close the steps array and object
                    snprintf(last_brace + 1, AI_MAX_SCRIPT_BUF - (int)(last_brace - ai->script_buf) - 1, "]}");
                }
            }
        }
        // Remove script block from text for dialogue
        if (script_end) {
            memmove(script_start, script_end + 9, strlen(script_end + 9) + 1);
        } else {
            *script_start = '\0'; // truncated — just cut it
        }
    }

    // What remains is the dialogue (trim whitespace)
    char* trimmed = full_text;
    while (*trimmed == ' ' || *trimmed == '\n' || *trimmed == '\r') trimmed++;
    int tlen = (int)strlen(trimmed);
    while (tlen > 0 && (trimmed[tlen-1] == ' ' || trimmed[tlen-1] == '\n' || trimmed[tlen-1] == '\r'))
        trimmed[--tlen] = '\0';
    snprintf(ai->dialogue, AI_MAX_DIALOGUE, "%s", trimmed);

    ai->action_count = 0;

    // Add cleaned dialogue to history (not the script)
    if (ai->history_count < AI_MAX_HISTORY) {
        AiMessage* msg = &ai->history[ai->history_count++];
        snprintf(msg->role, sizeof(msg->role), "assistant");
        snprintf(msg->content, sizeof(msg->content), "%s", ai->dialogue);
    }

    ai->state = AI_SUCCESS;
}

static DWORD WINAPI ai_claude_thread_func(LPVOID param) {
    AiHttpThread* data = (AiHttpThread*)param;
    AiNpc* ai = data->ai;

    wchar_t whost[] = L"api.anthropic.com";

    HINTERNET hSession = WinHttpOpen(L"HexPlanets-AI/1.0",
        WINHTTP_ACCESS_TYPE_DEFAULT_PROXY,
        WINHTTP_NO_PROXY_NAME, WINHTTP_NO_PROXY_BYPASS, 0);
    if (!hSession) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpOpen failed");
        ai->state = AI_ERROR;
        return 0;
    }

    HINTERNET hConnect = WinHttpConnect(hSession, whost, INTERNET_DEFAULT_HTTPS_PORT, 0);
    if (!hConnect) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpConnect failed for Claude API");
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hSession);
        return 0;
    }

    HINTERNET hRequest = WinHttpOpenRequest(hConnect, L"POST", L"/v1/messages",
        NULL, WINHTTP_NO_REFERER, WINHTTP_DEFAULT_ACCEPT_TYPES, WINHTTP_FLAG_SECURE);
    if (!hRequest) {
        snprintf(ai->error, sizeof(ai->error), "WinHttpOpenRequest failed");
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    // Add headers
    WinHttpAddRequestHeaders(hRequest,
        L"Content-Type: application/json", -1,
        WINHTTP_ADDREQ_FLAG_ADD | WINHTTP_ADDREQ_FLAG_REPLACE);
    WinHttpAddRequestHeaders(hRequest,
        L"anthropic-version: 2023-06-01", -1,
        WINHTTP_ADDREQ_FLAG_ADD | WINHTTP_ADDREQ_FLAG_REPLACE);

    // Add API key header
    wchar_t api_key_header[256];
    wchar_t wkey[128];
    MultiByteToWideChar(CP_UTF8, 0, ai->claude_api_key, -1, wkey, 128);
    swprintf(api_key_header, 256, L"x-api-key: %s", wkey);
    WinHttpAddRequestHeaders(hRequest, api_key_header, -1,
        WINHTTP_ADDREQ_FLAG_ADD | WINHTTP_ADDREQ_FLAG_REPLACE);

    WinHttpSetTimeouts(hRequest, 60000, 60000, 60000, 60000);

    DWORD body_len = (DWORD)strlen(data->post_body);
    BOOL ok = WinHttpSendRequest(hRequest,
        WINHTTP_NO_ADDITIONAL_HEADERS, 0,
        (LPVOID)data->post_body, body_len, body_len, 0);

    if (!ok) {
        snprintf(ai->error, sizeof(ai->error), "Claude API request failed (%lu)", GetLastError());
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hRequest);
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    ok = WinHttpReceiveResponse(hRequest, NULL);
    if (!ok) {
        snprintf(ai->error, sizeof(ai->error), "Claude API receive failed");
        ai->state = AI_ERROR;
        WinHttpCloseHandle(hRequest);
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    DWORD status_code = 0;
    DWORD sz = sizeof(status_code);
    WinHttpQueryHeaders(hRequest,
        WINHTTP_QUERY_STATUS_CODE | WINHTTP_QUERY_FLAG_NUMBER,
        WINHTTP_HEADER_NAME_BY_INDEX, &status_code, &sz, WINHTTP_NO_HEADER_INDEX);

    char body[16384];
    int body_pos = 0;
    DWORD bytes_read = 0;
    while (body_pos < (int)sizeof(body) - 1) {
        DWORD avail = 0;
        WinHttpQueryDataAvailable(hRequest, &avail);
        if (avail == 0) break;
        if (avail > (DWORD)(sizeof(body) - 1 - body_pos))
            avail = (DWORD)(sizeof(body) - 1 - body_pos);
        WinHttpReadData(hRequest, body + body_pos, avail, &bytes_read);
        body_pos += (int)bytes_read;
    }
    body[body_pos] = '\0';

    WinHttpCloseHandle(hRequest);
    WinHttpCloseHandle(hConnect);
    WinHttpCloseHandle(hSession);

    printf("[AI] Claude HTTP %d (model=%s, body_len=%d)\n",
           (int)status_code, ai->claude_model, body_pos);
    if (status_code >= 200 && status_code < 300) {
        parse_claude_response(ai, body);
        if (ai->state == AI_ERROR) {
            printf("[AI] Claude parse FAILED: %s\n", ai->error);
            printf("[AI] Claude full body: %.*s\n", body_pos < 1000 ? body_pos : 1000, body);
        } else {
            printf("[AI] Claude response: \"%.200s\"\n", ai->dialogue);
        }
    } else {
        ai->state = AI_ERROR;
        // Try to extract error message from response
        const char* emsg = strstr(body, "\"message\"");
        if (emsg) {
            emsg = strchr(emsg + 9, '"');
            if (emsg) {
                emsg++; // skip opening quote
                char errbuf[120];
                int ei = 0;
                while (*emsg && *emsg != '"' && ei < 119) errbuf[ei++] = *emsg++;
                errbuf[ei] = '\0';
                snprintf(ai->error, sizeof(ai->error), "Claude: %s", errbuf);
            } else {
                snprintf(ai->error, sizeof(ai->error), "Claude HTTP %d", (int)status_code);
            }
        } else {
            snprintf(ai->error, sizeof(ai->error), "Claude HTTP %d", (int)status_code);
        }
        printf("[AI] Claude error: HTTP %d\n[AI] Response: %.*s\n",
               (int)status_code, body_pos < 1000 ? body_pos : 1000, body);
    }
    fflush(stdout);
    return 0;
}

// ---- Public API ----

void ai_npc_init(AiNpc* ai, const char* model_path, const char* grammar_path, int port) {
    memset(ai, 0, sizeof(AiNpc));
    ai->server_port = port;

    // Load provider config
    ai_load_config(ai);

    snprintf(g_model_path, sizeof(g_model_path), "%s", model_path);
    snprintf(g_grammar_path, sizeof(g_grammar_path), "%s", grammar_path);

    // Load intrinsic directives from shared config file
    FILE* df = fopen("config/ai_directives.txt", "r");
    if (df) {
        char directives[512];
        int len = (int)fread(directives, 1, sizeof(directives) - 1, df);
        directives[len] = '\0';
        fclose(df);
        snprintf(ai->system_prompt, sizeof(ai->system_prompt),
            "%s\nRespond with dialogue and a list of actions in JSON format.", directives);
        printf("[AI] Loaded directives from config/ai_directives.txt\n");
    } else {
        snprintf(ai->system_prompt, sizeof(ai->system_prompt),
            "You are a worker NPC in a hex-planet voxel game. "
            "You help the player build structures by placing and breaking blocks. "
            "Respond with dialogue and a list of actions. "
            "Be helpful, concise, and stay in character.");
        printf("[AI] config/ai_directives.txt not found, using defaults.\n");
    }
    fflush(stdout);

    // Launch local server if using local provider
    if (ai->provider == AI_PROVIDER_LOCAL) {
        FILE* f = fopen("tools/ai/llama-server.exe", "rb");
        if (f) {
            fclose(f);
            f = fopen(model_path, "rb");
            if (f) {
                fclose(f);
                launch_server(ai);
            } else {
                printf("[AI] Model not found: %s\n", model_path);
                printf("[AI] Run tools\\ai-install.bat to download the AI model.\n");
            }
        } else {
            printf("[AI] llama-server not found. Run tools\\ai-install.bat to install.\n");
        }
    } else {
        ai->server_running = true; // Claude is always "running"
        printf("[AI] Using Claude API — no local server needed.\n");
    }
    fflush(stdout);
}

void ai_npc_shutdown(AiNpc* ai) {
    // Wait for per-instance pending request
    if (ai->_http_data) {
        HANDLE t = (HANDLE)ai->_http_thread;
        if (t) {
            WaitForSingleObject(t, 3000);
            CloseHandle(t);
        }
        free(ai->_http_data);
        ai->_http_data = NULL;
        ai->_http_thread = NULL;
    }
    kill_server(ai);
}

void ai_npc_set_system_prompt(AiNpc* ai, const char* prompt) {
    snprintf(ai->system_prompt, sizeof(ai->system_prompt), "%s", prompt);
}

void ai_npc_send(AiNpc* ai, const char* player_message) {
    if (ai->state == AI_PENDING) {
        printf("[AI] BLOCKED: ai_npc_send called while state=AI_PENDING (still waiting for response)\n");
        fflush(stdout);
        ai->send_blocked = true;
        return;
    }

    if (!ai->server_running) {
        printf("[AI] BLOCKED: server_running=false\n");
        fflush(stdout);
        ai->state = AI_ERROR;
        snprintf(ai->error, sizeof(ai->error), "llama-server not running");
        return;
    }

    // Add player message to history (strip context injection — save only the human text)
    if (ai->history_count >= AI_MAX_HISTORY) {
        memmove(&ai->history[0], &ai->history[2],
                (AI_MAX_HISTORY - 2) * sizeof(AiMessage));
        ai->history_count -= 2;
    }
    AiMessage* msg = &ai->history[ai->history_count++];
    snprintf(msg->role, sizeof(msg->role), "user");
    // Strip [PLAYER SAYS] prefix if present — only save the actual player text
    const char* clean_msg = strstr(player_message, "[PLAYER SAYS]\n");
    if (clean_msg) {
        clean_msg += 14; // skip "[PLAYER SAYS]\n"
    } else {
        clean_msg = player_message;
    }
    snprintf(msg->content, sizeof(msg->content), "%s", clean_msg);

    // Clean up previous per-instance HTTP thread (non-blocking)
    if (ai->_http_data) {
        AiHttpThread* prev = (AiHttpThread*)ai->_http_data;
        HANDLE prev_thread = (HANDLE)ai->_http_thread;
        if (prev_thread) {
            DWORD result = WaitForSingleObject(prev_thread, 0);
            if (result == WAIT_TIMEOUT) {
                printf("[AI] Previous request still pending, please wait.\n");
                fflush(stdout);
                return;
            }
            CloseHandle(prev_thread);
        }
        free(prev);
        ai->_http_data = NULL;
        ai->_http_thread = NULL;
    }

    AiHttpThread* http = (AiHttpThread*)calloc(1, sizeof(AiHttpThread));
    http->ai = ai;
    ai->_http_data = http;

    // Build request body
    if (ai->provider == AI_PROVIDER_CLAUDE) {
        build_claude_request_body(ai, http->post_body, sizeof(http->post_body));
    } else {
        build_request_body(ai, http->post_body, sizeof(http->post_body));
    }
    // Log what we're sending
    {
        // Find max_tokens in the request body
        const char* mt = strstr(http->post_body, "max_tokens");
        char mt_val[16] = "?";
        if (mt) { mt += 12; int mi = 0; while (mt[mi] >= '0' && mt[mi] <= '9' && mi < 15) { mt_val[mi] = mt[mi]; mi++; } mt_val[mi] = '\0'; }
        printf("[AI] Request: script_mode=%d, max_tokens=%s, body_len=%d\n",
               ai->script_mode, mt_val, (int)strlen(http->post_body));
        fflush(stdout);
    }

    ai->state = AI_PENDING;
    ai->dialogue[0] = '\0';
    ai->action_count = 0;

    printf("[AI] Sending via %s: \"%s\"\n",
           ai->provider == AI_PROVIDER_CLAUDE ? "Claude" : "local", player_message);
    fflush(stdout);

    // Log player message
    ai_log_message("agent", "player_msg", player_message);

    HANDLE thread = CreateThread(NULL, 0,
        ai->provider == AI_PROVIDER_CLAUDE ? ai_claude_thread_func : ai_http_thread_func,
        http, 0, NULL);
    ai->_http_thread = (void*)thread;
}

void ai_npc_poll(AiNpc* ai) {
    // State is updated by the background thread — nothing to do here.
    (void)ai;
}

bool ai_npc_server_healthy(AiNpc* ai) {
    if (!ai->server_running) return false;

    // Quick health check: try connecting to the health endpoint
    wchar_t whost[] = L"127.0.0.1";
    INTERNET_PORT port = (INTERNET_PORT)ai->server_port;

    HINTERNET hSession = WinHttpOpen(L"HexPlanets-AI/1.0",
        WINHTTP_ACCESS_TYPE_DEFAULT_PROXY,
        WINHTTP_NO_PROXY_NAME, WINHTTP_NO_PROXY_BYPASS, 0);
    if (!hSession) return false;

    HINTERNET hConnect = WinHttpConnect(hSession, whost, port, 0);
    if (!hConnect) { WinHttpCloseHandle(hSession); return false; }

    HINTERNET hRequest = WinHttpOpenRequest(hConnect, L"GET", L"/health",
        NULL, WINHTTP_NO_REFERER, WINHTTP_DEFAULT_ACCEPT_TYPES, 0);
    if (!hRequest) {
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return false;
    }

    // Short timeout for health check
    WinHttpSetTimeouts(hRequest, 2000, 2000, 2000, 2000);

    BOOL ok = WinHttpSendRequest(hRequest,
        WINHTTP_NO_ADDITIONAL_HEADERS, 0,
        WINHTTP_NO_REQUEST_DATA, 0, 0, 0);
    if (!ok) {
        WinHttpCloseHandle(hRequest);
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return false;
    }

    ok = WinHttpReceiveResponse(hRequest, NULL);
    DWORD status_code = 0;
    DWORD sz = sizeof(status_code);
    if (ok) {
        WinHttpQueryHeaders(hRequest,
            WINHTTP_QUERY_STATUS_CODE | WINHTTP_QUERY_FLAG_NUMBER,
            WINHTTP_HEADER_NAME_BY_INDEX, &status_code, &sz, WINHTTP_NO_HEADER_INDEX);
    }

    WinHttpCloseHandle(hRequest);
    WinHttpCloseHandle(hConnect);
    WinHttpCloseHandle(hSession);

    return (status_code == 200);
}

#else
// Linux/macOS: stub for now
void ai_npc_init(AiNpc* ai, const char* model_path, const char* grammar_path, int port) {
    (void)model_path; (void)grammar_path; (void)port;
    memset(ai, 0, sizeof(AiNpc));
    printf("[AI] Not yet implemented on this platform.\n");
}
void ai_npc_shutdown(AiNpc* ai) { (void)ai; }
void ai_npc_set_system_prompt(AiNpc* ai, const char* prompt) { (void)ai; (void)prompt; }
void ai_npc_send(AiNpc* ai, const char* player_message) {
    (void)player_message;
    ai->state = AI_ERROR;
    snprintf(ai->error, sizeof(ai->error), "AI not implemented on this platform");
}
void ai_npc_poll(AiNpc* ai) { (void)ai; }
bool ai_npc_server_healthy(AiNpc* ai) { (void)ai; return false; }
#endif

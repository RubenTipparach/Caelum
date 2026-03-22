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
        char escaped[1024];
        json_escape(ai->history[i].content, escaped, sizeof(escaped));
        pos += snprintf(buf + pos, buf_size - pos,
            "{\"role\":\"%s\",\"content\":\"%s\"}%s",
            ai->history[i].role, escaped,
            (i < ai->history_count - 1) ? "," : "");
    }

    pos += snprintf(buf + pos, buf_size - pos,
        "],\"temperature\":0.7,\"max_tokens\":1024,\"stream\":false}");

    return pos;
}

// ---- Parse the LLM response ----

// Extract the "content" field from the chat completions response,
// then parse the inner JSON (dialogue + actions).
static void parse_response(AiNpc* ai, const char* body) {
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

    // Now parse the inner JSON (our grammar-constrained output)
    // Extract dialogue
    json_find_string(content, "dialogue", ai->dialogue, sizeof(ai->dialogue));

    // Extract actions array
    ai->action_count = 0;
    const char* actions_start = strstr(content, "\"actions\"");
    if (!actions_start) return;

    // Find the opening bracket
    const char* bracket = strchr(actions_start, '[');
    if (!bracket) return;

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

    // Add assistant response to history
    if (ai->history_count < AI_MAX_HISTORY) {
        AiMessage* msg = &ai->history[ai->history_count++];
        snprintf(msg->role, sizeof(msg->role), "assistant");
        snprintf(msg->content, sizeof(msg->content), "%s", ai->dialogue);
    }

    ai->state = AI_SUCCESS;
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
        " -ngl 99"
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

    // Set longer timeout for LLM inference (60 seconds)
    WinHttpSetTimeouts(hRequest, 60000, 60000, 60000, 60000);

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
    char body[8192];
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
        printf("[AI] Response received (%d actions)\n", ai->action_count);
    } else {
        ai->state = AI_ERROR;
        snprintf(ai->error, sizeof(ai->error), "HTTP %d from llama-server", (int)status_code);
        printf("[AI] Error: HTTP %d: %.200s\n", (int)status_code, body);
    }
    fflush(stdout);
    return 0;
}

// ---- Public API ----

void ai_npc_init(AiNpc* ai, const char* model_path, const char* grammar_path, int port) {
    memset(ai, 0, sizeof(AiNpc));
    ai->server_port = port;

    snprintf(g_model_path, sizeof(g_model_path), "%s", model_path);
    snprintf(g_grammar_path, sizeof(g_grammar_path), "%s", grammar_path);

    // Default system prompt
    snprintf(ai->system_prompt, sizeof(ai->system_prompt),
        "You are a worker NPC in a hex-planet voxel game. "
        "You help the player build structures by placing and breaking blocks. "
        "Respond with dialogue and a list of actions. "
        "Be helpful, concise, and stay in character.");

    // Check if llama-server binary exists
    FILE* f = fopen("tools/ai/llama-server.exe", "rb");
    if (f) {
        fclose(f);
        // Check if model exists
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
    fflush(stdout);
}

void ai_npc_shutdown(AiNpc* ai) {
    // Wait for any pending request
    if (g_active_ai_http) {
        if (g_active_ai_http->thread) {
            WaitForSingleObject(g_active_ai_http->thread, 5000);
            CloseHandle(g_active_ai_http->thread);
        }
        free(g_active_ai_http);
        g_active_ai_http = NULL;
    }
    kill_server(ai);
}

void ai_npc_set_system_prompt(AiNpc* ai, const char* prompt) {
    snprintf(ai->system_prompt, sizeof(ai->system_prompt), "%s", prompt);
}

void ai_npc_send(AiNpc* ai, const char* player_message) {
    if (ai->state == AI_PENDING) return; // already in flight

    if (!ai->server_running) {
        ai->state = AI_ERROR;
        snprintf(ai->error, sizeof(ai->error), "llama-server not running");
        return;
    }

    // Add player message to history
    if (ai->history_count >= AI_MAX_HISTORY) {
        // Shift history: drop oldest 2 messages (1 exchange)
        memmove(&ai->history[0], &ai->history[2],
                (AI_MAX_HISTORY - 2) * sizeof(AiMessage));
        ai->history_count -= 2;
    }
    AiMessage* msg = &ai->history[ai->history_count++];
    snprintf(msg->role, sizeof(msg->role), "user");
    snprintf(msg->content, sizeof(msg->content), "%s", player_message);

    // Clean up previous HTTP thread
    if (g_active_ai_http) {
        if (g_active_ai_http->thread) {
            WaitForSingleObject(g_active_ai_http->thread, 5000);
            CloseHandle(g_active_ai_http->thread);
        }
        free(g_active_ai_http);
        g_active_ai_http = NULL;
    }

    g_active_ai_http = (AiHttpThread*)calloc(1, sizeof(AiHttpThread));
    g_active_ai_http->ai = ai;

    build_request_body(ai, g_active_ai_http->post_body,
                       sizeof(g_active_ai_http->post_body));

    ai->state = AI_PENDING;
    ai->dialogue[0] = '\0';
    ai->action_count = 0;

    printf("[AI] Sending: \"%s\"\n", player_message);
    fflush(stdout);

    g_active_ai_http->thread = CreateThread(NULL, 0,
        ai_http_thread_func, g_active_ai_http, 0, NULL);
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

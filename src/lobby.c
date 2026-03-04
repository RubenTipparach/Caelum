#include "lobby.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char g_api_base[256] = "";

void lobby_init(const char* api_base_url) {
    snprintf(g_api_base, sizeof(g_api_base), "%s", api_base_url);
    printf("[LOBBY] Initialized: %s\n", g_api_base);
    fflush(stdout);
}

void lobby_shutdown(void) {
    g_api_base[0] = '\0';
}

// ---- Minimal JSON field extraction ----
// Extracts the value of "key":"value" from a JSON string.
// Only handles simple string values. Returns false if not found.
static bool json_get_string(const char* json, const char* key, char* out, int out_size) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    // Skip whitespace and colon
    while (*p == ' ' || *p == ':' || *p == '\t') p++;
    if (*p != '"') return false;
    p++; // skip opening quote
    int i = 0;
    while (*p && *p != '"' && i < out_size - 1) {
        out[i++] = *p++;
    }
    out[i] = '\0';
    return true;
}

// Extract integer value: "key":123
static bool json_get_int(const char* json, const char* key, int* out) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "\"%s\"", key);
    const char* p = strstr(json, pattern);
    if (!p) return false;
    p += strlen(pattern);
    while (*p == ' ' || *p == ':' || *p == '\t') p++;
    *out = atoi(p);
    return true;
}

// Parse common fields from lobby response
static void parse_lobby_response(LobbyRequest* req) {
    json_get_string(req->body, "lobby_id", req->lobby_id, sizeof(req->lobby_id));
    json_get_int(req->body, "client_id", &req->client_id);
}

// ---- Platform-specific HTTP ----

#ifdef __EMSCRIPTEN__

#include <emscripten/fetch.h>

typedef struct FetchContext {
    LobbyRequest* req;
} FetchContext;

static void on_fetch_success(emscripten_fetch_t* fetch) {
    FetchContext* ctx = (FetchContext*)fetch->userData;
    LobbyRequest* req = ctx->req;

    req->http_status = fetch->status;
    int len = fetch->numBytes;
    if (len > (int)sizeof(req->body) - 1) len = (int)sizeof(req->body) - 1;
    memcpy(req->body, fetch->data, len);
    req->body[len] = '\0';
    req->body_len = len;

    if (fetch->status >= 200 && fetch->status < 300) {
        req->state = LOBBY_SUCCESS;
        parse_lobby_response(req);
        printf("[LOBBY] Success: %s\n", req->body);
    } else {
        req->state = LOBBY_ERROR;
        snprintf(req->error, sizeof(req->error), "HTTP %d", fetch->status);
        printf("[LOBBY] Error: HTTP %d: %s\n", fetch->status, req->body);
    }
    fflush(stdout);

    free(ctx);
    emscripten_fetch_close(fetch);
}

static void on_fetch_error(emscripten_fetch_t* fetch) {
    FetchContext* ctx = (FetchContext*)fetch->userData;
    LobbyRequest* req = ctx->req;

    req->state = LOBBY_ERROR;
    snprintf(req->error, sizeof(req->error), "Network error (status %d)", fetch->status);
    req->http_status = fetch->status;
    printf("[LOBBY] Fetch error: %s\n", req->error);
    fflush(stdout);

    free(ctx);
    emscripten_fetch_close(fetch);
}

static void do_post(const char* endpoint, const char* json_body, LobbyRequest* req) {
    char url[512];
    snprintf(url, sizeof(url), "%s%s", g_api_base, endpoint);

    FetchContext* ctx = (FetchContext*)malloc(sizeof(FetchContext));
    ctx->req = req;

    emscripten_fetch_attr_t attr;
    emscripten_fetch_attr_init(&attr);
    strcpy(attr.requestMethod, "POST");
    attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
    attr.onsuccess = on_fetch_success;
    attr.onerror = on_fetch_error;
    attr.userData = ctx;

    const char* headers[] = {"Content-Type", "application/json", NULL};
    attr.requestHeaders = headers;
    attr.requestData = json_body;
    attr.requestDataSize = strlen(json_body);

    req->state = LOBBY_PENDING;
    emscripten_fetch(&attr, url);
}

#elif defined(_WIN32)

#include <windows.h>
#include <winhttp.h>

// Background thread for blocking WinHTTP call
typedef struct HttpThreadData {
    LobbyRequest* req;
    char url[512];
    char post_body[1024];
    HANDLE thread;
} HttpThreadData;

// We keep one active thread data at a time (single request in flight)
static HttpThreadData* g_active_http = NULL;

static DWORD WINAPI http_thread_func(LPVOID param) {
    HttpThreadData* data = (HttpThreadData*)param;
    LobbyRequest* req = data->req;

    // Parse URL: expecting https://host/path
    // Find host and path from URL
    const char* url = data->url;
    const char* host_start = url;
    if (strncmp(url, "https://", 8) == 0) host_start = url + 8;
    else if (strncmp(url, "http://", 7) == 0) host_start = url + 7;

    char host[256] = "";
    char path[256] = "/";
    const char* slash = strchr(host_start, '/');
    if (slash) {
        int host_len = (int)(slash - host_start);
        if (host_len > 255) host_len = 255;
        memcpy(host, host_start, host_len);
        host[host_len] = '\0';
        snprintf(path, sizeof(path), "%s", slash);
    } else {
        snprintf(host, sizeof(host), "%s", host_start);
    }

    // Convert host to wide string
    wchar_t whost[256];
    MultiByteToWideChar(CP_UTF8, 0, host, -1, whost, 256);
    wchar_t wpath[256];
    MultiByteToWideChar(CP_UTF8, 0, path, -1, wpath, 256);

    HINTERNET hSession = WinHttpOpen(L"HexPlanets/1.0",
        WINHTTP_ACCESS_TYPE_DEFAULT_PROXY,
        WINHTTP_NO_PROXY_NAME, WINHTTP_NO_PROXY_BYPASS, 0);
    if (!hSession) {
        snprintf(req->error, sizeof(req->error), "WinHttpOpen failed");
        req->state = LOBBY_ERROR;
        return 0;
    }

    HINTERNET hConnect = WinHttpConnect(hSession, whost, INTERNET_DEFAULT_HTTPS_PORT, 0);
    if (!hConnect) {
        snprintf(req->error, sizeof(req->error), "WinHttpConnect failed");
        req->state = LOBBY_ERROR;
        WinHttpCloseHandle(hSession);
        return 0;
    }

    HINTERNET hRequest = WinHttpOpenRequest(hConnect, L"POST", wpath,
        NULL, WINHTTP_NO_REFERER, WINHTTP_DEFAULT_ACCEPT_TYPES, WINHTTP_FLAG_SECURE);
    if (!hRequest) {
        snprintf(req->error, sizeof(req->error), "WinHttpOpenRequest failed");
        req->state = LOBBY_ERROR;
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    // Add Content-Type header
    WinHttpAddRequestHeaders(hRequest,
        L"Content-Type: application/json", -1,
        WINHTTP_ADDREQ_FLAG_ADD | WINHTTP_ADDREQ_FLAG_REPLACE);

    // Send request
    DWORD body_len = (DWORD)strlen(data->post_body);
    BOOL ok = WinHttpSendRequest(hRequest,
        WINHTTP_NO_ADDITIONAL_HEADERS, 0,
        (LPVOID)data->post_body, body_len, body_len, 0);

    if (!ok) {
        snprintf(req->error, sizeof(req->error), "WinHttpSendRequest failed (%lu)", GetLastError());
        req->state = LOBBY_ERROR;
        WinHttpCloseHandle(hRequest);
        WinHttpCloseHandle(hConnect);
        WinHttpCloseHandle(hSession);
        return 0;
    }

    ok = WinHttpReceiveResponse(hRequest, NULL);
    if (!ok) {
        snprintf(req->error, sizeof(req->error), "WinHttpReceiveResponse failed");
        req->state = LOBBY_ERROR;
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
    req->http_status = (int)status_code;

    // Read response body
    req->body_len = 0;
    DWORD bytes_read = 0;
    while (req->body_len < (int)sizeof(req->body) - 1) {
        DWORD avail = 0;
        WinHttpQueryDataAvailable(hRequest, &avail);
        if (avail == 0) break;
        if (avail > (DWORD)(sizeof(req->body) - 1 - req->body_len))
            avail = (DWORD)(sizeof(req->body) - 1 - req->body_len);
        WinHttpReadData(hRequest, req->body + req->body_len, avail, &bytes_read);
        req->body_len += (int)bytes_read;
    }
    req->body[req->body_len] = '\0';

    WinHttpCloseHandle(hRequest);
    WinHttpCloseHandle(hConnect);
    WinHttpCloseHandle(hSession);

    if (status_code >= 200 && status_code < 300) {
        req->state = LOBBY_SUCCESS;
        parse_lobby_response(req);
        printf("[LOBBY] Success: %s\n", req->body);
    } else {
        req->state = LOBBY_ERROR;
        snprintf(req->error, sizeof(req->error), "HTTP %d", (int)status_code);
        printf("[LOBBY] Error: HTTP %d: %s\n", (int)status_code, req->body);
    }
    fflush(stdout);
    return 0;
}

static void do_post(const char* endpoint, const char* json_body, LobbyRequest* req) {
    // Clean up previous request if any
    if (g_active_http) {
        if (g_active_http->thread) {
            WaitForSingleObject(g_active_http->thread, 5000);
            CloseHandle(g_active_http->thread);
        }
        free(g_active_http);
        g_active_http = NULL;
    }

    g_active_http = (HttpThreadData*)calloc(1, sizeof(HttpThreadData));
    g_active_http->req = req;
    snprintf(g_active_http->url, sizeof(g_active_http->url), "%s%s", g_api_base, endpoint);
    snprintf(g_active_http->post_body, sizeof(g_active_http->post_body), "%s", json_body);

    req->state = LOBBY_PENDING;
    g_active_http->thread = CreateThread(NULL, 0, http_thread_func, g_active_http, 0, NULL);
}

#else
// Linux/macOS fallback: stub (not implemented yet)
static void do_post(const char* endpoint, const char* json_body, LobbyRequest* req) {
    (void)endpoint; (void)json_body;
    req->state = LOBBY_ERROR;
    snprintf(req->error, sizeof(req->error), "HTTP not implemented on this platform");
}
#endif

// ---- Public API ----

void lobby_create(const char* host_name, LobbyRequest* req) {
    memset(req, 0, sizeof(LobbyRequest));
    char body[512];
    snprintf(body, sizeof(body), "{\"host_name\":\"%s\"}", host_name);
    printf("[LOBBY] Creating lobby (host=%s)...\n", host_name);
    fflush(stdout);
    do_post("/api/lobby/create", body, req);
}

void lobby_join(const char* lobby_id, const char* player_name, LobbyRequest* req) {
    memset(req, 0, sizeof(LobbyRequest));
    char body[512];
    snprintf(body, sizeof(body), "{\"lobby_id\":\"%s\",\"player_name\":\"%s\"}", lobby_id, player_name);
    printf("[LOBBY] Joining lobby %s (player=%s)...\n", lobby_id, player_name);
    fflush(stdout);
    do_post("/api/lobby/join", body, req);
}

void lobby_poll(LobbyRequest* req) {
    // For WinHTTP: the background thread sets req->state directly.
    // For Emscripten: the fetch callback sets req->state directly.
    // Nothing to do here — state is updated asynchronously.
    (void)req;
}

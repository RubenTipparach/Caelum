#ifndef AI_NPC_H
#define AI_NPC_H

#include <stdbool.h>
#include <stdint.h>

// --- AI NPC System ---
// Supports local llama-server or Claude API (configured via config/ai_config.yaml).
// Non-blocking: call ai_npc_poll() each frame while a request is pending.

typedef enum {
    AI_PROVIDER_LOCAL,
    AI_PROVIDER_CLAUDE,
} AiProvider;

typedef enum {
    AI_IDLE,        // No request in flight
    AI_PENDING,     // Waiting for LLM response
    AI_SUCCESS,     // Response received and parsed
    AI_ERROR,       // Request failed
} AiState;

// Single action from the NPC
typedef enum {
    AI_ACTION_PLACE_BLOCK,
    AI_ACTION_BREAK_BLOCK,
    AI_ACTION_SAY,
    AI_ACTION_MOVE_TO,
    AI_ACTION_LOOK_AT,
} AiActionType;

typedef struct {
    AiActionType type;
    int q, r, layer;        // hex coords (place/break)
    char block[16];         // block type name (place)
    char text[256];         // say text
} AiAction;

#define AI_MAX_ACTIONS 32
#define AI_MAX_DIALOGUE 512
#define AI_MAX_SCRIPT_BUF 8192
#define AI_MAX_HISTORY 32

typedef struct {
    char role[16];          // "user" or "assistant"
    char content[512];
} AiMessage;

typedef struct {
    AiState state;
    int     http_status;
    char    error[128];

    // Parsed response
    char      dialogue[AI_MAX_DIALOGUE];
    char      script_buf[AI_MAX_SCRIPT_BUF]; // extracted [SCRIPT]...[/SCRIPT] content
    AiAction  actions[AI_MAX_ACTIONS];
    int       action_count;

    // Conversation history (rolling window)
    AiMessage history[AI_MAX_HISTORY];
    int       history_count;

    // System prompt (personality/directives)
    char      system_prompt[1024];

    // Provider config
    AiProvider provider;
    bool      server_running;
    int       server_port;
    char      claude_api_key[128];
    char      claude_model[64];

    // Transient context — prepended to the last user message in the request body only
    // Not saved to history. Set before ai_npc_send(), cleared after.
    char      context_prefix[2048];

    // Script generation mode — when true, response is parsed as a JSON script
    bool      script_mode;

    // Vision: base64 JPEG to send with next message (set before ai_npc_send)
    const char* vision_base64;   // pointer to base64 data (not owned)
    int         vision_base64_len;

    // Per-instance async HTTP (not shared global)
    void*     _http_thread;   // HANDLE, managed internally
    void*     _http_data;     // AiHttpThread*, managed internally

    bool      send_blocked;   // set when user tried to send while AI_PENDING
} AiNpc;

// Load provider config from config/ai_config.yaml into an AiNpc struct.
void ai_load_config(AiNpc* ai);

// Initialize AI system. Call once at startup.
// model_path: path to GGUF model file (e.g. "tools/ai/Qwen3-4B-Q4_K_M.gguf")
// grammar_path: path to GBNF grammar file (e.g. "tools/ai/grammar.gbnf")
// port: localhost port for llama-server (default 8080)
void ai_npc_init(AiNpc* ai, const char* model_path, const char* grammar_path, int port);

// Shut down: kill llama-server process, clean up.
void ai_npc_shutdown(AiNpc* ai);

// Set the system prompt (personality, directives, world context).
void ai_npc_set_system_prompt(AiNpc* ai, const char* prompt);

// Send a player message to the NPC (non-blocking).
// Poll with ai_npc_poll() until state != AI_PENDING.
void ai_npc_send(AiNpc* ai, const char* player_message);

// Poll for async HTTP completion. Call each frame.
void ai_npc_poll(AiNpc* ai);

// Check if llama-server is available (health endpoint).
bool ai_npc_server_healthy(AiNpc* ai);

#endif

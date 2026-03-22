#ifndef AI_MEMORY_H
#define AI_MEMORY_H

#include <stdbool.h>

// --- AI Persistent Memory ---
// Each agent has a memory directory: cache/agents/<name>/
//
// Files:
//   memory.md      — Long-term knowledge. Human-editable markdown.
//                     Read into LLM system prompt on every message.
//                     The LLM can request memory writes via [REMEMBER] tags.
//
//   journal.jsonl   — Action outcomes, discoveries, events.
//                     Appended automatically after scripts/actions complete.
//                     Recent entries loaded into context for continuity.
//
//   chat_log.jsonl  — Raw conversation log (already exists in ai_npc.c).
//                     Last N entries loaded on startup for conversation continuity.
//
// All files are plain text and can be reviewed/edited outside the game.

#define AI_MEMORY_MAX       4096    // max memory.md content loaded
#define AI_JOURNAL_MAX      2048    // max recent journal entries loaded
#define AI_CHAT_HISTORY_MAX 10      // number of recent chat entries to load

typedef struct {
    char agent_dir[256];            // e.g. "cache/agents/waltery"
    char memory[AI_MEMORY_MAX];     // contents of memory.md
    char journal[AI_JOURNAL_MAX];   // recent journal entries
    bool loaded;
} AiMemory;

// Initialize memory for an agent. agent_dir is the base directory.
void ai_memory_init(AiMemory* mem, const char* agent_dir);

// Load memory.md and recent journal entries from disk.
void ai_memory_load(AiMemory* mem);

// Append a line to memory.md (persistent knowledge).
void ai_memory_remember(AiMemory* mem, const char* text);

// Append a journal entry (timestamped event/outcome).
void ai_memory_journal(AiMemory* mem, const char* event);

// Load recent chat history from chat_log.jsonl into an AiNpc's history buffer.
// Populates npc->history[] and npc->history_count.
#include "ai_npc.h"
void ai_memory_load_chat_history(AiMemory* mem, AiNpc* npc);

// Build the full context string (memory + recent journal) for injection
// into the system prompt. Returns chars written.
int ai_memory_build_context(const AiMemory* mem, char* buf, int buf_size);

// Scan an LLM response for [REMEMBER] tags and save them to memory.md.
// Format: [REMEMBER]something to remember[/REMEMBER]
// Returns number of memories extracted.
int ai_memory_extract_from_response(AiMemory* mem, const char* response);

#endif

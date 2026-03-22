#include "ai_memory.h"
#include "ai_npc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#define mkdir_p(path) _mkdir(path)
#else
#include <sys/stat.h>
#define mkdir_p(path) mkdir(path, 0755)
#endif

void ai_memory_init(AiMemory* mem, const char* agent_dir) {
    memset(mem, 0, sizeof(AiMemory));
    snprintf(mem->agent_dir, sizeof(mem->agent_dir), "%s", agent_dir);

    // Ensure directory exists
    mkdir_p(agent_dir);

    // Create memory.md if it doesn't exist
    char path[512];
    snprintf(path, sizeof(path), "%s/memory.md", agent_dir);
    FILE* f = fopen(path, "r");
    if (!f) {
        f = fopen(path, "w");
        if (f) {
            fprintf(f, "# Walter's Memory\n\n"
                       "<!-- This file is Walter's persistent knowledge. -->\n"
                       "<!-- Edit freely — Walter reads this on every conversation. -->\n"
                       "<!-- The game appends new memories below. -->\n\n"
                       "## Things I Know\n\n"
                       "## Things the Player Told Me\n\n"
                       "## Places I've Been\n\n"
                       "## Tasks Completed\n\n");
            fclose(f);
            printf("[MEMORY] Created %s\n", path);
        }
    } else {
        fclose(f);
    }

    fflush(stdout);
}

void ai_memory_load(AiMemory* mem) {
    // Load memory.md
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/memory.md", mem->agent_dir);
        FILE* f = fopen(path, "rb");
        if (f) {
            int len = (int)fread(mem->memory, 1, AI_MEMORY_MAX - 1, f);
            mem->memory[len] = '\0';
            fclose(f);
            printf("[MEMORY] Loaded memory.md (%d bytes)\n", len);
        } else {
            mem->memory[0] = '\0';
        }
    }

    // Load recent journal entries (last N lines)
    {
        char path[512];
        snprintf(path, sizeof(path), "%s/journal.jsonl", mem->agent_dir);
        FILE* f = fopen(path, "rb");
        if (f) {
            // Read entire file, then take last ~2KB
            fseek(f, 0, SEEK_END);
            long size = ftell(f);
            if (size > AI_JOURNAL_MAX - 1) {
                fseek(f, size - (AI_JOURNAL_MAX - 1), SEEK_SET);
                // Skip to next newline to avoid partial line
                int c;
                while ((c = fgetc(f)) != EOF && c != '\n') {}
            } else {
                fseek(f, 0, SEEK_SET);
            }
            int len = (int)fread(mem->journal, 1, AI_JOURNAL_MAX - 1, f);
            mem->journal[len] = '\0';
            fclose(f);
            printf("[MEMORY] Loaded journal (%d bytes)\n", len);
        } else {
            mem->journal[0] = '\0';
        }
    }

    mem->loaded = true;
    fflush(stdout);
}

void ai_memory_remember(AiMemory* mem, const char* text) {
    char path[512];
    snprintf(path, sizeof(path), "%s/memory.md", mem->agent_dir);
    FILE* f = fopen(path, "a");
    if (f) {
        // Get timestamp
        time_t now = time(NULL);
        struct tm* t = localtime(&now);
        char ts[32];
        strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M", t);

        fprintf(f, "- [%s] %s\n", ts, text);
        fclose(f);

        // Reload into buffer
        ai_memory_load(mem);

        printf("[MEMORY] Remembered: %s\n", text);
        fflush(stdout);
    }
}

void ai_memory_journal(AiMemory* mem, const char* event) {
    char path[512];
    snprintf(path, sizeof(path), "%s/journal.jsonl", mem->agent_dir);
    FILE* f = fopen(path, "a");
    if (f) {
        time_t now = time(NULL);
        struct tm* t = localtime(&now);
        char ts[32];
        strftime(ts, sizeof(ts), "%Y-%m-%dT%H:%M:%S", t);

        // Simple JSON line
        fprintf(f, "{\"ts\":\"%s\",\"event\":\"%s\"}\n", ts, event);
        fclose(f);
    }
}

void ai_memory_load_chat_history(AiMemory* mem, AiNpc* npc) {
    char path[512];
    snprintf(path, sizeof(path), "%s/chat_log.jsonl", mem->agent_dir);
    FILE* f = fopen(path, "rb");
    if (!f) return;

    // Read all lines, keep last AI_CHAT_HISTORY_MAX pairs
    char lines[AI_CHAT_HISTORY_MAX * 2][512];
    int line_count = 0;
    char line_buf[512];

    while (fgets(line_buf, sizeof(line_buf), f)) {
        // Circular buffer of recent lines
        snprintf(lines[line_count % (AI_CHAT_HISTORY_MAX * 2)], 512, "%s", line_buf);
        line_count++;
    }
    fclose(f);

    // Extract messages from the most recent lines
    int start = (line_count > AI_CHAT_HISTORY_MAX * 2) ?
                line_count - AI_CHAT_HISTORY_MAX * 2 : 0;
    int loaded = 0;

    npc->history_count = 0;
    for (int i = start; i < line_count && npc->history_count < AI_MAX_HISTORY; i++) {
        const char* line = lines[i % (AI_CHAT_HISTORY_MAX * 2)];

        // Parse JSONL: {"ts":"...","type":"player_msg"|"agent_msg","agent":"...","text":"..."}
        const char* type_ptr = strstr(line, "\"type\":\"");
        const char* text_ptr = strstr(line, "\"text\":\"");
        if (!type_ptr || !text_ptr) continue;

        type_ptr += 8;
        text_ptr += 8;

        // Determine role
        const char* role = NULL;
        if (strncmp(type_ptr, "player_msg", 10) == 0) role = "user";
        else if (strncmp(type_ptr, "agent_msg", 9) == 0) role = "assistant";
        if (!role) continue;

        // Extract text value
        char text[512];
        int ti = 0;
        const char* p = text_ptr;
        while (*p && *p != '"' && ti < 510) {
            if (*p == '\\' && *(p+1) == 'n') { text[ti++] = ' '; p += 2; continue; }
            if (*p == '\\' && *(p+1) == '"') { text[ti++] = '"'; p += 2; continue; }
            if (*p == '\\' && *(p+1)) { p++; }
            text[ti++] = *p++;
        }
        text[ti] = '\0';

        AiMessage* msg = &npc->history[npc->history_count++];
        snprintf(msg->role, sizeof(msg->role), "%s", role);
        snprintf(msg->content, sizeof(msg->content), "%s", text);
        loaded++;
    }

    if (loaded > 0) {
        printf("[MEMORY] Loaded %d chat history entries\n", loaded);
        fflush(stdout);
    }
}

int ai_memory_build_context(const AiMemory* mem, char* buf, int buf_size) {
    int pos = 0;

    if (mem->memory[0]) {
        pos += snprintf(buf + pos, buf_size - pos,
            "\n=== YOUR MEMORY (persistent knowledge) ===\n%s\n",
            mem->memory);
    }

    if (mem->journal[0]) {
        pos += snprintf(buf + pos, buf_size - pos,
            "\n=== RECENT EVENTS ===\n%s\n",
            mem->journal);
    }

    pos += snprintf(buf + pos, buf_size - pos,
        "\n=== MEMORY INSTRUCTIONS ===\n"
        "To remember something important, include [REMEMBER]text[/REMEMBER] in your response.\n"
        "This saves it to your persistent memory file for future conversations.\n"
        "Only remember things that are important or that the player asks you to remember.\n");

    return pos;
}

int ai_memory_extract_from_response(AiMemory* mem, const char* response) {
    int count = 0;
    const char* p = response;

    while ((p = strstr(p, "[REMEMBER]")) != NULL) {
        p += 10; // skip [REMEMBER]
        const char* end = strstr(p, "[/REMEMBER]");
        if (!end) break;

        int len = (int)(end - p);
        if (len > 0 && len < 500) {
            char text[512];
            memcpy(text, p, len);
            text[len] = '\0';

            // Trim whitespace
            char* t = text;
            while (*t == ' ' || *t == '\n') t++;
            int tlen = (int)strlen(t);
            while (tlen > 0 && (t[tlen-1] == ' ' || t[tlen-1] == '\n')) tlen--;
            t[tlen] = '\0';

            if (tlen > 0) {
                ai_memory_remember(mem, t);
                count++;
            }
        }

        p = end + 11; // skip [/REMEMBER]
    }

    return count;
}

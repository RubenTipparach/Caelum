#include "ai_emotions.h"
#include "ai_agent.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

void ai_emotions_init(AiEmotions* emo) {
    memset(emo, 0, sizeof(AiEmotions));
    // Start neutral
    emo->joy = 0.3f;
    emo->boredom = 0.0f;
    emo->annoyance = 0.0f;
    emo->curiosity = 0.4f;
    emo->satisfaction = 0.2f;

    // Personality-driven rates
    emo->joy_decay = 0.005f;          // joy fades slowly
    emo->boredom_rate = 0.02f;        // gets bored after ~50s idle
    emo->annoyance_decay = 0.01f;     // cools down over ~100s
    emo->curiosity_decay = 0.008f;    // curiosity fades moderately
    emo->satisfaction_decay = 0.003f; // satisfaction lingers
}

static float clampf(float v, float lo, float hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

void ai_emotions_update(AiEmotions* emo, float dt, int agent_state) {
    // Decay emotions toward neutral
    emo->joy = clampf(emo->joy - emo->joy_decay * dt, 0.0f, 1.0f);
    emo->annoyance = clampf(emo->annoyance - emo->annoyance_decay * dt, 0.0f, 1.0f);
    emo->curiosity = clampf(emo->curiosity - emo->curiosity_decay * dt, 0.0f, 1.0f);
    emo->satisfaction = clampf(emo->satisfaction - emo->satisfaction_decay * dt, 0.0f, 1.0f);

    // Boredom: grows when idle, shrinks when active
    if (agent_state == AGENT_STATE_IDLE || agent_state == AGENT_STATE_SLEEPING) {
        emo->idle_time += dt;
        emo->boredom = clampf(emo->boredom + emo->boredom_rate * dt, 0.0f, 1.0f);
    } else {
        emo->idle_time = 0.0f;
        // Active states reduce boredom
        emo->boredom = clampf(emo->boredom - 0.05f * dt, 0.0f, 1.0f);
    }

    // Building states boost joy
    if (agent_state == AGENT_STATE_WORKING || agent_state == AGENT_STATE_PLACING) {
        emo->joy = clampf(emo->joy + 0.03f * dt, 0.0f, 1.0f);
    }

    // Walking/exploring boosts curiosity slightly
    if (agent_state == AGENT_STATE_WALKING) {
        emo->curiosity = clampf(emo->curiosity + 0.01f * dt, 0.0f, 1.0f);
    }

    // Message spam tracking — window resets every 30 seconds
    emo->msg_window_timer += dt;
    if (emo->msg_window_timer > 30.0f) {
        emo->msg_window_timer = 0.0f;
        emo->msgs_in_window = 0;
    }

    emo->msg_cooldown += dt;
}

void ai_emotions_on_message_received(AiEmotions* emo) {
    emo->msgs_in_window++;

    // Rapid-fire messages increase annoyance
    if (emo->msg_cooldown < 3.0f) {
        // Less than 3 seconds since last message — spamming
        emo->annoyance = clampf(emo->annoyance + 0.1f, 0.0f, 1.0f);
    }

    // Many messages in window — getting pestered
    if (emo->msgs_in_window > 5) {
        emo->annoyance = clampf(emo->annoyance + 0.05f * (emo->msgs_in_window - 5), 0.0f, 1.0f);
    }

    // Conversation reduces boredom
    emo->boredom = clampf(emo->boredom - 0.15f, 0.0f, 1.0f);

    // Reset cooldown
    emo->msg_cooldown = 0.0f;
}

void ai_emotions_on_task_complete(AiEmotions* emo) {
    emo->satisfaction = clampf(emo->satisfaction + 0.3f, 0.0f, 1.0f);
    emo->joy = clampf(emo->joy + 0.2f, 0.0f, 1.0f);
    emo->boredom = clampf(emo->boredom - 0.3f, 0.0f, 1.0f);
}

void ai_emotions_on_task_failed(AiEmotions* emo) {
    emo->annoyance = clampf(emo->annoyance + 0.15f, 0.0f, 1.0f);
    emo->satisfaction = clampf(emo->satisfaction - 0.1f, 0.0f, 1.0f);
}

void ai_emotions_on_building(AiEmotions* emo) {
    // Walter LOVES building — big joy boost
    emo->joy = clampf(emo->joy + 0.1f, 0.0f, 1.0f);
    emo->boredom = clampf(emo->boredom - 0.2f, 0.0f, 1.0f);
    emo->satisfaction = clampf(emo->satisfaction + 0.05f, 0.0f, 1.0f);
}

void ai_emotions_on_exploring(AiEmotions* emo) {
    emo->curiosity = clampf(emo->curiosity + 0.15f, 0.0f, 1.0f);
    emo->boredom = clampf(emo->boredom - 0.1f, 0.0f, 1.0f);
}

void ai_emotions_on_new_order(AiEmotions* emo) {
    emo->curiosity = clampf(emo->curiosity + 0.2f, 0.0f, 1.0f);
    emo->boredom = clampf(emo->boredom - 0.3f, 0.0f, 1.0f);
    emo->joy = clampf(emo->joy + 0.1f, 0.0f, 1.0f);
    // Getting a new task reduces annoyance (purpose!)
    emo->annoyance = clampf(emo->annoyance - 0.1f, 0.0f, 1.0f);
}

int ai_emotions_describe(const AiEmotions* emo, char* buf, int buf_size) {
    int pos = 0;
    pos += snprintf(buf + pos, buf_size - pos, "=== YOUR CURRENT MOOD ===\n");

    // Primary emotion
    if (emo->annoyance > 0.7f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're irritated. The player keeps bugging you with repeated stuff. "
            "Be short, sarcastic, or tell them to knock it off.\n");
    } else if (emo->annoyance > 0.4f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're getting a bit annoyed. Keep it together but you can show some edge.\n");
    } else if (emo->boredom > 0.7f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're really bored. You want something to do — building, exploring, anything. "
            "Complain about having nothing to do or suggest a project.\n");
    } else if (emo->boredom > 0.4f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're getting a little restless. Wouldn't mind having something to work on.\n");
    } else if (emo->joy > 0.7f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're in a great mood! Building is going well, life is good. Be upbeat and enthusiastic.\n");
    } else if (emo->satisfaction > 0.6f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You feel accomplished. Recent work went well. Be confident and content.\n");
    } else if (emo->curiosity > 0.6f) {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're curious about something — the terrain, a new area, or what the player is planning.\n");
    } else {
        pos += snprintf(buf + pos, buf_size - pos,
            "You're in a neutral mood. Calm, ready for whatever.\n");
    }

    return pos;
}

int ai_emotions_wants_to_build(const AiEmotions* emo) {
    // Bored enough + has enough satisfaction from past building = wants to build
    return emo->boredom > 0.6f && emo->idle_time > 30.0f;
}

int ai_emotions_is_annoyed(const AiEmotions* emo) {
    return emo->annoyance > 0.5f;
}

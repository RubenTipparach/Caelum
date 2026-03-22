#ifndef AI_EMOTIONS_H
#define AI_EMOTIONS_H

// --- AI Emotion System ---
// Simple emotional state that drifts over time based on activities.
// Injected into LLM context so responses reflect mood.
// Drives idle behavior (bored → build, happy → whistle/celebrate).
//
// All values are 0.0 to 1.0.

typedef struct {
    float joy;          // 0=neutral, 1=ecstatic. Building, completing tasks, good conversation
    float boredom;      // 0=engaged, 1=bored stiff. Increases while idle, decreases with activity
    float annoyance;    // 0=chill, 1=fed up. Rapid-fire messages, interruptions, repeated failures
    float curiosity;    // 0=uninterested, 1=fascinated. New terrain, new orders, exploration
    float satisfaction; // 0=unfulfilled, 1=deeply satisfied. Completing tasks, building things

    // Decay rates (per second)
    float joy_decay;
    float boredom_rate;     // how fast boredom grows when idle
    float annoyance_decay;
    float curiosity_decay;
    float satisfaction_decay;

    // Tracking
    float idle_time;         // seconds spent idle (for boredom)
    float msg_cooldown;      // seconds since last player message (for annoyance)
    int   msgs_in_window;    // messages received in last 30 seconds
    float msg_window_timer;  // 30-second window tracker
} AiEmotions;

// Initialize with default personality values
void ai_emotions_init(AiEmotions* emo);

// Update emotions each frame (decay, boredom growth, etc.)
void ai_emotions_update(AiEmotions* emo, float dt, int agent_state);

// Event triggers — call these when things happen
void ai_emotions_on_message_received(AiEmotions* emo);  // player talked to us
void ai_emotions_on_task_complete(AiEmotions* emo);      // finished a script/order
void ai_emotions_on_task_failed(AiEmotions* emo);        // action failed
void ai_emotions_on_building(AiEmotions* emo);           // placed a block
void ai_emotions_on_exploring(AiEmotions* emo);          // walking to new area
void ai_emotions_on_new_order(AiEmotions* emo);          // received a new order

// Build a mood description string for LLM context injection.
// Returns chars written.
int ai_emotions_describe(const AiEmotions* emo, char* buf, int buf_size);

// Should the agent try to self-initiate building? (boredom threshold)
// Returns true when bored enough to start a creative project
int ai_emotions_wants_to_build(const AiEmotions* emo);

// Is the agent too annoyed to chat nicely?
int ai_emotions_is_annoyed(const AiEmotions* emo);

#endif

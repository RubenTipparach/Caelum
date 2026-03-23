#ifndef AI_VISION_H
#define AI_VISION_H

#include <stdbool.h>

// --- AI Vision System ---
// Captures what the agent "sees" (the current game frame) and encodes it
// for sending to Claude's vision API. Images are saved to the agent's
// cache folder for review.

#define AI_VISION_MAX_BASE64 (256 * 1024)  // 256KB base64 budget

typedef struct {
    char agent_dir[256];
    char last_image_path[512];
    char base64_buf[AI_VISION_MAX_BASE64];
    int base64_len;
    bool capture_requested;
    bool capture_ready;
    int capture_id;
} AiVision;

// Initialize vision for an agent
void ai_vision_init(AiVision* vis, const char* agent_dir);

// Request a capture on the next frame
void ai_vision_request_capture(AiVision* vis);

// Capture the current backbuffer, downscale, save as JPEG, base64 encode.
// Call AFTER render_frame / sg_commit (when backbuffer has the frame).
// Returns true if capture succeeded.
bool ai_vision_capture(AiVision* vis);

// Is a capture ready to send?
bool ai_vision_ready(const AiVision* vis);

// Clear the capture after sending
void ai_vision_clear(AiVision* vis);

#endif

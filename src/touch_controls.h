#ifndef TOUCH_CONTROLS_H
#define TOUCH_CONTROLS_H

#include <stdbool.h>
#include <stdint.h>
#include "sokol_app.h"

/* Forward declaration — camera.h owns the full definition */
typedef struct Camera Camera;

typedef struct TouchJoystick {
    bool active;
    uintptr_t touch_id;
    float origin_x, origin_y;      /* Where the touch started (center of joystick) */
    float current_x, current_y;    /* Current finger position */
    float dx, dy;                  /* Normalized output (-1..1) */
    float radius;                  /* Joystick visual radius (pixels) */
} TouchJoystick;

typedef struct TouchButton {
    float x, y, w, h;      /* Screen rect (computed from layout) */
    bool pressed;
    bool just_released;     /* Set true on finger-up, cleared by caller */
    uintptr_t touch_id;
} TouchButton;

typedef struct TouchControls {
    bool is_touch_device;       /* Auto-detected on first TOUCHES_BEGAN */

    TouchJoystick move_stick;   /* Left side — floating joystick */
    TouchJoystick look_stick;   /* Right side — trackpad (delta-based) */

    TouchButton btn_jump;       /* Green */
    TouchButton btn_crouch;     /* Blue */
    TouchButton btn_break;      /* Red */
    TouchButton btn_place;      /* Yellow */

    int hotbar_touch_slot;      /* -1 if no hotbar tap this event */
} TouchControls;

void touch_init(TouchControls* tc);
bool touch_handle_event(TouchControls* tc, const sapp_event* ev, int hotbar_count);
void touch_apply_to_camera(TouchControls* tc, Camera* cam, float dt);
void touch_render(const TouchControls* tc, int hotbar_selected, int hotbar_count);
void touch_render_crosshair(void);

#endif

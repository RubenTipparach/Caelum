// platform_macos.m — macOS-specific platform fixes
#ifdef __APPLE__

#import <AppKit/AppKit.h>
#import <CoreGraphics/CoreGraphics.h>
#import <objc/runtime.h>
#include <stdio.h>
#include <stdatomic.h>

// ---- Fix: windowDidChangeOcclusionState crash ----
// When the window becomes occluded (alt-tab, another window covers it),
// AppKit's NSWindow.occlusionState property getter uses an internal
// os_unfair_lock. If this notification fires during a dispatch queue drain
// (e.g., from a display-link callback context), the lock can be in a bad
// state, causing "Unlock of an os_unfair_lock not owned by current thread".
//
// Fix: swizzle sokol's windowDidChangeOcclusionState: to defer execution
// via dispatch_async, so it runs outside any re-entrant context.

static IMP g_original_occlusion_imp;

static void swizzled_windowDidChangeOcclusionState(id self, SEL _cmd, NSNotification* notification) {
    // Defer to next main queue iteration to avoid re-entrant lock issues
    NSNotification* notifCopy = notification;  // block captures
    dispatch_async(dispatch_get_main_queue(), ^{
        ((void(*)(id, SEL, NSNotification*))g_original_occlusion_imp)(self, _cmd, notifCopy);
    });
}

static void _platform_macos_fix_occlusion_crash(void) {
    // Find sokol's _sapp_macos_window_delegate class (created at runtime by sokol)
    Class delegateClass = NSClassFromString(@"_sapp_macos_window_delegate");
    if (!delegateClass) {
        printf("[PLATFORM] macOS: Could not find sokol window delegate class\n");
        return;
    }

    SEL sel = @selector(windowDidChangeOcclusionState:);
    Method method = class_getInstanceMethod(delegateClass, sel);
    if (!method) {
        printf("[PLATFORM] macOS: windowDidChangeOcclusionState: not found, skipping swizzle\n");
        return;
    }

    g_original_occlusion_imp = method_setImplementation(method, (IMP)swizzled_windowDidChangeOcclusionState);
    printf("[PLATFORM] macOS: Swizzled windowDidChangeOcclusionState: (crash fix)\n");
}

// ---- CGEventTap mouse polling ----
// On macOS with Metal, rendering runs on a display-link thread while
// mouseMoved: events are delivered by the main run loop.  When the run
// loop is starved (display-link monopolises the thread), mouse events
// stop arriving for hundreds of milliseconds.
//
// Fix: install a CGEventTap that catches mouse-moved events at the HID
// level and accumulates deltas into a lock-free buffer.  camera code
// polls this buffer once per frame via platform_macos_poll_mouse().

static atomic_int g_tap_dx;
static atomic_int g_tap_dy;
static atomic_bool g_tap_active;  // true when pointer-lock polling is wanted
static atomic_bool g_tap_available; // true if CGEventTap was successfully created
static CFMachPortRef g_tap_port;
static CFRunLoopSourceRef g_tap_source;

static CGEventRef mouse_tap_callback(CGEventTapProxy proxy, CGEventType type,
                                      CGEventRef event, void* refcon) {
    (void)proxy; (void)refcon;
    if (type == kCGEventMouseMoved || type == kCGEventLeftMouseDragged ||
        type == kCGEventRightMouseDragged || type == kCGEventOtherMouseDragged) {
        if (atomic_load(&g_tap_active)) {
            int64_t dx = CGEventGetIntegerValueField(event, kCGMouseEventDeltaX);
            int64_t dy = CGEventGetIntegerValueField(event, kCGMouseEventDeltaY);
            atomic_fetch_add(&g_tap_dx, (int)dx);
            atomic_fetch_add(&g_tap_dy, (int)dy);
        }
    }
    // If the tap is disabled by the system (timeout), re-enable it
    if (type == kCGEventTapDisabledByTimeout && g_tap_port) {
        CGEventTapEnable(g_tap_port, true);
    }
    return event;  // pass event through unchanged
}

// ---- Manual NSEvent polling ----
// Manually drain pending mouse-move NSEvents to work around display-link
// run loop starvation. Only drains mouse movement events to avoid
// interfering with other event types. Safe now that the occlusion crash
// is fixed via the swizzle above.

void platform_macos_poll_events(void) {
    if (!atomic_load(&g_tap_active)) return;  // only poll when mouse is locked
    @autoreleasepool {
        NSEvent* event;
        while ((event = [NSApp nextEventMatchingMask:NSEventMaskMouseMoved |
                                                     NSEventMaskLeftMouseDragged |
                                                     NSEventMaskRightMouseDragged |
                                                     NSEventMaskOtherMouseDragged
                                           untilDate:nil
                                              inMode:NSDefaultRunLoopMode
                                             dequeue:YES])) {
            [NSApp sendEvent:event];
        }
    }
}

void platform_macos_init(void) {
    [NSEvent setMouseCoalescingEnabled:NO];

    // Fix occlusion crash first (must happen after sokol creates the delegate)
    _platform_macos_fix_occlusion_crash();

    // Create a passive event tap for mouse movement
    CGEventMask mask = CGEventMaskBit(kCGEventMouseMoved) |
                       CGEventMaskBit(kCGEventLeftMouseDragged) |
                       CGEventMaskBit(kCGEventRightMouseDragged) |
                       CGEventMaskBit(kCGEventOtherMouseDragged);

    g_tap_port = CGEventTapCreate(kCGSessionEventTap,
                                   kCGHeadInsertEventTap,
                                   kCGEventTapOptionListenOnly,  // passive, no stealing
                                   mask,
                                   mouse_tap_callback,
                                   NULL);
    if (g_tap_port) {
        g_tap_source = CFMachPortCreateRunLoopSource(kCFAllocatorDefault, g_tap_port, 0);
        // Add to the COMMON mode so it fires even when a display-link callback is running
        CFRunLoopAddSource(CFRunLoopGetMain(), g_tap_source, kCFRunLoopCommonModes);
        CGEventTapEnable(g_tap_port, true);
        atomic_store(&g_tap_available, true);
        printf("[PLATFORM] macOS: CGEventTap installed for mouse polling\n");
    } else {
        atomic_store(&g_tap_available, false);
        printf("[PLATFORM] macOS: CGEventTap unavailable, using sokol events + manual polling\n");
    }
    fflush(stdout);
}

// Returns true if CGEventTap is active and working
bool platform_macos_has_tap(void) {
    return atomic_load(&g_tap_available);
}

// Activate/deactivate the tap polling (call when pointer lock changes)
void platform_macos_set_mouse_tap(bool active) {
    atomic_store(&g_tap_active, active);
    if (!active) {
        // Clear any stale deltas
        atomic_store(&g_tap_dx, 0);
        atomic_store(&g_tap_dy, 0);
    }
}

// Drain accumulated deltas (call once per frame from camera_update)
void platform_macos_poll_mouse(float* out_dx, float* out_dy) {
    int dx = atomic_exchange(&g_tap_dx, 0);
    int dy = atomic_exchange(&g_tap_dy, 0);
    *out_dx = (float)dx;
    *out_dy = (float)dy;
}

#endif // __APPLE__

#ifndef SCREENSHOT_H
#define SCREENSHOT_H

// Capture the current backbuffer and save as BMP.
// Currently only implemented for D3D11; no-op on other backends.
void screenshot_capture(void);

#endif

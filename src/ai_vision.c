#include "ai_vision.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define COBJMACROS
#include <d3d11.h>
#include "sokol_app.h"
#include "sokol_gfx.h"
#include "sokol_glue.h"

// stb_image_write for JPEG encoding (already compiled in screenshot.c)
extern int stbi_write_jpg_to_func(void (*func)(void*, void*, int), void* context,
                                   int w, int h, int comp,
                                   const void* data, int quality);

// Base64 encoding
static const char b64_table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static int base64_encode(const unsigned char* src, int src_len, char* dst, int dst_size) {
    int i = 0, j = 0;
    while (i < src_len && j < dst_size - 4) {
        unsigned int a = (i < src_len) ? src[i++] : 0;
        unsigned int b = (i < src_len) ? src[i++] : 0;
        unsigned int c = (i < src_len) ? src[i++] : 0;
        unsigned int triple = (a << 16) | (b << 8) | c;
        dst[j++] = b64_table[(triple >> 18) & 0x3F];
        dst[j++] = b64_table[(triple >> 12) & 0x3F];
        dst[j++] = (i > src_len + 1) ? '=' : b64_table[(triple >> 6) & 0x3F];
        dst[j++] = (i > src_len) ? '=' : b64_table[triple & 0x3F];
    }
    dst[j] = '\0';
    return j;
}

// JPEG write callback — accumulates into a buffer
typedef struct {
    unsigned char* buf;
    int len;
    int cap;
} JpegBuf;

static void jpeg_write_cb(void* context, void* data, int size) {
    JpegBuf* jb = (JpegBuf*)context;
    if (jb->len + size <= jb->cap) {
        memcpy(jb->buf + jb->len, data, size);
        jb->len += size;
    }
}

void ai_vision_init(AiVision* vis, const char* agent_dir) {
    memset(vis, 0, sizeof(AiVision));
    snprintf(vis->agent_dir, sizeof(vis->agent_dir), "%s", agent_dir);
    // Create vision directory
    char path[512];
    snprintf(path, sizeof(path), "%s/vision", agent_dir);
    _mkdir(path);
}

void ai_vision_request_capture(AiVision* vis) {
    vis->capture_requested = true;
}

bool ai_vision_capture(AiVision* vis) {
    if (!vis->capture_requested) return false;
    vis->capture_requested = false;

    sapp_environment env = sapp_get_environment();
    ID3D11Device* device = (ID3D11Device*)env.d3d11.device;
    ID3D11DeviceContext* context = (ID3D11DeviceContext*)env.d3d11.device_context;
    IDXGISwapChain* swapchain = (IDXGISwapChain*)sapp_d3d11_get_swap_chain();

    if (!device || !context || !swapchain) return false;

    ID3D11Texture2D* backbuffer = NULL;
    HRESULT hr = IDXGISwapChain_GetBuffer(swapchain, 0, &IID_ID3D11Texture2D, (void**)&backbuffer);
    if (FAILED(hr) || !backbuffer) return false;

    D3D11_TEXTURE2D_DESC desc;
    ID3D11Texture2D_GetDesc(backbuffer, &desc);

    // Create staging texture
    D3D11_TEXTURE2D_DESC staging_desc = desc;
    staging_desc.Usage = D3D11_USAGE_STAGING;
    staging_desc.BindFlags = 0;
    staging_desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
    staging_desc.MiscFlags = 0;
    staging_desc.SampleDesc.Count = 1;
    staging_desc.SampleDesc.Quality = 0;

    ID3D11Texture2D* staging = NULL;
    hr = ID3D11Device_CreateTexture2D(device, &staging_desc, NULL, &staging);
    if (FAILED(hr) || !staging) {
        ID3D11Texture2D_Release(backbuffer);
        return false;
    }

    if (desc.SampleDesc.Count > 1) {
        D3D11_TEXTURE2D_DESC resolve_desc = desc;
        resolve_desc.SampleDesc.Count = 1;
        resolve_desc.SampleDesc.Quality = 0;
        resolve_desc.Usage = D3D11_USAGE_DEFAULT;
        resolve_desc.BindFlags = 0;
        resolve_desc.CPUAccessFlags = 0;

        ID3D11Texture2D* resolve_tex = NULL;
        hr = ID3D11Device_CreateTexture2D(device, &resolve_desc, NULL, &resolve_tex);
        if (FAILED(hr)) {
            ID3D11Texture2D_Release(staging);
            ID3D11Texture2D_Release(backbuffer);
            return false;
        }
        ID3D11DeviceContext_ResolveSubresource(context,
            (ID3D11Resource*)resolve_tex, 0,
            (ID3D11Resource*)backbuffer, 0, desc.Format);
        ID3D11DeviceContext_CopyResource(context,
            (ID3D11Resource*)staging, (ID3D11Resource*)resolve_tex);
        ID3D11Texture2D_Release(resolve_tex);
    } else {
        ID3D11DeviceContext_CopyResource(context,
            (ID3D11Resource*)staging, (ID3D11Resource*)backbuffer);
    }
    ID3D11Texture2D_Release(backbuffer);

    D3D11_MAPPED_SUBRESOURCE mapped;
    hr = ID3D11DeviceContext_Map(context, (ID3D11Resource*)staging, 0,
        D3D11_MAP_READ, 0, &mapped);
    if (FAILED(hr)) {
        ID3D11Texture2D_Release(staging);
        return false;
    }

    int w = (int)desc.Width;
    int h = (int)desc.Height;

    // Downscale to max 512px wide for API efficiency
    int scale = 1;
    while (w / scale > 512) scale++;
    int sw = w / scale;
    int sh = h / scale;

    // Convert to RGB (downscaled)
    unsigned char* rgb = (unsigned char*)malloc((size_t)sw * sh * 3);
    const unsigned char* src = (const unsigned char*)mapped.pData;
    bool is_rgba = (desc.Format == DXGI_FORMAT_R8G8B8A8_UNORM ||
                    desc.Format == DXGI_FORMAT_R8G8B8A8_UNORM_SRGB);

    for (int y = 0; y < sh; y++) {
        const unsigned char* row = src + (y * scale) * mapped.RowPitch;
        for (int x = 0; x < sw; x++) {
            int sx = x * scale;
            unsigned char* dst = &rgb[(y * sw + x) * 3];
            if (is_rgba) {
                dst[0] = row[sx * 4 + 0];
                dst[1] = row[sx * 4 + 1];
                dst[2] = row[sx * 4 + 2];
            } else {
                dst[0] = row[sx * 4 + 2];
                dst[1] = row[sx * 4 + 1];
                dst[2] = row[sx * 4 + 0];
            }
        }
    }

    ID3D11DeviceContext_Unmap(context, (ID3D11Resource*)staging, 0);
    ID3D11Texture2D_Release(staging);

    // Encode to JPEG in memory
    unsigned char* jpeg_buf = (unsigned char*)malloc(sw * sh * 3);
    JpegBuf jb = { jpeg_buf, 0, sw * sh * 3 };
    stbi_write_jpg_to_func(jpeg_write_cb, &jb, sw, sh, 3, rgb, 70);

    // Save JPEG to file
    vis->capture_id++;
    snprintf(vis->last_image_path, sizeof(vis->last_image_path),
             "%s/vision/view_%04d.jpg", vis->agent_dir, vis->capture_id);
    FILE* f = fopen(vis->last_image_path, "wb");
    if (f) {
        fwrite(jpeg_buf, 1, jb.len, f);
        fclose(f);
    }

    // Base64 encode for API
    vis->base64_len = base64_encode(jpeg_buf, jb.len,
                                     vis->base64_buf, AI_VISION_MAX_BASE64);
    vis->capture_ready = true;

    printf("[VISION] Captured %dx%d → %dx%d, JPEG %d bytes, base64 %d bytes\n",
           w, h, sw, sh, jb.len, vis->base64_len);
    printf("[VISION] Saved to %s\n", vis->last_image_path);
    fflush(stdout);

    free(rgb);
    free(jpeg_buf);
    return true;
}

#else
// Non-Windows stubs
void ai_vision_init(AiVision* vis, const char* agent_dir) {
    memset(vis, 0, sizeof(AiVision));
    snprintf(vis->agent_dir, sizeof(vis->agent_dir), "%s", agent_dir);
}
void ai_vision_request_capture(AiVision* vis) { (void)vis; }
bool ai_vision_capture(AiVision* vis) { (void)vis; return false; }
#endif

bool ai_vision_ready(const AiVision* vis) { return vis->capture_ready; }
void ai_vision_clear(AiVision* vis) { vis->capture_ready = false; vis->base64_len = 0; }

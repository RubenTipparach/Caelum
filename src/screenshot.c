#include "screenshot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <dirent.h>
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

static int screenshot_counter = -1;  // -1 = not initialized

// Ensure screenshots directory exists
static void ensure_screenshots_dir(void) {
    struct stat st;
    if (stat("screenshots", &st) != 0) {
        #ifdef _WIN32
        _mkdir("screenshots");
        #else
        mkdir("screenshots", 0755);
        #endif
    }
}

// Scan screenshots/ for the highest existing number
static int find_highest_screenshot(void) {
    int highest = 0;
#ifdef _WIN32
    WIN32_FIND_DATAA fd;
    HANDLE hFind = FindFirstFileA("screenshots\\screenshot_*.png", &fd);
    if (hFind != INVALID_HANDLE_VALUE) {
        do {
            int num = 0;
            if (sscanf(fd.cFileName, "screenshot_%d.png", &num) == 1) {
                if (num > highest) highest = num;
            }
        } while (FindNextFileA(hFind, &fd));
        FindClose(hFind);
    }
#else
    DIR* dir = opendir("screenshots");
    if (dir) {
        struct dirent* ent;
        while ((ent = readdir(dir)) != NULL) {
            int num = 0;
            if (sscanf(ent->d_name, "screenshot_%d.png", &num) == 1) {
                if (num > highest) highest = num;
            }
        }
        closedir(dir);
    }
#endif
    return highest;
}

// Save RGB buffer to numbered PNG file
static void save_screenshot_png(unsigned char* rgb, int w, int h) {
    ensure_screenshots_dir();

    if (screenshot_counter < 0) {
        screenshot_counter = find_highest_screenshot();
    }

    char filename[128];
    screenshot_counter++;
    snprintf(filename, sizeof(filename), "screenshots/screenshot_%04d.png", screenshot_counter);

    int stride = w * 3;
    if (stbi_write_png(filename, w, h, 3, rgb, stride)) {
        printf("[SCREENSHOT] Saved %s (%dx%d)\n", filename, w, h);
    } else {
        printf("[SCREENSHOT] Failed: could not write %s\n", filename);
    }
    fflush(stdout);
}

// ---- D3D11 backend (Windows) ----
#ifdef SOKOL_D3D11

#include "sokol_app.h"

#define COBJMACROS
#include <d3d11.h>

void screenshot_capture(void) {
    sapp_environment env = sapp_get_environment();
    ID3D11Device* device = (ID3D11Device*)env.d3d11.device;
    ID3D11DeviceContext* context = (ID3D11DeviceContext*)env.d3d11.device_context;

    const void* rtv_ptr = sapp_d3d11_get_swap_chain();

    if (!device || !context) {
        printf("[SCREENSHOT] Failed: could not get D3D11 device/context\n");
        fflush(stdout);
        return;
    }

    IDXGISwapChain* swapchain = (IDXGISwapChain*)rtv_ptr;
    if (!swapchain) {
        printf("[SCREENSHOT] Failed: could not get swap chain\n");
        fflush(stdout);
        return;
    }

    ID3D11Texture2D* backbuffer = NULL;
    HRESULT hr = IDXGISwapChain_GetBuffer(swapchain, 0, &IID_ID3D11Texture2D, (void**)&backbuffer);
    if (FAILED(hr) || !backbuffer) {
        printf("[SCREENSHOT] Failed: could not get backbuffer (0x%08lx)\n", hr);
        fflush(stdout);
        return;
    }

    D3D11_TEXTURE2D_DESC desc;
    ID3D11Texture2D_GetDesc(backbuffer, &desc);

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
        printf("[SCREENSHOT] Failed: could not create staging texture (0x%08lx)\n", hr);
        fflush(stdout);
        ID3D11Texture2D_Release(backbuffer);
        return;
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
        if (FAILED(hr) || !resolve_tex) {
            printf("[SCREENSHOT] Failed: could not create resolve texture\n");
            fflush(stdout);
            ID3D11Texture2D_Release(staging);
            ID3D11Texture2D_Release(backbuffer);
            return;
        }

        ID3D11DeviceContext_ResolveSubresource(context,
            (ID3D11Resource*)resolve_tex, 0,
            (ID3D11Resource*)backbuffer, 0,
            desc.Format);
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
        printf("[SCREENSHOT] Failed: could not map staging texture (0x%08lx)\n", hr);
        fflush(stdout);
        ID3D11Texture2D_Release(staging);
        return;
    }

    int w = (int)desc.Width;
    int h = (int)desc.Height;

    unsigned char* rgb = (unsigned char*)malloc((size_t)w * h * 3);
    const unsigned char* src = (const unsigned char*)mapped.pData;

    for (int y = 0; y < h; y++) {
        const unsigned char* row = src + y * mapped.RowPitch;
        for (int x = 0; x < w; x++) {
            unsigned char* dst = &rgb[(y * w + x) * 3];
            if (desc.Format == DXGI_FORMAT_R8G8B8A8_UNORM ||
                desc.Format == DXGI_FORMAT_R8G8B8A8_UNORM_SRGB) {
                dst[0] = row[x * 4 + 0];
                dst[1] = row[x * 4 + 1];
                dst[2] = row[x * 4 + 2];
            } else {
                dst[0] = row[x * 4 + 2];
                dst[1] = row[x * 4 + 1];
                dst[2] = row[x * 4 + 0];
            }
        }
    }

    ID3D11DeviceContext_Unmap(context, (ID3D11Resource*)staging, 0);
    ID3D11Texture2D_Release(staging);

    save_screenshot_png(rgb, w, h);
    free(rgb);
}

// ---- Metal backend (macOS) ----
#elif defined(SOKOL_METAL)

// Implemented in screenshot_metal.m (Objective-C)
extern unsigned char* screenshot_metal_capture(int* out_w, int* out_h);

void screenshot_capture(void) {
    int w = 0, h = 0;
    unsigned char* rgb = screenshot_metal_capture(&w, &h);
    if (!rgb) {
        printf("[SCREENSHOT] Failed: could not capture Metal framebuffer\n");
        fflush(stdout);
        return;
    }

    save_screenshot_png(rgb, w, h);
    free(rgb);
}

// ---- Other backends: no-op ----
#else
void screenshot_capture(void) {
    printf("[SCREENSHOT] Not implemented on this backend\n");
    fflush(stdout);
}
#endif

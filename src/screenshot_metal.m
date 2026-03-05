// screenshot_metal.m — Metal backend screenshot capture for macOS
#ifdef __APPLE__

#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#include "sokol_app.h"
#include <stdlib.h>

// Called from screenshot.c — captures the current Metal drawable to an RGB buffer.
// Returns malloc'd RGB buffer (caller must free), or NULL on failure.
// Sets *out_w and *out_h to pixel dimensions.
unsigned char* screenshot_metal_capture(int* out_w, int* out_h) {
    sapp_swapchain swapchain = sapp_get_swapchain();
    id<CAMetalDrawable> drawable = (__bridge id<CAMetalDrawable>)swapchain.metal.current_drawable;
    if (!drawable) return NULL;

    id<MTLTexture> texture = drawable.texture;
    if (!texture) return NULL;

    NSUInteger w = texture.width;
    NSUInteger h = texture.height;
    *out_w = (int)w;
    *out_h = (int)h;

    // If MSAA, use the msaa_color_texture resolve target instead
    if (swapchain.metal.msaa_color_texture) {
        id<MTLTexture> msaa_tex = (__bridge id<MTLTexture>)swapchain.metal.msaa_color_texture;
        // We need to resolve MSAA to a non-MSAA texture first
        sapp_environment env = sapp_get_environment();
        id<MTLDevice> device = (__bridge id<MTLDevice>)env.metal.device;

        MTLTextureDescriptor* desc = [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:msaa_tex.pixelFormat
                                                                                       width:w
                                                                                      height:h
                                                                                   mipmapped:NO];
        desc.usage = MTLTextureUsageShaderRead;
        desc.storageMode = MTLStorageModeShared;

        id<MTLTexture> resolve_tex = [device newTextureWithDescriptor:desc];
        if (!resolve_tex) return NULL;

        id<MTLCommandQueue> queue = [device newCommandQueue];
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
        [blit copyFromTexture:texture
                  sourceSlice:0 sourceLevel:0
                 sourceOrigin:MTLOriginMake(0, 0, 0)
                   sourceSize:MTLSizeMake(w, h, 1)
                    toTexture:resolve_tex
             destinationSlice:0 destinationLevel:0
            destinationOrigin:MTLOriginMake(0, 0, 0)];
        [blit endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];

        texture = resolve_tex;
    }

    // Read pixels from the texture
    NSUInteger bytes_per_row = w * 4;
    unsigned char* rgba = (unsigned char*)malloc(bytes_per_row * h);
    if (!rgba) return NULL;

    // For shared/managed textures we can read directly
    // For private textures, we need a blit to a shared texture first
    if (texture.storageMode == MTLStorageModePrivate) {
        sapp_environment env = sapp_get_environment();
        id<MTLDevice> device = (__bridge id<MTLDevice>)env.metal.device;

        MTLTextureDescriptor* desc = [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:texture.pixelFormat
                                                                                       width:w
                                                                                      height:h
                                                                                   mipmapped:NO];
        desc.storageMode = MTLStorageModeShared;
        desc.usage = MTLTextureUsageShaderRead;

        id<MTLTexture> shared_tex = [device newTextureWithDescriptor:desc];
        if (!shared_tex) { free(rgba); return NULL; }

        id<MTLCommandQueue> queue = [device newCommandQueue];
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLBlitCommandEncoder> blit = [cmd blitCommandEncoder];
        [blit copyFromTexture:texture
                  sourceSlice:0 sourceLevel:0
                 sourceOrigin:MTLOriginMake(0, 0, 0)
                   sourceSize:MTLSizeMake(w, h, 1)
                    toTexture:shared_tex
             destinationSlice:0 destinationLevel:0
            destinationOrigin:MTLOriginMake(0, 0, 0)];
        [blit endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];

        texture = shared_tex;
    }

    [texture getBytes:rgba bytesPerRow:bytes_per_row
           fromRegion:MTLRegionMake2D(0, 0, w, h) mipmapLevel:0];

    // Convert BGRA to RGB
    unsigned char* rgb = (unsigned char*)malloc(w * h * 3);
    if (!rgb) { free(rgba); return NULL; }

    for (NSUInteger i = 0; i < w * h; i++) {
        MTLPixelFormat fmt = texture.pixelFormat;
        if (fmt == MTLPixelFormatBGRA8Unorm || fmt == MTLPixelFormatBGRA8Unorm_sRGB) {
            rgb[i * 3 + 0] = rgba[i * 4 + 2]; // R from B
            rgb[i * 3 + 1] = rgba[i * 4 + 1]; // G
            rgb[i * 3 + 2] = rgba[i * 4 + 0]; // B from R
        } else {
            rgb[i * 3 + 0] = rgba[i * 4 + 0]; // R
            rgb[i * 3 + 1] = rgba[i * 4 + 1]; // G
            rgb[i * 3 + 2] = rgba[i * 4 + 2]; // B
        }
    }

    free(rgba);
    return rgb;
}

#endif // __APPLE__

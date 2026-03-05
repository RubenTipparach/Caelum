#!/usr/bin/env bash
set -euo pipefail

echo "=== Hex Planets: Native Build ==="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Compile shaders if sokol-shdc is available
if [ -x "$SCRIPT_DIR/tools/sokol-shdc" ]; then
    SHADER_SLANG="hlsl5:glsl430:glsl300es:metal_macos:wgsl"
    NEEDS_COMPILE=false
    for f in shaders/*.glsl; do
        if [ "$f" -nt "${f}.h" ]; then
            NEEDS_COMPILE=true
            break
        fi
    done
    if $NEEDS_COMPILE; then
        echo "[0/3] Recompiling shaders..."
        for f in shaders/*.glsl; do
            "$SCRIPT_DIR/tools/sokol-shdc" --input "$f" --output "${f}.h" --slang "$SHADER_SLANG"
        done
    fi
fi

# Configure if needed
if [ ! -f build/CMakeCache.txt ]; then
    echo "[1/3] Configuring CMake..."
    cmake -B build -DCMAKE_BUILD_TYPE=Release
else
    echo "[1/3] CMake already configured"
fi

# Build
echo "[2/3] Building (Release)..."
cmake --build build --config Release

# Run
echo "[3/3] Launching hex-planets..."
if [ -f build/hex-planets ]; then
    ./build/hex-planets
elif [ -f build/Release/hex-planets ]; then
    ./build/Release/hex-planets
else
    echo "ERROR: hex-planets binary not found in build/"
    exit 1
fi

#!/usr/bin/env bash
set -euo pipefail

echo "=== Hex Planets: Web Build (Emscripten) ==="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check for emcc, try to activate local emsdk if not found
if ! command -v emcc &>/dev/null; then
    if [ -f "$SCRIPT_DIR/emsdk/emsdk_env.sh" ]; then
        echo "Activating local emsdk..."
        source "$SCRIPT_DIR/emsdk/emsdk_env.sh"
    fi
fi

if ! command -v emcc &>/dev/null; then
    echo "ERROR: emcc not found. Run ./setup.sh first."
    exit 1
fi

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

# Configure
echo "[1/3] Configuring CMake with Emscripten..."
emcmake cmake -B build-web -DCMAKE_BUILD_TYPE=Release

# Build
echo "[2/3] Building..."
cmake --build build-web

echo "[3/3] Web build complete! Output in dist/"
echo ""

# Serve locally
echo "Starting server on http://localhost:8080 ..."
echo "Press Ctrl+C to stop."
echo ""
cd dist

# Try to open browser
case "$(uname -s)" in
    Darwin) open "http://localhost:8080" ;;
    Linux)  xdg-open "http://localhost:8080" 2>/dev/null || true ;;
esac

python3 -m http.server 8080

#!/usr/bin/env bash
set -euo pipefail

echo "=== Hex Planets: Dependencies, Shader Tools & emsdk Setup ==="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- System dependencies ---
echo "[1/5] Checking system dependencies..."

install_deps_macos() {
    if ! command -v brew &>/dev/null; then
        echo "ERROR: Homebrew not found. Install from https://brew.sh"
        exit 1
    fi
    local pkgs=()
    command -v cmake &>/dev/null || pkgs+=(cmake)
    command -v git   &>/dev/null || pkgs+=(git)
    command -v python3 &>/dev/null || pkgs+=(python3)
    if [ ${#pkgs[@]} -gt 0 ]; then
        echo "  Installing: ${pkgs[*]}"
        brew install "${pkgs[@]}"
    else
        echo "  All system dependencies present."
    fi
}

install_deps_linux() {
    local pkgs=()
    command -v cmake   &>/dev/null || pkgs+=(cmake)
    command -v git     &>/dev/null || pkgs+=(git)
    command -v python3 &>/dev/null || pkgs+=(python3)

    # X11/GL dev libs needed for Sokol on Linux
    local dev_pkgs=(libx11-dev libxi-dev libxcursor-dev libgl-dev)

    if [ ${#pkgs[@]} -gt 0 ] || [ ${#dev_pkgs[@]} -gt 0 ]; then
        echo "  Installing: ${pkgs[*]} ${dev_pkgs[*]}"
        sudo apt-get update -qq
        sudo apt-get install -y "${pkgs[@]}" "${dev_pkgs[@]}"
    else
        echo "  All system dependencies present."
    fi
}

case "$(uname -s)" in
    Darwin) install_deps_macos ;;
    Linux)  install_deps_linux ;;
    *)      echo "WARNING: Unsupported OS — install cmake, git, python3 manually." ;;
esac

# --- Git submodules ---
echo "[2/5] Initializing git submodules..."
git submodule update --init --recursive

# --- sokol-shdc (shader compiler) ---
echo "[3/5] Setting up sokol-shdc shader compiler..."

TOOLS_DIR="$SCRIPT_DIR/tools"
mkdir -p "$TOOLS_DIR"

if [ -x "$TOOLS_DIR/sokol-shdc" ]; then
    echo "  sokol-shdc already installed."
else
    SHDC_BASE="https://raw.githubusercontent.com/floooh/sokol-tools-bin/master/bin"
    case "$(uname -s)-$(uname -m)" in
        Darwin-arm64)  SHDC_URL="$SHDC_BASE/osx_arm64/sokol-shdc" ;;
        Darwin-x86_64) SHDC_URL="$SHDC_BASE/osx/sokol-shdc" ;;
        Linux-x86_64)  SHDC_URL="$SHDC_BASE/linux/sokol-shdc" ;;
        Linux-aarch64) SHDC_URL="$SHDC_BASE/linux_arm64/sokol-shdc" ;;
        *)             echo "WARNING: No sokol-shdc binary for $(uname -s)-$(uname -m), skipping."; SHDC_URL="" ;;
    esac
    if [ -n "${SHDC_URL:-}" ]; then
        echo "  Downloading sokol-shdc..."
        curl -L -o "$TOOLS_DIR/sokol-shdc" "$SHDC_URL"
        chmod +x "$TOOLS_DIR/sokol-shdc"
    fi
fi

# Compile shaders if sokol-shdc is available
if [ -x "$TOOLS_DIR/sokol-shdc" ]; then
    SHADER_SLANG="hlsl5:glsl430:glsl300es:metal_macos:wgsl"
    echo "  Compiling shaders (slang: $SHADER_SLANG)..."
    for f in "$SCRIPT_DIR"/shaders/*.glsl; do
        echo "    $(basename "$f")"
        "$TOOLS_DIR/sokol-shdc" --input "$f" --output "${f}.h" --slang "$SHADER_SLANG"
    done
    echo "  All shaders compiled."
fi

# --- emsdk ---
echo "[4/5] Setting up Emscripten SDK..."

if [ -f "$SCRIPT_DIR/emsdk/emsdk" ]; then
    echo "  emsdk already installed at $SCRIPT_DIR/emsdk"
    echo "  Activating..."
    cd "$SCRIPT_DIR/emsdk"
    ./emsdk activate latest
    source ./emsdk_env.sh
else
    echo "  Cloning emsdk..."
    git clone https://github.com/emscripten-core/emsdk.git "$SCRIPT_DIR/emsdk"
    cd "$SCRIPT_DIR/emsdk"
    echo "  Installing latest Emscripten..."
    ./emsdk install latest
    echo "  Activating..."
    ./emsdk activate latest
    source ./emsdk_env.sh
fi

cd "$SCRIPT_DIR"

# --- Verify ---
echo "[5/5] Verifying installation..."
echo "  cmake:      $(cmake --version | head -1)"
echo "  git:        $(git --version)"
if [ -x "$SCRIPT_DIR/tools/sokol-shdc" ]; then
    echo "  sokol-shdc: $("$SCRIPT_DIR/tools/sokol-shdc" --help 2>&1 | head -1)"
fi
echo "  emcc:       $(emcc --version | head -1)"

echo ""
echo "=== Setup complete! ==="
echo "  Native build:  ./build-run.sh"
echo "  Web build:     ./build-web.sh"

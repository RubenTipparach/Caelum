@echo off
setlocal

REM ============================================================
REM ai-install.bat — Downloads llama-server and Qwen3-4B Q4_K_M
REM Run this once to set up local AI for hex-planets.
REM Files are placed in tools/ai/ (git-ignored).
REM ============================================================

set AI_DIR=%~dp0ai
set LLAMA_VERSION=b8468
set MODEL_REPO=Qwen
set MODEL_NAME=Qwen3-4B-GGUF
set MODEL_FILE=Qwen3-4B-Q4_K_M.gguf

REM llama.cpp release binary (Windows CUDA x64 — requires NVIDIA GPU)
set LLAMA_ZIP=llama-%LLAMA_VERSION%-bin-win-cuda-12.4-x64.zip
set LLAMA_URL=https://github.com/ggml-org/llama.cpp/releases/download/%LLAMA_VERSION%/%LLAMA_ZIP%

REM CUDA runtime DLLs (required alongside the CUDA build)
set CUDART_ZIP=cudart-llama-bin-win-cuda-12.4-x64.zip
set CUDART_URL=https://github.com/ggml-org/llama.cpp/releases/download/%LLAMA_VERSION%/%CUDART_ZIP%

REM Hugging Face model URL
set MODEL_URL=https://huggingface.co/%MODEL_REPO%/%MODEL_NAME%/resolve/main/%MODEL_FILE%

echo.
echo === Hex Planets AI Setup ===
echo.
echo This will download:
echo   1. llama-server (%LLAMA_VERSION%)    ~50 MB
echo   2. Qwen3-4B Q4_K_M model           ~2.8 GB
echo.
echo Install directory: %AI_DIR%
echo.

if not exist "%AI_DIR%" mkdir "%AI_DIR%"

REM ---- Step 1: Download llama.cpp release ----
if exist "%AI_DIR%\llama-server.exe" (
    echo [OK] llama-server.exe already exists, skipping.
    goto step2
)

echo [1/2] Downloading llama.cpp %LLAMA_VERSION%...
curl -L -o "%AI_DIR%\%LLAMA_ZIP%" "%LLAMA_URL%"
if errorlevel 1 (
    echo [ERROR] Failed to download llama.cpp. Check your internet connection.
    exit /b 1
)

echo Extracting...
if not exist "%AI_DIR%\llama-release" mkdir "%AI_DIR%\llama-release"
tar -xf "%AI_DIR%\%LLAMA_ZIP%" -C "%AI_DIR%\llama-release"

REM Copy llama-server.exe to ai/ root (flat zip layout since b8xxx)
copy "%AI_DIR%\llama-release\llama-server.exe" "%AI_DIR%\llama-server.exe" >nul 2>&1
if not exist "%AI_DIR%\llama-server.exe" (
    REM Try nested path layout (older releases)
    for /r "%AI_DIR%\llama-release" %%f in (llama-server.exe) do (
        copy "%%f" "%AI_DIR%\llama-server.exe" >nul
    )
)

REM Copy required DLLs next to the exe
for /r "%AI_DIR%\llama-release" %%f in (*.dll) do (
    copy "%%f" "%AI_DIR%\" >nul 2>&1
)

REM Cleanup zip and extracted folder
del "%AI_DIR%\%LLAMA_ZIP%" >nul 2>&1
rmdir /s /q "%AI_DIR%\llama-release" >nul 2>&1

if not exist "%AI_DIR%\llama-server.exe" (
    echo [ERROR] Could not find llama-server.exe in release archive.
    echo         You may need to extract manually from: %LLAMA_URL%
    exit /b 1
)
echo [OK] llama-server.exe installed.

REM ---- Step 1b: Download CUDA runtime DLLs ----
echo [1b/2] Downloading CUDA runtime DLLs...
curl -L -o "%AI_DIR%\%CUDART_ZIP%" "%CUDART_URL%"
if errorlevel 1 (
    echo [WARN] Failed to download CUDA runtime. GPU acceleration may not work.
    goto step2
)
echo Extracting CUDA runtime...
if not exist "%AI_DIR%\cudart-release" mkdir "%AI_DIR%\cudart-release"
tar -xf "%AI_DIR%\%CUDART_ZIP%" -C "%AI_DIR%\cudart-release"
for /r "%AI_DIR%\cudart-release" %%f in (*.dll) do (
    copy "%%f" "%AI_DIR%\" >nul 2>&1
)
del "%AI_DIR%\%CUDART_ZIP%" >nul 2>&1
rmdir /s /q "%AI_DIR%\cudart-release" >nul 2>&1
echo [OK] CUDA runtime installed.

REM ---- Step 2: Download model ----
:step2
if exist "%AI_DIR%\%MODEL_FILE%" (
    echo [OK] %MODEL_FILE% already exists, skipping.
    goto step3
)

echo [2/2] Downloading %MODEL_FILE% (~2.8 GB, this may take a while)...
curl -L -o "%AI_DIR%\%MODEL_FILE%" "%MODEL_URL%"
if errorlevel 1 (
    echo [ERROR] Failed to download model. Check your internet connection.
    exit /b 1
)
echo [OK] Model downloaded.

REM ---- Step 3: Copy GBNF grammar ----
:step3
if exist "%~dp0ai-grammar.gbnf" (
    copy "%~dp0ai-grammar.gbnf" "%AI_DIR%\grammar.gbnf" >nul
    echo [OK] Grammar file copied.
)

echo.
echo === Setup Complete ===
echo.
echo To test manually:
echo   cd tools\ai
echo   llama-server.exe -m %MODEL_FILE% --port 8080 -ngl 99
echo.
echo The game will launch this automatically when AI is enabled.
echo.

endlocal

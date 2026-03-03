@echo off
setlocal enabledelayedexpansion

echo === Hex Planets: Web Build (Emscripten) ===

:: Check for emcc, try to activate local emsdk if not found
where emcc >nul 2>&1
if errorlevel 1 (
    if exist "%~dp0emsdk\upstream\emscripten\emcc.bat" (
        echo Activating local emsdk...
        set "EMSDK=%~dp0emsdk"
        set "EM_CONFIG=%~dp0emsdk\.emscripten"
        set "EMSDK_PYTHON=%~dp0emsdk\python\3.13.3_64bit\python.exe"
        set "EMSDK_NODE=%~dp0emsdk\node\22.16.0_64bit\bin\node.exe"
        set "PATH=%~dp0emsdk;%~dp0emsdk\upstream\emscripten;%~dp0emsdk\upstream\bin;%~dp0emsdk\node\22.16.0_64bit\bin;%~dp0emsdk\python\3.13.3_64bit;%PATH%"
    )
)
where emcc >nul 2>&1
if errorlevel 1 (
    echo ERROR: emcc not found. Run setup-emsdk.bat first.
    exit /b 1
)

:: Configure
echo [1/3] Configuring CMake with Emscripten...
call emcmake cmake -B build-web -DCMAKE_BUILD_TYPE=Release
if !ERRORLEVEL! neq 0 goto :error

:: Build
echo [2/3] Building...
cmake --build build-web
if !ERRORLEVEL! neq 0 goto :error

echo [3/3] Web build complete! Output in dist/
echo.

:: Serve locally
echo Starting server on http://localhost:8080 ...
echo Press Ctrl+C to stop.
echo.
cd dist
start "" http://localhost:8080
python -m http.server 8080
goto :done

:error
echo.
echo BUILD FAILED
exit /b 1

:done
echo Done.

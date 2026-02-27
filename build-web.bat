@echo off
setlocal

echo === Hex Planets: Web Build (Emscripten) ===

:: Check for emsdk
where emcc >nul 2>&1
if errorlevel 1 (
    echo ERROR: emcc not found. Please activate emsdk first:
    echo   emsdk activate latest
    echo   emsdk_env.bat
    exit /b 1
)

:: Configure
echo [1/3] Configuring CMake with Emscripten...
emcmake cmake -B build-web -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 goto :error

:: Build
echo [2/3] Building...
cmake --build build-web
if errorlevel 1 goto :error

echo [3/3] Web build complete! Output in dist/
echo.

:: Offer to serve locally
echo To test locally, run:
echo   cd dist ^&^& python -m http.server 8080
echo Then open http://localhost:8080 in your browser.
echo.

set /p SERVE="Start local server now? [Y/n] "
if /i "%SERVE%"=="n" goto :done
if /i "%SERVE%"=="N" goto :done

echo Starting server on http://localhost:8080 ...
cd dist
python -m http.server 8080
goto :done

:error
echo.
echo BUILD FAILED
exit /b 1

:done
echo Done.

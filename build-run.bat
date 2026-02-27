@echo off
setlocal

echo === Hex Planets: Windows Build ===

:: Configure if needed
if not exist build\CMakeCache.txt (
    echo [1/3] Configuring CMake...
    cmake -B build -G "Visual Studio 17 2022"
    if errorlevel 1 goto :error
) else (
    echo [1/3] CMake already configured
)

:: Build
echo [2/3] Building (Release)...
cmake --build build --config Release
if errorlevel 1 goto :error

:: Run - app allocates its own console window for logging
echo [3/3] Launching hex-planets.exe
echo   (A separate console window will open with debug logs)
start "" build\Release\hex-planets.exe
goto :done

:error
echo.
echo BUILD FAILED
exit /b 1

:done
echo Done.

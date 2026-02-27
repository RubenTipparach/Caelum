@echo off
setlocal

echo === Hex Planets: Test Suite ===

:: Configure
if not exist build\CMakeCache.txt (
    echo [1/3] Configuring CMake...
    cmake -B build -G "Visual Studio 17 2022"
    if errorlevel 1 goto :error
) else (
    echo [1/3] CMake already configured
)

:: Build test target
echo [2/3] Building tests...
cmake --build build --config Release --target hex-tests
if errorlevel 1 goto :error

:: Run tests
echo [3/3] Running tests...
echo.
build\Release\hex-tests.exe
set TEST_RESULT=%ERRORLEVEL%

if %TEST_RESULT% NEQ 0 (
    echo.
    echo TESTS FAILED
    exit /b 1
)

echo.
echo ALL TESTS PASSED
exit /b 0

:error
echo.
echo BUILD FAILED
exit /b 1

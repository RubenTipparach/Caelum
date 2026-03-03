@echo off
setlocal

echo === Emscripten SDK Setup ===

:: Check if already installed
if exist "%~dp0emsdk\emsdk.bat" (
    echo emsdk already installed at %~dp0emsdk
    echo Activating...
    cd /d "%~dp0emsdk"
    call emsdk activate latest
    call emsdk_env.bat
    echo.
    echo Ready! You can now run build-web.bat
    goto :done
)

:: Check for git
where git >nul 2>&1
if errorlevel 1 (
    echo ERROR: git not found. Please install git first.
    exit /b 1
)

:: Clone and install
echo Cloning emsdk...
git clone https://github.com/emscripten-core/emsdk.git "%~dp0emsdk"
if errorlevel 1 goto :error

cd /d "%~dp0emsdk"

echo Installing latest Emscripten...
call emsdk install latest
if errorlevel 1 goto :error

echo Activating...
call emsdk activate latest
if errorlevel 1 goto :error

call emsdk_env.bat

echo.
echo === Setup complete! ===
echo Run build-web.bat to build the web version.
goto :done

:error
echo.
echo SETUP FAILED
exit /b 1

:done

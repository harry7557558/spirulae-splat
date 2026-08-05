@echo off
setlocal
rem Development build for Windows. Extra arguments are passed to CMake, e.g.
rem   build_develop.bat -DSS_BUILD_CLI=ON -DSS_BUILD_GUI=ON
rem Builds the standalone spirula.exe.
rem Works from a plain cmd prompt: locates VS via vswhere and calls vcvars64,
rem uses the VS-bundled CMake/Ninja when none are on PATH, and picks the
rem newest installed CUDA toolkit unless CUDA_PATH is already set.

cd /d "%~dp0"

rem ---------------------------------------------------------------------------
rem Codegen (optional -- the generated files are committed)
rem ---------------------------------------------------------------------------
where python >nul 2>&1 || goto :skip_codegen
python tools\codegen\generate_headers.py || goto :codegen_warn
python tools\codegen\generate_kernel_instantiation.py || goto :codegen_warn
goto :msvc_env
:codegen_warn
echo codegen failed -- continuing with committed generated files
goto :msvc_env
:skip_codegen
echo python not found -- skipping codegen (using committed generated files)

rem ---------------------------------------------------------------------------
rem MSVC environment. vcvars64 is called even when cl is already on PATH: a
rem machine-wide cl/INCLUDE (e.g. exported from an old developer prompt) can
rem point at an uninstalled Windows SDK, which breaks compiles with missing
rem stddef.h/stdlib.h. Re-initializing is harmless in a developer prompt.
rem ---------------------------------------------------------------------------
:msvc_env
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSROOT="
if exist "%VSWHERE%" for /f "usebackq delims=" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSROOT=%%i"
if not defined VSROOT set "VSROOT=%ProgramFiles%\Microsoft Visual Studio\2022\Community"
if exist "%VSROOT%\VC\Auxiliary\Build\vcvars64.bat" (
    call "%VSROOT%\VC\Auxiliary\Build\vcvars64.bat" >nul
    goto :have_cl
)
where cl >nul 2>&1 && goto :have_cl
echo ERROR: vcvars64.bat not found -- install Visual Studio with C++ tools,
echo        or run this script from a developer command prompt.
exit /b 1
:have_cl

rem ---------------------------------------------------------------------------
rem CMake + Ninja (fall back to the copies bundled with Visual Studio)
rem ---------------------------------------------------------------------------
where cmake >nul 2>&1
if errorlevel 1 set "PATH=%PATH%;%VSROOT%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin"
where ninja >nul 2>&1
if errorlevel 1 set "PATH=%PATH%;%VSROOT%\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"

rem ---------------------------------------------------------------------------
rem CUDA toolkit: newest installed version. An ambient CUDA_PATH is NOT
rem trusted by default (it may pin a toolkit older than the MSVC in use);
rem pass a trailing -DCMAKE_CUDA_COMPILER=... to pick a specific one.
rem ---------------------------------------------------------------------------
set "_CUDA_ROOT="
for /d %%d in ("%ProgramFiles%\NVIDIA GPU Computing Toolkit\CUDA\v*") do set "_CUDA_ROOT=%%d"
if not defined _CUDA_ROOT set "_CUDA_ROOT=%CUDA_PATH%"
set "CUDAARG="
if defined _CUDA_ROOT (
    set "CUDA_PATH=%_CUDA_ROOT%"
    set "PATH=%_CUDA_ROOT%\bin;%PATH%"
    set CUDAARG=-DCMAKE_CUDA_COMPILER="%_CUDA_ROOT%\bin\nvcc.exe"
)

rem ---------------------------------------------------------------------------
rem Configure + build (RAM-aware job count, mirrors build_develop.bash)
rem ---------------------------------------------------------------------------
cmake -G Ninja -B build -DCMAKE_BUILD_TYPE=Release %CUDAARG% %*
if errorlevel 1 exit /b 1

echo.
set JOB_RAM_MB=750
set AVAILABLE_MB=0
for /f %%m in ('powershell -NoProfile -Command "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1024)"') do set AVAILABLE_MB=%%m
set /a JOBS=AVAILABLE_MB/JOB_RAM_MB
if %JOBS% gtr %NUMBER_OF_PROCESSORS% set JOBS=%NUMBER_OF_PROCESSORS%
if %JOBS% lss 1 set JOBS=1
echo Available RAM : %AVAILABLE_MB% MB
echo CPU cores     : %NUMBER_OF_PROCESSORS%
echo Using jobs    : %JOBS%
echo.

cmake --build build -j %JOBS%
if errorlevel 1 exit /b 1

echo.
echo Build complete: build\spirula.exe

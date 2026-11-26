@echo off
rem Run with CALL in a command prompt; this script preserves the compiler environment.
if defined OCFEM_INTEL_COMPILER_ROOT goto compiler_found
set "OCFEM_INTEL_COMPILER_ROOT=C:\Program Files (x86)\Intel\oneAPI\compiler\2025.2"
if exist "%OCFEM_INTEL_COMPILER_ROOT%\env\vars.bat" goto compiler_found
set "OCFEM_INTEL_COMPILER_ROOT=C:\Program Files (x86)\Intel\oneAPI\compiler\latest"
:compiler_found
if not exist "%OCFEM_INTEL_COMPILER_ROOT%\env\vars.bat" (
  echo Set OCFEM_INTEL_COMPILER_ROOT to the installed Intel compiler directory.
  exit /b 1
)
call "%OCFEM_INTEL_COMPILER_ROOT%\env\vars.bat" vs2022
if errorlevel 1 exit /b 1
if not defined MKLROOT set "MKLROOT=C:\Program Files (x86)\Intel\oneAPI\mkl\2025.2"
if not exist "%MKLROOT%\env\vars.bat" set "MKLROOT=C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
if not exist "%MKLROOT%\env\vars.bat" exit /b 1
call "%MKLROOT%\env\vars.bat" intel64
exit /b %errorlevel%

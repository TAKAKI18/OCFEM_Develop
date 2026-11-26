@echo off
setlocal
cd /d "%~dp0\.."
call "%~dp0setup_env_windows.cmd"
if errorlevel 1 exit /b 1
if not defined OCFEM_PYTHON set "OCFEM_PYTHON=python"
set OMP_NUM_THREADS=2
set MKL_NUM_THREADS=2
"%OCFEM_PYTHON%" tools\check_runtime.py
exit /b %errorlevel%

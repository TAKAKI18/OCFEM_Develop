@echo off
setlocal
cd /d "%~dp0\.."
call "%~dp0setup_env_windows.cmd"
if errorlevel 1 exit /b 1
if not defined OCFEM_PYTHON set "OCFEM_PYTHON=python"
"%OCFEM_PYTHON%" tools\verify_sources.py
if errorlevel 1 exit /b 1
"%OCFEM_PYTHON%" tools\build.py %*
exit /b %errorlevel%

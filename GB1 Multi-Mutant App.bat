@echo off
setlocal
cd /d "%~dp0"
title GB1 Multi-Mutant App
echo Starting GB1 Multi-Mutant App...
echo Keep this window open while using the app.
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0server.ps1"
echo.
echo App server stopped.
pause

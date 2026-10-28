@echo off
REM run_fig7c.bat -- Fig 7C in open-source PyMOL (no watermark); log in figures\final\fig7c_log.txt
call conda activate chromatin
cd /d %~dp0
if not exist ..\figures\final mkdir ..\figures\final
pymol -cq fig7c_render.pml > ..\figures\final\fig7c_log.txt 2>&1
echo done, see ..\figures\final\fig7c_log.txt

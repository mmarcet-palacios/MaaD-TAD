@echo off
REM Run once in "x64 Native Tools Command Prompt for VS 2022", from MaaD-TAD\control_loci
cl /nologo /O2 /openmp /EHsc seqsim_windows_fast.cpp /Fe:seqsim_windows_fast.exe
if exist seqsim_windows_fast.exe (echo Built OK) else (echo Build failed)

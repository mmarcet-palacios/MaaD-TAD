@echo off
rem Build the five representative ensemble members to atoms (relax, Builder, Checker).
rem Run from the MaaD-TAD folder after unzipping ensemble_paths_frames.zip here.
for %%m in (m032 m005 m002 m030 m019) do (
  python build_member.py ensemble\%%m --workers 24 %*
)

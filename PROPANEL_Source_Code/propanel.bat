@echo off
cls
setlocal

rem ===== Settings =====
rem Use ifort (classic) by default. To use ifx, replace "ifort" with "ifx" below.
set FC=ifx

rem Flags mapped from Linux: -check bounds -traceback -fltconsistency -fpe0
rem Closest Windows equivalents:
rem   /check:bounds  /traceback  /fp:consistent  /fpe:0
set FFLAGS=/check:all /traceback /fpe:0 /warn:all /Od /fp:precise 

rem ===== Clean =====
if not exist Code mkdir Code
del /q Code\*.* 2>nul
del /q ProPanel2026_v1.0_debug.exe 2>nul

rem ===== Compile (Base) =====
%FC% /c %FFLAGS% Base\propanel_mod.f90 || goto :err
%FC% /c %FFLAGS% Base\ProPanel2026_v1.0.f90 || goto :err
%FC% /c %FFLAGS% Base\delvars.f90 || goto :err
%FC% /c %FFLAGS% Base\progress.f90 || goto :err

rem ===== Compile (Grids) =====
%FC% /c %FFLAGS% Grids\bladegrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\bladewakegrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\nozzlegrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\nozzlewakegrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\nozzledef.f90 || goto :err
%FC% /c %FFLAGS% Grids\hubgrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\geoduct37.f90 || goto :err

rem ===== Compile (Grape) =====
%FC% /c %FFLAGS% Grape\angri.f90 || goto :err
%FC% /c %FFLAGS% Grape\bord.f90 || goto :err
%FC% /c %FFLAGS% Grape\calcb.f90 || goto :err
%FC% /c %FFLAGS% Grape\calphi.f90 || goto :err
%FC% /c %FFLAGS% Grape\coef.f90 || goto :err
%FC% /c %FFLAGS% Grape\grape.f90 || goto :err
%FC% /c %FFLAGS% Grape\guessa.f90 || goto :err
%FC% /c %FFLAGS% Grape\rhs.f90 || goto :err
%FC% /c %FFLAGS% Grape\sip.f90 || goto :err
%FC% /c %FFLAGS% Grape\splin.f90 || goto :err

rem ===== Compile (Calc) =====
%FC% /c %FFLAGS% Calc\linint.f90 || goto :err
%FC% /c %FFLAGS% Calc\intk1.f90 || goto :err
%FC% /c %FFLAGS% Calc\splint.f90 || goto :err
%FC% /c %FFLAGS% Calc\stret_choice.f90 || goto :err
%FC% /c %FFLAGS% Calc\dscal.f90 || goto :err
%FC% /c %FFLAGS% Calc\stret.f90 || goto :err
%FC% /c %FFLAGS% Calc\stret2.f90 || goto :err
%FC% /c %FFLAGS% Calc\sxx.f90 || goto :err
%FC% /c %FFLAGS% Calc\shxx.f90 || goto :err
%FC% /c %FFLAGS% Calc\spline.f90 || goto :err
%FC% /c %FFLAGS% Calc\ispline.f90 || goto :err

rem ===== Compile (Linpack; .f is fixed-form and is detected by extension) =====
%FC% /c %FFLAGS% Linpack\cubspl.f || goto :err
%FC% /c %FFLAGS% Linpack\ppvalu.f || goto :err
%FC% /c %FFLAGS% Linpack\interv.f || goto :err

rem ===== Link to EXE =====
%FC% /nologo /Fe:ProPanel2026_v1.0_debug.exe *.obj || goto :err

rem ===== Stage outputs =====
move /y *.obj Code >nul
move /y *.mod Code >nul 2>nul
move /y *_genmod.f90 Code >nul 2>nul
if exist ProPanel2026_v1.0_debug.exe copy /y ProPanel2026_v1.0_debug.exe Code >nul

echo.
echo Build succeeded: ProPanel2026_v1.0_debug.exe
exit /b 0

:err
echo.
echo Build failed with errorlevel %errorlevel%.
exit /b 1

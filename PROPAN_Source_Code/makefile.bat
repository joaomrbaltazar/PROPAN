@echo off
setlocal

rem ===== Settings =====
rem Use ifort (classic) by default. To use ifx, replace "ifort" with "ifx" below.
set FC=ifx

rem Flags mapped from Linux: -check bounds -traceback -fltconsistency -fpe0
rem Closest Windows equivalents:
rem   /check:bounds  /traceback  /fp:consistent  /fpe:0
rem   /O2 /check:bounds /warn:interfaces /traceback /fpe:0 /Qinit:snan,arrays
set FFLAGS=/check:bounds /traceback /fpe:0 /warn:interfaces /Od /fp:precise

rem ===== Clean =====
if not exist Code mkdir Code
del /q Code\*.* 2>nul
del /q ProPan2025_v1.0.exe 2>nul

rem ===== Compile (Base) =====
%FC% /c %FFLAGS% Base\propan_mod.f90 || goto :err
%FC% /c %FFLAGS% Base\ProPan2025_v1.0.f90 || goto :err
%FC% /c %FFLAGS% Base\progress.f90 || goto :err
%FC% /c %FFLAGS% Base\input.f90 || goto :err
%FC% /c %FFLAGS% Base\inivars.f90 || goto :err
%FC% /c %FFLAGS% Base\delvars.f90 || goto :err
%FC% /c %FFLAGS% Base\counters.f90 || goto :err
%FC% /c %FFLAGS% Base\output.f90 || goto :err

rem ===== Compile (Geom) =====
%FC% /c %FFLAGS% Geom\bladegeom.f90 || goto :err
%FC% /c %FFLAGS% Geom\bladewakegeom.f90 || goto :err
%FC% /c %FFLAGS% Geom\nozzlegeom.f90 || goto :err
%FC% /c %FFLAGS% Geom\nozzlewakegeom.f90 || goto :err
%FC% /c %FFLAGS% Geom\hubgeom.f90 || goto :err
%FC% /c %FFLAGS% Geom\panel.f90 || goto :err
%FC% /c %FFLAGS% Geom\pancoor.f90 || goto :err

rem ===== Compile (Numdif) =====
%FC% /c %FFLAGS% Numdif\numdifblade.f90 || goto :err
%FC% /c %FFLAGS% Numdif\numdifbladewake.f90 || goto :err
%FC% /c %FFLAGS% Numdif\numdifnozzle.f90 || goto :err
%FC% /c %FFLAGS% Numdif\numdifhub.f90 || goto :err

rem ===== Compile (InfCoef) =====
%FC% /c %FFLAGS% InfCoef\InfCoefMatrixStd.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\InfCoefMatrixUnStd.f90 || goto :err

rem ===== Compile (PotCoef) =====
%FC% /c %FFLAGS% InfCoef\PotCoef\bladecoef.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\bladewakecoef.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\nozzlecoef.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\nozzlewakecoef.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\hubcoef.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\potpan.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\potpan_num.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\PotCoef\potlinsub.f90 || goto :err

rem ===== Compile (VelCoef) =====
%FC% /c %FFLAGS% InfCoef\VelCoef\bladevelo.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\bladewakevelo.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\nozzlevelo.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\nozzlewakevelo.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\hubvelo.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\imagehubvelo.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\velpan.f90 || goto :err
%FC% /c %FFLAGS% InfCoef\VelCoef\gaussvel.f90 || goto :err

rem ===== Compile (Solve) =====
%FC% /c %FFLAGS% Solve\SolveRhsWet.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveRhsCav.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveRhsCavRed.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveLkWet.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveLkCav.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveIpkcWet.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveIpkcCav.f90 || goto :err
%FC% /c %FFLAGS% Solve\SolveWake.f90 || goto :err

rem ===== Compile (WakeAlign) =====
%FC% /c %FFLAGS% WakeAlign\wakealign1.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\wakealign2.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\nozzledef.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\geoduct37.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\bladewakedisp.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\nozzledisp.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\nozzlewakedisp.f90 || goto :err
%FC% /c %FFLAGS% WakeAlign\ff.f90 || goto :err

rem ===== Compile (Cav) =====
%FC% /c %FFLAGS% Cav\cavcheck.f90 || goto :err
%FC% /c %FFLAGS% Cav\cavprop.f90 || goto :err
%FC% /c %FFLAGS% Cav\caverrc.f90 || goto :err
%FC% /c %FFLAGS% Cav\cavpotp.f90 || goto :err
%FC% /c %FFLAGS% Cav\cavpots.f90 || goto :err
%FC% /c %FFLAGS% Cav\cavthickp.f90 || goto :err
%FC% /c %FFLAGS% Cav\cavthicks.f90 || goto :err
%FC% /c %FFLAGS% Cav\cavrecover.f90 || goto :err

rem ===== Compile (VtCp) =====
%FC% /c %FFLAGS% VtCp\velp.f90 || goto :err
%FC% /c %FFLAGS% VtCp\velpw.f90 || goto :err
%FC% /c %FFLAGS% VtCp\presp.f90 || goto :err
%FC% /c %FFLAGS% VtCp\veln.f90 || goto :err
%FC% /c %FFLAGS% VtCp\presn.f90 || goto :err
%FC% /c %FFLAGS% VtCp\presnte.f90 || goto :err
%FC% /c %FFLAGS% VtCp\velh.f90 || goto :err
%FC% /c %FFLAGS% VtCp\presh.f90 || goto :err

rem ===== Compile (Wake) =====
%FC% /c %FFLAGS% Calc\Wake\vwake.f90 || goto :err
%FC% /c %FFLAGS% Calc\Wake\fourier_coef.f90 || goto :err
%FC% /c %FFLAGS% Calc\Wake\fourier_function.f90 || goto :err

rem ===== Compile (Grids) =====
%FC% /c %FFLAGS% Grids\hubgrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\nozzlegrid.f90 || goto :err
%FC% /c %FFLAGS% Grids\nozzlewakegrid.f90 || goto :err

rem ===== Compile (Grape) =====
%FC% /c %FFLAGS% Grids\Grape\angri.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\bord.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\calcb.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\calphi.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\coef.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\grape.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\guessa.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\rhs.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\sip.f90 || goto :err
%FC% /c %FFLAGS% Grids\Grape\splin.f90 || goto :err

rem ===== Compile (Field) =====
%FC% /c %FFLAGS% Field\bladecoeff.f90 || goto :err
%FC% /c %FFLAGS% Field\bladewakecoeff.f90 || goto :err
%FC% /c %FFLAGS% Field\hubcoeff.f90 || goto :err
%FC% /c %FFLAGS% Field\nozzlecoeff.f90 || goto :err
%FC% /c %FFLAGS% Field\nozzlewakecoeff.f90 || goto :err
%FC% /c %FFLAGS% Field\velfStd.f90 || goto :err
%FC% /c %FFLAGS% Field\velfUnStd.f90 || goto :err
%FC% /c %FFLAGS% Field\presfStd.f90 || goto :err
%FC% /c %FFLAGS% Field\presfUnStd.f90 || goto :err

rem ===== Compile (Calc) =====
%FC% /c %FFLAGS% Calc\JacobianWet.f90 || goto :err
%FC% /c %FFLAGS% Calc\JacobianCav.f90 || goto :err
%FC% /c %FFLAGS% Calc\linint.f90 || goto :err
%FC% /c %FFLAGS% Calc\intk1.f90 || goto :err
%FC% /c %FFLAGS% Calc\splint.f90 || goto :err
%FC% /c %FFLAGS% Calc\spline.f90 || goto :err
%FC% /c %FFLAGS% Calc\ispline.f90 || goto :err
%FC% /c %FFLAGS% Calc\gaussint.f90 || goto :err
%FC% /c %FFLAGS% Calc\bisof.f || goto :err
%FC% /c %FFLAGS% Calc\CtCq.f90 || goto :err
%FC% /c %FFLAGS% Calc\gaperrg.f90 || goto :err
%FC% /c %FFLAGS% Calc\periodicflow.f90 || goto :err
%FC% /c %FFLAGS% Calc\velinf.f90 || goto :err
%FC% /c %FFLAGS% Calc\stret2.f90 || goto :err
%FC% /c %FFLAGS% Calc\sxx.f90 || goto :err
%FC% /c %FFLAGS% Calc\shxx.f90 || goto :err
%FC% /c %FFLAGS% Calc\sdiv.f90 || goto :err
%FC% /c %FFLAGS% Calc\frame.f90 || goto :err

rem ===== Compile (Linpack) =====
%FC% /c %FFLAGS% Linpack\cubspl.f || goto :err
%FC% /c %FFLAGS% Linpack\daxpy.f || goto :err
%FC% /c %FFLAGS% Linpack\dcopy.f || goto :err
%FC% /c %FFLAGS% Linpack\ddot.f || goto :err
%FC% /c %FFLAGS% Linpack\dgedi.f || goto :err
%FC% /c %FFLAGS% Linpack\dgefa.f || goto :err
%FC% /c %FFLAGS% Linpack\dgemm.f || goto :err
%FC% /c %FFLAGS% Linpack\dger.f || goto :err
%FC% /c %FFLAGS% Linpack\dgesl.f || goto :err
%FC% /c %FFLAGS% Linpack\drotm.f || goto :err
%FC% /c %FFLAGS% Linpack\drotmg.f || goto :err
%FC% /c %FFLAGS% Linpack\dscal.f || goto :err
%FC% /c %FFLAGS% Linpack\dsort.f || goto :err
%FC% /c %FFLAGS% Linpack\dswap.f || goto :err
%FC% /c %FFLAGS% Linpack\dtrsv.f || goto :err
%FC% /c %FFLAGS% Linpack\fdump.f || goto :err
%FC% /c %FFLAGS% Linpack\i1mach.f || goto :err
%FC% /c %FFLAGS% Linpack\idamax.f || goto :err
%FC% /c %FFLAGS% Linpack\interv.f || goto :err
%FC% /c %FFLAGS% Linpack\j4save.f || goto :err
%FC% /c %FFLAGS% Linpack\lsame.f || goto :err
%FC% /c %FFLAGS% Linpack\ppvalu.f || goto :err
%FC% /c %FFLAGS% Linpack\xerbla.f || goto :err
%FC% /c %FFLAGS% Linpack\xercnt.f || goto :err
%FC% /c %FFLAGS% Linpack\xerhlt.f || goto :err
%FC% /c %FFLAGS% Linpack\xermsg.f || goto :err
%FC% /c %FFLAGS% Linpack\xerprn.f || goto :err
%FC% /c %FFLAGS% Linpack\xersve.f || goto :err
%FC% /c %FFLAGS% Linpack\xgetua.f || goto :err

rem ===== Compile (IMSL) =====
%FC% /c %FFLAGS% IMSL\c1dim.f || goto :err
%FC% /c %FFLAGS% IMSL\c1iarg.f || goto :err
%FC% /c %FFLAGS% IMSL\c1ind.f || goto :err
%FC% /c %FFLAGS% IMSL\c1tci.f || goto :err
%FC% /c %FFLAGS% IMSL\c1tic.f || goto :err
%FC% /c %FFLAGS% IMSL\c12ile.f || goto :err
%FC% /c %FFLAGS% IMSL\da1ot.f || goto :err
%FC% /c %FFLAGS% IMSL\dc1r.f || goto :err
%FC% /c %FFLAGS% IMSL\dc1trg.f || goto :err
%FC% /c %FFLAGS% IMSL\dc1wfr.f || goto :err
%FC% /c %FFLAGS% IMSL\dcsfrg.f || goto :err
%FC% /c %FFLAGS% IMSL\df2lsq.f || goto :err
%FC% /c %FFLAGS% IMSL\dfnlsq.f || goto :err
%FC% /c %FFLAGS% IMSL\dgirts.f || goto :err
%FC% /c %FFLAGS% IMSL\difnan.f || goto :err
%FC% /c %FFLAGS% IMSL\dmach.f || goto :err
%FC% /c %FFLAGS% IMSL\dr2ivn.f || goto :err
%FC% /c %FFLAGS% IMSL\dr3ivn.f || goto :err
%FC% /c %FFLAGS% IMSL\dset.f || goto :err
%FC% /c %FFLAGS% IMSL\dxyz.f || goto :err
%FC% /c %FFLAGS% IMSL\e1init.f || goto :err
%FC% /c %FFLAGS% IMSL\e1inpl.f || goto :err
%FC% /c %FFLAGS% IMSL\e1mes.f || goto :err
%FC% /c %FFLAGS% IMSL\e1pop.f || goto :err
%FC% /c %FFLAGS% IMSL\e1pos.f || goto :err
%FC% /c %FFLAGS% IMSL\e1prt.f || goto :err
%FC% /c %FFLAGS% IMSL\e1psh.f || goto :err
%FC% /c %FFLAGS% IMSL\e1std.f || goto :err
%FC% /c %FFLAGS% IMSL\e1sti.f || goto :err
%FC% /c %FFLAGS% IMSL\e1stl.f || goto :err
%FC% /c %FFLAGS% IMSL\e1ucs.f || goto :err
%FC% /c %FFLAGS% IMSL\e1usr.f || goto :err
%FC% /c %FFLAGS% IMSL\e3prt.f || goto :err
%FC% /c %FFLAGS% IMSL\i1cstr.f || goto :err
%FC% /c %FFLAGS% IMSL\i1dx.f || goto :err
%FC% /c %FFLAGS% IMSL\i1erif.f || goto :err
%FC% /c %FFLAGS% IMSL\i1kgt.f || goto :err
%FC% /c %FFLAGS% IMSL\i1kqu.f || goto :err
%FC% /c %FFLAGS% IMSL\i1krl.f || goto :err
%FC% /c %FFLAGS% IMSL\i1kst.f || goto :err
%FC% /c %FFLAGS% IMSL\i1x.f || goto :err
%FC% /c %FFLAGS% IMSL\iachar.f || goto :err
%FC% /c %FFLAGS% IMSL\icase.f || goto :err
%FC% /c %FFLAGS% IMSL\idanan.f || goto :err
%FC% /c %FFLAGS% IMSL\imach.f || goto :err
%FC% /c %FFLAGS% IMSL\iwkin.f || goto :err
%FC% /c %FFLAGS% IMSL\m1ve.f || goto :err
%FC% /c %FFLAGS% IMSL\m1vech.f || goto :err
%FC% /c %FFLAGS% IMSL\n1rcd.f || goto :err
%FC% /c %FFLAGS% IMSL\n1rgb.f || goto :err
%FC% /c %FFLAGS% IMSL\n1rty.f || goto :err
%FC% /c %FFLAGS% IMSL\s1anum.f || goto :err
%FC% /c %FFLAGS% IMSL\umach.f || goto :err

rem ===== Link to EXE =====
%FC% /nologo /Fe:ProPan2025_v1.0.exe *.obj /link /STACK:16777216,65536 || goto :err

rem ===== Stage outputs =====
move /y *.obj Code >nul
move /y *.mod Code >nul 2>nul
move /y *_genmod.f90 Code >nul 2>nul
if exist ProPan2025_v1.0.exe copy /y ProPan2025_v1.0.exe Code >nul

echo.
echo Build succeeded: ProPan2025_v1.0.exe
exit /b 0

:err
echo.
echo Build failed with errorlevel %errorlevel%.
exit /b 1

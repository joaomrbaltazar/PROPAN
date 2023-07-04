#!/bin/bash
# clean
mkdir -p Code
rm -f Code/*.*
rm -f ProPan2021_v1.0.out
# Source Folder
ifort -c -fltconsistency -fpe0 Base/propan_mod.f90
ifort -c -fltconsistency -fpe0 Base/ProPan2020_v1.2.f90
ifort -c -fltconsistency -fpe0 Base/progress.f90
ifort -c -fltconsistency -fpe0 Base/input.f90
ifort -c -fltconsistency -fpe0 Base/inivars.f90
ifort -c -fltconsistency -fpe0 Base/delvars.f90
ifort -c -fltconsistency -fpe0 Base/counters.f90
ifort -c -fltconsistency -fpe0 Base/output.f90
# Geom Folder
ifort -c -fltconsistency -fpe0 Geom/bladegeom.f90
ifort -c -fltconsistency -fpe0 Geom/bladewakegeom.f90
ifort -c -fltconsistency -fpe0 Geom/nozzlegeom.f90
ifort -c -fltconsistency -fpe0 Geom/nozzlewakegeom.f90
ifort -c -fltconsistency -fpe0 Geom/hubgeom.f90
ifort -c -fltconsistency -fpe0 Geom/panel.f90
ifort -c -fltconsistency -fpe0 Geom/pancoor.f90
# Numdif Folder
ifort -c -fltconsistency -fpe0 Numdif/numdifblade.f90
ifort -c -fltconsistency -fpe0 Numdif/numdifbladewake.f90
ifort -c -fltconsistency -fpe0 Numdif/numdifnozzle.f90
ifort -c -fltconsistency -fpe0 Numdif/numdifhub.f90
# InfCoef Folder
ifort -c -fltconsistency -fpe0 InfCoef/InfCoefMatrixStd.f90
ifort -c -fltconsistency -fpe0 InfCoef/InfCoefMatrixUnStd.f90
# PotCoef Folder
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/bladecoef.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/bladewakecoef.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/nozzlecoef.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/nozzlewakecoef.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/hubcoef.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/potpan.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/potpan_num.f90
ifort -c -fltconsistency -fpe0 InfCoef/PotCoef/potlinsub.f90
# VelCoef Folder
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/bladevelo.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/bladewakevelo.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/nozzlevelo.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/nozzlewakevelo.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/hubvelo.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/imagehubvelo.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/velpan.f90
ifort -c -fltconsistency -fpe0 InfCoef/VelCoef/gaussvel.f90
# Solve Folder
ifort -c -fltconsistency -fpe0 Solve/SolveRhsWet.f90
ifort -c -fltconsistency -fpe0 Solve/SolveRhsCav.f90
ifort -c -fltconsistency -fpe0 Solve/SolveRhsCavRed.f90
ifort -c -fltconsistency -fpe0 Solve/SolveLkWet.f90
ifort -c -fltconsistency -fpe0 Solve/SolveLkCav.f90
ifort -c -fltconsistency -fpe0 Solve/SolveIpkcWet.f90
ifort -c -fltconsistency -fpe0 Solve/SolveIpkcCav.f90
ifort -c -fltconsistency -fpe0 Solve/SolveWake.f90
# WakeAlign Folder
ifort -c -fltconsistency -fpe0 WakeAlign/wakealign1.f90
ifort -c -fltconsistency -fpe0 WakeAlign/nozzledef.f90
ifort -c -fltconsistency -fpe0 WakeAlign/geoduct37.f90
ifort -c -fltconsistency -fpe0 WakeAlign/bladewakedisp.f90
ifort -c -fltconsistency -fpe0 WakeAlign/nozzledisp.f90
ifort -c -fltconsistency -fpe0 WakeAlign/nozzlewakedisp.f90
ifort -c -fltconsistency -fpe0 WakeAlign/ff.f90
cp WakeAlign/wakealign2.o .
# Cav Folder
ifort -c -fltconsistency -fpe0 Cav/cavcheck.f90
ifort -c -fltconsistency -fpe0 Cav/cavprop.f90
ifort -c -fltconsistency -fpe0 Cav/caverrc.f90
ifort -c -fltconsistency -fpe0 Cav/cavpotp.f90
ifort -c -fltconsistency -fpe0 Cav/cavpots.f90
ifort -c -fltconsistency -fpe0 Cav/cavthickp.f90
ifort -c -fltconsistency -fpe0 Cav/cavthicks.f90
ifort -c -fltconsistency -fpe0 Cav/cavrecover.f90
# VtCp Folder
ifort -c -fltconsistency -fpe0 VtCp/velp.f90
ifort -c -fltconsistency -fpe0 VtCp/velpw.f90
ifort -c -fltconsistency -fpe0 VtCp/presp.f90
ifort -c -fltconsistency -fpe0 VtCp/veln.f90
ifort -c -fltconsistency -fpe0 VtCp/presn.f90
ifort -c -fltconsistency -fpe0 VtCp/presnte.f90
ifort -c -fltconsistency -fpe0 VtCp/velh.f90
ifort -c -fltconsistency -fpe0 VtCp/presh.f90
# Wake Folder
ifort -c -fltconsistency -fpe0 Calc/Wake/vwake.f90
ifort -c -fltconsistency -fpe0 Calc/Wake/fourier_coef.f90
ifort -c -fltconsistency -fpe0 Calc/Wake/fourier_function.f90
# Grids Folder
ifort -c -fltconsistency -fpe0 Grids/hubgrid.f90
ifort -c -fltconsistency -fpe0 Grids/nozzlegrid.f90
ifort -c -fltconsistency -fpe0 Grids/nozzlewakegrid.f90
# Grape Folder
ifort -c -fltconsistency -fpe0 Grids/Grape/angri.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/bord.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/calcb.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/calphi.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/coef.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/grape.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/guessa.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/rhs.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/sip.f90
ifort -c -fltconsistency -fpe0 Grids/Grape/splin.f90
# Field Folder
ifort -c -fltconsistency -fpe0 Field/bladecoeff.f90
ifort -c -fltconsistency -fpe0 Field/bladewakecoeff.f90
ifort -c -fltconsistency -fpe0 Field/hubcoeff.f90
ifort -c -fltconsistency -fpe0 Field/nozzlecoeff.f90
ifort -c -fltconsistency -fpe0 Field/nozzlewakecoeff.f90
ifort -c -fltconsistency -fpe0 Field/velfStd.f90
ifort -c -fltconsistency -fpe0 Field/velfUnStd.f90
ifort -c -fltconsistency -fpe0 Field/presfStd.f90
ifort -c -fltconsistency -fpe0 Field/presfUnStd.f90
# Calc Folder
ifort -c -fltconsistency -fpe0 Calc/JacobianWet.f90
ifort -c -fltconsistency -fpe0 Calc/JacobianCav.f90
ifort -c -fltconsistency -fpe0 Calc/linint.f90
ifort -c -fltconsistency -fpe0 Calc/intk1.f90
ifort -c -fltconsistency -fpe0 Calc/splint.f90
ifort -c -fltconsistency -fpe0 Calc/spline.f90
ifort -c -fltconsistency -fpe0 Calc/ispline.f90
ifort -c -fltconsistency -fpe0 Calc/gaussint.f90
ifort -c -fltconsistency -fpe0 Calc/bisof.f
ifort -c -fltconsistency -fpe0 Calc/CtCq.f90
ifort -c -fltconsistency -fpe0 Calc/gaperrg.f90
ifort -c -fltconsistency -fpe0 Calc/periodicflow.f90
ifort -c -fltconsistency -fpe0 Calc/velinf.f90
ifort -c -fltconsistency -fpe0 Calc/stret2.f90
ifort -c -fltconsistency -fpe0 Calc/sxx.f90
ifort -c -fltconsistency -fpe0 Calc/shxx.f90
ifort -c -fltconsistency -fpe0 Calc/sdiv.f90
ifort -c -fltconsistency -fpe0 Calc/frame.f90
# Linpack Folder
ifort -c -fltconsistency -fpe0 Linpack/cubspl.f
ifort -c -fltconsistency -fpe0 Linpack/daxpy.f
ifort -c -fltconsistency -fpe0 Linpack/dcopy.f
ifort -c -fltconsistency -fpe0 Linpack/ddot.f
ifort -c -fltconsistency -fpe0 Linpack/dgedi.f
ifort -c -fltconsistency -fpe0 Linpack/dgefa.f
ifort -c -fltconsistency -fpe0 Linpack/dgemm.f
ifort -c -fltconsistency -fpe0 Linpack/dger.f
ifort -c -fltconsistency -fpe0 Linpack/dgesl.f
ifort -c -fltconsistency -fpe0 Linpack/drotm.f
ifort -c -fltconsistency -fpe0 Linpack/drotmg.f
ifort -c -fltconsistency -fpe0 Linpack/dscal.f
ifort -c -fltconsistency -fpe0 Linpack/dsort.f
ifort -c -fltconsistency -fpe0 Linpack/dswap.f
ifort -c -fltconsistency -fpe0 Linpack/dtrsv.f
ifort -c -fltconsistency -fpe0 Linpack/fdump.f
ifort -c -fltconsistency -fpe0 Linpack/i1mach.f
ifort -c -fltconsistency -fpe0 Linpack/idamax.f
ifort -c -fltconsistency -fpe0 Linpack/interv.f
ifort -c -fltconsistency -fpe0 Linpack/j4save.f
ifort -c -fltconsistency -fpe0 Linpack/lsame.f
ifort -c -fltconsistency -fpe0 Linpack/ppvalu.f
ifort -c -fltconsistency -fpe0 Linpack/xerbla.f
ifort -c -fltconsistency -fpe0 Linpack/xercnt.f
ifort -c -fltconsistency -fpe0 Linpack/xerhlt.f
ifort -c -fltconsistency -fpe0 Linpack/xermsg.f
ifort -c -fltconsistency -fpe0 Linpack/xerprn.f
ifort -c -fltconsistency -fpe0 Linpack/xersve.f
ifort -c -fltconsistency -fpe0 Linpack/xgetua.f
# IMSL Folder
ifort -c -fltconsistency -fpe0 IMSL/c1dim.f
ifort -c -fltconsistency -fpe0 IMSL/c1iarg.f
ifort -c -fltconsistency -fpe0 IMSL/c1ind.f
ifort -c -fltconsistency -fpe0 IMSL/c1tci.f
ifort -c -fltconsistency -fpe0 IMSL/c1tic.f
ifort -c -fltconsistency -fpe0 IMSL/c12ile.f
ifort -c -fltconsistency -fpe0 IMSL/da1ot.f
ifort -c -fltconsistency -fpe0 IMSL/dc1r.f
ifort -c -fltconsistency -fpe0 IMSL/dc1trg.f
ifort -c -fltconsistency -fpe0 IMSL/dc1wfr.f
ifort -c -fltconsistency -fpe0 IMSL/dcsfrg.f
ifort -c -fltconsistency -fpe0 IMSL/df2lsq.f
ifort -c -fltconsistency -fpe0 IMSL/dfnlsq.f
ifort -c -fltconsistency -fpe0 IMSL/dgirts.f
ifort -c -fltconsistency -fpe0 IMSL/difnan.f
ifort -c -fltconsistency -fpe0 IMSL/dmach.f
ifort -c -fltconsistency -fpe0 IMSL/dr2ivn.f
ifort -c -fltconsistency -fpe0 IMSL/dr3ivn.f
ifort -c -fltconsistency -fpe0 IMSL/dset.f
ifort -c -fltconsistency -fpe0 IMSL/dxyz.f
ifort -c -fltconsistency -fpe0 IMSL/e1init.f
ifort -c -fltconsistency -fpe0 IMSL/e1inpl.f
ifort -c -fltconsistency -fpe0 IMSL/e1mes.f
ifort -c -fltconsistency -fpe0 IMSL/e1pop.f
ifort -c -fltconsistency -fpe0 IMSL/e1pos.f
ifort -c -fltconsistency -fpe0 IMSL/e1prt.f
ifort -c -fltconsistency -fpe0 IMSL/e1psh.f
ifort -c -fltconsistency -fpe0 IMSL/e1std.f
ifort -c -fltconsistency -fpe0 IMSL/e1sti.f
ifort -c -fltconsistency -fpe0 IMSL/e1stl.f
ifort -c -fltconsistency -fpe0 IMSL/e1ucs.f
ifort -c -fltconsistency -fpe0 IMSL/e1usr.f
ifort -c -fltconsistency -fpe0 IMSL/e3prt.f
ifort -c -fltconsistency -fpe0 IMSL/i1cstr.f
ifort -c -fltconsistency -fpe0 IMSL/i1dx.f
ifort -c -fltconsistency -fpe0 IMSL/i1erif.f
ifort -c -fltconsistency -fpe0 IMSL/i1kgt.f
ifort -c -fltconsistency -fpe0 IMSL/i1kqu.f
ifort -c -fltconsistency -fpe0 IMSL/i1krl.f
ifort -c -fltconsistency -fpe0 IMSL/i1kst.f
ifort -c -fltconsistency -fpe0 IMSL/i1x.f
ifort -c -fltconsistency -fpe0 IMSL/iachar.f
ifort -c -fltconsistency -fpe0 IMSL/icase.f
ifort -c -fltconsistency -fpe0 IMSL/idanan.f
ifort -c -fltconsistency -fpe0 IMSL/imach.f
ifort -c -fltconsistency -fpe0 IMSL/iwkin.f
ifort -c -fltconsistency -fpe0 IMSL/m1ve.f
ifort -c -fltconsistency -fpe0 IMSL/m1vech.f
ifort -c -fltconsistency -fpe0 IMSL/n1rcd.f
ifort -c -fltconsistency -fpe0 IMSL/n1rgb.f
ifort -c -fltconsistency -fpe0 IMSL/n1rty.f
ifort -c -fltconsistency -fpe0 IMSL/s1anum.f
ifort -c -fltconsistency -fpe0 IMSL/umach.f
# Executable
ifort -o ProPan2021_v1.0.out *.o
mv *.o Code
mv *.mod Code
if [ -f ProPan2021_v1.0.out ]; then
   cp ProPan2021_v1.0.out Code
fi

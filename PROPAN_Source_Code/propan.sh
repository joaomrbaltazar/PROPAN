#!/bin/bash
# clean
rm -f Code/*.*
##rm -f ProPan2020_v1.2.out
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
##ifort -c -fltconsistency -fpe0 Source/Geom/bladegeom.f90
##ifort -c -fltconsistency -fpe0 Source/Geom/bladewakegeom.f90
##ifort -c -fltconsistency -fpe0 Source/Geom/nozzlegeom.f90
##ifort -c -fltconsistency -fpe0 Source/Geom/nozzlewakegeom.f90
##ifort -c -fltconsistency -fpe0 Source/Geom/hubgeom.f90
##ifort -c -fltconsistency -fpe0 Source/Geom/panel.f90
##ifort -c -fltconsistency -fpe0 Source/Geom/pancoor.f90
# Numdif Folder
##ifort -c -fltconsistency -fpe0 Source/Numdif/numdifblade.f90
##ifort -c -fltconsistency -fpe0 Source/Numdif/numdifbladewake.f90
##ifort -c -fltconsistency -fpe0 Source/Numdif/numdifnozzle.f90
##ifort -c -fltconsistency -fpe0 Source/Numdif/numdifhub.f90
# InfCoef Folder
##ifort -c -fltconsistency -fpe0 Source/InfCoef/InfCoefMatrixStd.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/InfCoefMatrixUnStd.f90
# PotCoef Folder
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/bladecoef.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/bladewakecoef.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/nozzlecoef.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/nozzlewakecoef.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/hubcoef.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/potpan.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/potpan_num.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/PotCoef/potlinsub.f90
# VelCoef Folder
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/bladevelo.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/bladewakevelo.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/nozzlevelo.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/nozzlewakevelo.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/hubvelo.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/imagehubvelo.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/velpan.f90
##ifort -c -fltconsistency -fpe0 Source/InfCoef/VelCoef/gaussvel.f90
# Solve Folder
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveRhsWet.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveRhsCav.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveRhsCavRed.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveLkWet.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveLkCav.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveIpkcWet.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveIpkcCav.f90
##ifort -c -fltconsistency -fpe0 Source/Solve/SolveWake.f90
# WakeAlign Folder
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/wakealign1.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/wakealign2.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/nozzledef.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/geoduct37.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/bladewakedisp.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/nozzledisp.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/nozzlewakedisp.f90
##ifort -c -fltconsistency -fpe0 Source/WakeAlign/ff.f90
# Cav Folder
##ifort -c -fltconsistency -fpe0 Source/Cav/cavcheck.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/cavprop.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/caverrc.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/cavpotp.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/cavpots.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/cavthickp.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/cavthicks.f90
##ifort -c -fltconsistency -fpe0 Source/Cav/cavrecover.f90
# VtCp Folder
##ifort -c -fltconsistency -fpe0 Source/VtCp/velp.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/velpw.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/presp.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/veln.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/presn.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/presnte.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/velh.f90
##ifort -c -fltconsistency -fpe0 Source/VtCp/presh.f90
# Wake Folder
##ifort -c -fltconsistency -fpe0 Source/Calc/Wake/vwake.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/Wake/fourier_coef.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/Wake/fourier_function.f90
# Grids Folder
##ifort -c -fltconsistency -fpe0 Source/Grids/hubgrid.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/nozzlegrid.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/nozzlewakegrid.f90
# Grape Folder
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/angri.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/bord.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/calcb.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/calphi.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/coef.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/grape.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/guessa.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/rhs.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/sip.f90
##ifort -c -fltconsistency -fpe0 Source/Grids/Grape/splin.f90
# Field Folder
##ifort -c -fltconsistency -fpe0 Source/Field/bladecoeff.f90
##ifort -c -fltconsistency -fpe0 Source/Field/bladewakecoeff.f90
##ifort -c -fltconsistency -fpe0 Source/Field/hubcoeff.f90
##ifort -c -fltconsistency -fpe0 Source/Field/nozzlecoeff.f90
##ifort -c -fltconsistency -fpe0 Source/Field/nozzlewakecoeff.f90
##ifort -c -fltconsistency -fpe0 Source/Field/velfStd.f90
##ifort -c -fltconsistency -fpe0 Source/Field/velfUnStd.f90
##ifort -c -fltconsistency -fpe0 Source/Field/presfStd.f90
##ifort -c -fltconsistency -fpe0 Source/Field/presfUnStd.f90
# Calc Folder
##ifort -c -fltconsistency -fpe0 Source/Calc/JacobianWet.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/JacobianCav.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/linint.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/intk1.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/splint.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/spline.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/ispline.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/gaussint.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/bisof.f
##ifort -c -fltconsistency -fpe0 Source/Calc/CtCq.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/gaperrg.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/periodicflow.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/velinf.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/stret2.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/sxx.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/shxx.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/sdiv.f90
##ifort -c -fltconsistency -fpe0 Source/Calc/frame.f90
# Linpack Folder
##ifort -c -fltconsistency -fpe0 Source/Linpack/cubspl.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/daxpy.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dcopy.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/ddot.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dgedi.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dgefa.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dgemm.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dger.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dgesl.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/drotm.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/drotmg.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dscal.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dsort.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dswap.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/dtrsv.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/fdump.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/i1mach.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/idamax.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/interv.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/j4save.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/lsame.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/ppvalu.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xerbla.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xercnt.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xerhlt.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xermsg.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xerprn.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xersve.f
##ifort -c -fltconsistency -fpe0 Source/Linpack/xgetua.f
# IMSL Folder
##ifort -c -fltconsistency -fpe0 Source/IMSL/c1dim.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/c1iarg.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/c1ind.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/c1tci.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/c1tic.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/c12ile.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/da1ot.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dc1r.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dc1trg.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dc1wfr.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dcsfrg.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/df2lsq.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dfnlsq.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dgirts.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/difnan.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dmach.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dr2ivn.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dr3ivn.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dset.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/dxyz.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1init.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1inpl.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1mes.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1pop.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1pos.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1prt.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1psh.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1std.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1sti.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1stl.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1ucs.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e1usr.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/e3prt.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1cstr.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1dx.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1erif.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1kgt.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1kqu.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1krl.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1kst.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/i1x.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/iachar.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/icase.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/idanan.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/imach.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/iwkin.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/m1ve.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/m1vech.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/n1rcd.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/n1rgb.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/n1rty.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/s1anum.f
##ifort -c -fltconsistency -fpe0 Source/IMSL/umach.f
# Executable
##ifort -o ProPan2020_v1.2.out *.o
##mv *.o Code
##mv *.mod Code
##if [ -f ProPan2020_v1.2.out ]; then
##   cp ProPan2020_v1.2.out Code
##fi

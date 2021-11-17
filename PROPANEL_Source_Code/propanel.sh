#!/bin/bash
# clean
mkdir -p Code
rm -f Code/*.*
rm -f ProPanel2021_v1.0_debug.out
# Source Folder
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/propanel_mod.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/ProPanel2021_v1.0.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/delvars.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/progress.f90
# Grids Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/bladegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/bladewakegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/nozzlegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/nozzlewakegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/nozzledef.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/geoduct37.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/hubgrid.f90
# Grape Folder
cp Grape/angri.o .
cp Grape/bord.o .
cp Grape/calcb.o .
cp Grape/calphi.o .
cp Grape/coef.o .
cp Grape/grape.o .
cp Grape/guessa.o .
cp Grape/rhs.o .
cp Grape/sip.o .
cp Grape/splin.o .
# Calc Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/linint.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/intk1.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/splint.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/spline.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/ispline.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/stret_choice.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/dscal.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/stret.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/stret2.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/sxx.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/shxx.f90
# Linpack Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Linpack/cubspl.f
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Linpack/ppvalu.f
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Linpack/interv.f
# Executable
ifort -o ProPanel2021_v1.0_debug.out *.o
mv *.o Code
mv *.mod Code
if [ -f ProPanel2021_v1.0_debug.out ]; then
   cp ProPanel2021_v1.0_debug.out Code
fi
#!/bin/bash
# clean
mkdir -p Code
rm -f Code/*.*
rm -f ProPost2021_v1.0_debug.out
# Source Folder
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/propost_mod.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/ProPost2020_v1.0.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/progress.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/input.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Base/delvars.f90
# Plots Folder
ifort -c -check bounds -traceback -fltconsistency -fpe0 Plots/plotcp.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Plots/plotharm.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Plots/plotwake.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Plots/plotsol.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Plots/plotcav.f90
# Calc Folder
ifort -c -check bounds -traceback -fltconsistency -fpe0 Calc/linint.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Calc/intk1.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Calc/splint.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Calc/periodicflow.f90
cp Calc/spline.o .
cp Calc/ispline.o .
# Visc Folder
ifort -c -check bounds -traceback -fltconsistency -fpe0 Visc/viscor.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Visc/clcd.f90
ifort -c -check bounds -traceback -fltconsistency -fpe0 Visc/interval.f90
# Geom Folder
ifort -c -check bounds -traceback -fltconsistency -fpe0 Geom/panel.f90
# Executable
ifort -o ProPost2021_v1.0_debug.out *.o
mv *.o Code
mv *.mod Code
if [ -f ProPost2021_v1.0_debug.out ]; then
   cp ProPost2021_v1.0_debug.out Code
fi

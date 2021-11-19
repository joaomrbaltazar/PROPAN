#!/bin/bash
# clean
mkdir -p Code
rm -f Code/*.*
rm -f ProPost2021_v1.0_debug.out
# Source Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/propost_mod.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/ProPost2020_v1.0.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/progress.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/input.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/delvars.f90
# Plots Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Plots/plotcp.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Plots/plotharm.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Plots/plotwake.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Plots/plotsol.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Plots/plotcav.f90
# Calc Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/linint.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/intk1.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/splint.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/spline.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/ispline.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Calc/periodicflow.f90
# Viscor Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Visc/viscor.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Visc/clcd.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Visc/interval.f90
# Geom Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Geom/panel.f90
# Executable
ifort -o ProPost2021_v1.0_debug.out *.o
mv *.o Code
mv *.mod Code
if [ -f ProPost2021_v1.0_debug.out ]; then
   cp ProPost2021_v1.0_debug.out Code
fi

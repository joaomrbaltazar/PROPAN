#!/bin/bash
# clean
##rm -f Code/*.*
##rm -f ProPanel2020_v1.1_debug.out
# Source Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/propanel_mod.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/ProPanel2020_v1.1.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/delvars.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/progress.f90
# Grids Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/bladegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/bladewakegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/nozzlegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/nozzlewakegrid.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/nozzledef.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/geoduct37.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grids/hubgrid.f90
# Grape Folder
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/angri.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/bord.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/calcb.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/calphi.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/coef.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/grape.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/guessa.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/rhs.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/sip.f90
##ifort -c -check bounds -traceback -fltconsistency -fpe0 Source/Grape/splin.f90
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
##ifort -o ProPanel2020_v1.1_debug.out *.o
##mv *.o Code
##mv *.mod Code
##if [ -f ProPanel2020_v1.1_debug.out ]; then
##   cp ProPanel2020_v1.1_debug.out Code
##fi

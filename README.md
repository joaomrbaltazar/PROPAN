# PROPAN - Potential Flow Code for Foils and Rotors
PROPAN is short for Propeller Panel Method. PROPAN is a panel code for the calculation of steady and unsteady potential flow around foils, open and ducted propellers, and wind and marine current turbines. PROPAN was developed by MARETEC (Marine and Environmental Technology Research Centre) at Instituto Superior Técnico (IST) which belongs to Lisbon University.
# What is this repository?
This is the PROPAN Potential Flow Code OFFICIAL repository.
# Overview
PROPAN Panel Code is complemented by a pre-processor (PROPANEL) and a post-processor (PROPOST). PROPANEL is a 3D surface grid generation tool. PROPOST is a post-processing tool that generates several output and solution files.
# PROPANEL capabilities:
Input: ASCII text files

Wings:
* conventional and quasi-orthogonal grids
* conventional grid for delta wings
* trailing-edge wake grid (user-specify or empirical models)
* leading-edge wake grid (empirical model for delta wings)

Propellers & Turbines: (Open or Ducted)
* conventional and quasi-orthogonal blade grids
* hub grid
* duct grid
* blade and duct wake grids (user-specify or empirical models)

Output: grid output in ASCII text files with Tecplot format
# PROPAN capabilities:
Input: ASCII text files

Wings: (not tested in latest version)
* steady flow
* empirical rigid wake model
* wetted flow
* leading-edge flow separation for delta wings

Marine Propellers & Turbines: (Open or Ducted)
* steady and unsteady flows
* empirical rigid wake model and wake alignment model
* wetted Flow, sheet partial (both face and back sides) and super cavitation (one side only) on the blades
* leading-edge flow separation (not tested in latest version)

Wind Turbines:
* steady and unsteady flows
* empirical rigid wake model
* dynamic inflow

Output: solution output in ASCII text files with Tecplot format

# PROPOST capabilities:
Input: ASCII text files

Output: solution output in ASCII text files with Tecplot format
* 2D blade pressure distribution
* 2D duct pressure distribution
* wake geometry at downstream plane
* 2D cavitation patterns
* 2D viscous effects (section lift and drag corrections)
* harmonic analysis for Unsteady Flow

PROPANEL, PROPAN and PROPOST run on Windows and Linux workstations. All routines are written in FORTRAN 95, combined with LINPACK and IMSL FORTRAN 77 routines. The code is not parallelised.
# Help, Bugs, Feedback
If you need help with PROPAN, chat with developers or ask any other questions about PROPAN, you can hang out by email: propan.code@gmail.com. To report bugs, please create a GitHub issue or contact by email. More information consult: https://www.researchgate.net/project/PROPAN-potential-flow-code-for-foils-and-rotors/
# License
GNU Affero General Public License v3.0

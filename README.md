# PROPAN - Potential Flow Code for Foils and Rotors
PROPAN is short for Propeller Panel Method. PROPAN is a panel code for the calculation of steady and unsteady potential flow around foils, open and ducted propellers, and wind and marine current turbines. PROPAN was developed by MARETEC (Marine and Environmental Technology Research Centre) at Instituto Superior Técnico (IST) which belongs to Lisbon University.
# What is this repository?
This is the PROPAN Potential Flow Code OFFICIAL repository.
# Overview
Panel code PROPAN has the following key capabilities:

Input: ACII text files

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

PROPAN runs on Windows and Linux workstations. All routines are written in FORTRAN 95, combined with LINPACK and IMSL FORTRAN 77 routines. The code is not parallelised.
# Help, Bugs, Feedback
If you need help with PROPAN, chat with developers or ask any other questions about PROPAN, you can hang out by email: propan.code@gmail.com. To report bugs, please create a GitHub issue or contact by email. More information consult: https://www.researchgate.net/project/PROPAN-potential-flow-code-for-foils-and-rotors/
# License
GNU Affero General Public License v3.0

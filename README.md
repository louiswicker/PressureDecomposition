# README for p_decomp program

> written by Lou Wicker, October 2024

# Contents

1. GENERAL DESCRIPTION
2. IMPORTANT NOTES
3. HOW TO RUN
#

1. GENERAL DESCRIPTION

- This program was originally created by George Bryan (NCAR) and modified by Jeff Trapp (UIUC) and Geoff Marion (CIWRO). Similar code is implemented
in CM1, but it is not set up for runs using distributed memory (the majority of CM1 simulations). 

- The code has been reoganized by Lou Wicker starting in summer of 2024. Closer attention has been paid to the buoyant pressure retrieval code. Over the years and through various users, incorrect vertical boundary conditions crept into the specification of the vertical tridiagonal weights.

- The buoyant retrieval has been carefully compared to the methods and tests described in:

> Jeevanjee, N., and D. M. Romps, 2015: Effective Buoyancy, Inertial Pressure, and the Mechanical Generation of Boundary Layer Mass Flux by Cold Pools. J. Atmos. Sci., 72, 3199–3213, https://doi.org/10.1175/JAS-D-14-0349.1.

 &nbsp and the beta forcing tests codes use the tests from that paper.

- Their approach uses the non-base state approach first introduced by:

> Davies-Jones, R., 2003: An expression for effective buoyancy in surroundings with horizontal density gradients. J. Atmos. Sci., 60, 2922–2925, https://doi.org/10.1175/1520-0469(2003)060&#60;2922:aefebi&#62;2.0.co;2.

  and discussed further and used in 3D storm analysis by:

> Dawson, D. T., M. Xue, A. Shapiro, J. A. Milbrandt, and A. D. Schenkman, 2016: Sensitivity of Real-Data Simulations of the 3 May 1999 Oklahoma City Tornadic Supercell and Associated Tornadoes to Multimoment Microphysics. Part II: Analysis of Buoyancy and Dynamic Pressure Forces in Simulated Tornado-Like Vortices. J. Atmos. Sci., 73, 1039–1061, https://doi.org/10.1175/JAS-D-15-0114.1.

#

2. IMPORTANT NOTES
- Program is currently set up to read/write netcdf output. Original
  program only read/write grads.
- DO NOT RUN PDCOMP OVER A SUBDOMAIN IN CM1. pdcomp uses information
  taken from horizontally averaging some variables at the top of the
  analysis domain. If you run pdcomp for two subdomains of the same
  model run, they WILL NOT be comparable. It is especially bad to do this
  with the pdcomp domain not extending the full depth of the model domain.
  Just don't do it (unless you're prepared to make very extensive
  modifications of the pdcomp calculations, CM1 output reading subroutines,
  etc.)!
- If you want to run pdcomp for a model run that uses an unsupported
  (i.e., not Morrison or NSSL double-moment) microphysics scheme, it's
  relatively easy to add that functionality. Simply change the if-statements
  that check if the input ptype is supported and add an additional
  if-statement where qtot is calculated (for calculating the buoyancy
  pressure) to use whatever mixing ratios your scheme uses.
- If you are using pdcomp on a model run that includes boundary layer
  turbulence from random potential temperature perturbations, recomputing
  the base state variables (e.g., th0, p0) may be necessary in order to
  account for modifications to the base state by said turbulence. There is a
  loop in getpp.f90 (immediately following 'Checkpoint 1') that is commented
  out by default that does this by averaging the base state values over some
  subdomain (hard set variables avgstart,avglen). It's recommended to do the
  averaging over as large a subdomain as possible to ensure that it is
  representative of this new base state.

#

3. BEFORE YOU RUN:
- Modify range of desired CM1 output times to perform analysis in run_pdcomp.bash
- Modify the filepaths in the def.pdcomp.input file and run_pdcomp.bash script
  to reflect your pdcomp and CM1 output directories.
- If not already done, make sure to set the LD_LIBRARY_PATH environment
  variable within your .bashrc (or whatever is relevant for the shell
  your machine is using).
- In the run directory, create a symbolic link to the netcdf.mod module file
  on your machine.


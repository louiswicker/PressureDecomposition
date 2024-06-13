#---------------------------#
README for pdcomp program #
# 
Original README (needs updating)
#---------------------------#

Contents:


How to compile PDcomp 


It might be easiest to install the hpc conda environment.  I use this environment for a number of projects,
one can compile CM1 using it as Openmpi is also available.  The file "gcc.env" in this directory is included in the src/Makefile to help specify the fortran commands (debug or otherwise) and the netCDF libs.

Installation Using Conda Procedure
----------------------------------

To install the hpc environment into your conda system, type:

"conda env create -f hpc.yml"

If that works, then go into the environment"

"conda activate hpc"

and you should now be able to compile (one more thing):

I assume you are using miniconda3 for your python environment - if so - and if your miniconda3 is located at

${HOME}/miniconda3

then the gcc.env will work.  

If not, but the "hpc" is installed, please edit the gcc.env to point to the conda install directory, 

${HOME}/your_conda_install

in the gcc.env file.

Installation Using Linux Sys
----------------------------

Edit the gcc.env file to where the include and libs are for netCDF4, e.g., edit the 

OUTPUTINC = 

LINKOPTS = 

macro definitions, as these are the ones used by the Makefile.  


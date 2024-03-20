#---------------------#
# Makefile for pdcomp #
#---------------------#
#  Config variables   #

# Fortran compiler
F90=gfortran

# Locations of netcdf include files and libraries

# Locations of FFTW include files and libraries
FFTWINC=-I/usr/include # location of include files for fftw
FFTWLIB=-L/usr/local/lib64 # location of libraries for fftw

include gcc.site

NCINC=${OUTPUTINC}
NCLIB=${LINKOPTS}

#################################################################################
# YOU PROBABLY DON'T WANT TO CHANGE ANYTHING BELOW THIS                         #
#################################################################################

OBFILES = getpp.o bc.o singleton.o pdcomp.o
OUTPUTOPT = -DNETCDF
LINKOPTS  = -lnetcdff

FOPTS = -O2 -fcheck=bounds
FOPTS = -g -O0 -fbacktrace -fcheck=all -Wall

#################################################################################

pdcomp :        $(OBFILES)
#		$(F90) -O3 -o ../run/pdcomp.exe $(OBFILES) -fcheck=bounds $(FFTWINC) $(FFTWLIB) -lfftw3 $(NCLIB) -lm $(LINKOPTS)
		$(F90) $(FOPTS) -o ../run/pdcomp.exe $(OBFILES) $(FFTWINC) $(FFTWLIB) $(NCLIB) -lm $(LINKOPTS)

getpp.o :       getpp.f90
		$(F90) -c $(FOPTS) getpp.f90 $(NCINC) $(NCLIB)

bc.o :          bc.f90
		$(F90) -c $(FOPTS) bc.f90 $(NCINC) $(NCLIB)

singleton.o :   singleton.f90
		$(F90) -c $(FOPTS) singleton.f90 $(NCINC) $(NCLIB)

pdcomp.o : 	pdcomp.f90
		$(F90) -c $(FOPTS) pdcomp.f90 $(NCINC) $(NCLIB)

clean:
		rm -rf *o *mod

#---------------------#
# Makefile for retrieving pressure from cloud models
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

OBFILES = bc.o singleton.o solve_ellip.o array_utils.o ncdf_utils.o
OUTPUTOPT = -DNETCDF
LINKOPTS  = -lnetcdff

FOPTS = -O2 -fcheck=bounds -Wno-unused-parameter 

FOPTS = -g -O0 -fbacktrace -fcheck=all -w -Wall -Wno-unused-parameter 

#################################################################################

retrievebeta:   $(OBFILES) retrievebeta.o
		$(F90) $(FOPTS) -o retrievebeta.exe retrievebeta.o $(OBFILES) $(FFTWINC) $(FFTWLIB) $(NCLIB) -lm $(LINKOPTS)

testbeta:   $(OBFILES) testbeta.o
		$(F90) $(FOPTS) -o testbeta.exe testbeta.o $(OBFILES) $(FFTWINC) $(FFTWLIB) $(NCLIB) -lm $(LINKOPTS)

testpb:     $(OBFILES) testpb.o
		$(F90) $(FOPTS) -o testpb.exe testpb.o $(OBFILES) $(FFTWINC) $(FFTWLIB) $(NCLIB) -lm $(LINKOPTS)

retrievebeta.o:   retrievebeta.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

testbeta.o:   testbeta.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

testpb.o:   testpb.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

array_utils.o: array_utils.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

ncdf_utils.o:  ncdf_utils.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

bc.o:       bc.f90
		$(F90) -c $(FOPTS) bc.f90 

singleton.o:   singleton.f90
		$(F90) -c $(FOPTS) singleton.f90 

solve_ellip.o:	solve_ellip.f90
		$(F90) -c $(FOPTS) solve_ellip.f90 

clean:
		rm -rf *.o *.mod *.exe

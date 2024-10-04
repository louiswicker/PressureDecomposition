#---------------------#
# Makefile for retrieving pressure from cloud models
#---------------------#
#  Config variables   #

VPATH=main:test:tools

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
		cp retrievebeta.exe ../run/retrievebeta.exe

testbeta:   $(OBFILES) testbeta.o
		$(F90) $(FOPTS) -o testbeta.exe testbeta.o $(OBFILES) $(FFTWINC) $(FFTWLIB) $(NCLIB) -lm $(LINKOPTS)

testpb:     $(OBFILES) testpb.o
		$(F90) $(FOPTS) -o testpb.exe testpb.o $(OBFILES) $(FFTWINC) $(FFTWLIB) $(NCLIB) -lm $(LINKOPTS)

retrievebeta.o:  main/retrievebeta.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

testbeta.o: test/testbeta.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

testpb.o:   test/testpb.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

array_utils.o: tools/array_utils.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

ncdf_utils.o:  tools/ncdf_utils.f90
		$(F90) -c $(FOPTS) $< $(NCINC) $(NCLIB)

bc.o:       tools/bc.f90
		$(F90) -c $(FOPTS) $<

singleton.o:  main/singleton.f90
		$(F90) -c $(FOPTS) $<

solve_ellip.o:	solve_ellip.f90
		$(F90) -c $(FOPTS) $<

clean:
		rm -rf *.o *.mod *.exe

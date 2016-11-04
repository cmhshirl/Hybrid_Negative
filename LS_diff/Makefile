############################################################################
#
##  Module:          Makefile
#
##  Purpose:         Top-level Makefile
#
##  Creation date:   Nov 27, 2015
#
##  Modified:
#
##  Send bug reports, comments or suggestions to 
#
#############################################################################
#~~~> gfortran (GNU FORTRAN Compiler)
export FC = gfortran
export FFLAGS = -cpp -g -ffree-line-length-none -O3 -fcheck=bounds -Wall -fbacktrace

vpath %.f90 ./
OUTNAME = cell_lsodar

default: driver clean

modRand.o: modRand.f90
	$(FC) $(FFLAGS) -c $<

modPar.o: modPar.f90
	$(FC) $(FFLAGS) -c $<

modHelper.o: modHelper.f90 modPar.o
	$(FC) $(FFLAGS) -c $<

modODE.o: modODE.f90 modPar.o modRand.o
	$(FC) $(FFLAGS) -c $<

modSSA.o: modSSA.f90 modODE.o 
	$(FC) $(FFLAGS) -c $<

modSolver.o: modSolver.f90 modODE.o
	$(FC) $(FFLAGS) -c $< 

modSimulate.o: modSimulate.f90 modHelper.o modSSA.o modSolver.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 modSimulate.o
	$(FC) $(FFLAGS) -c $< 


default: driver

driver: main.o modRand.o modPar.o modSSA.o modODE.o modSolver.o dlsodar.a modHelper.o modSimulate.o
	$(FC) $(FFLAGS) -o $(OUTNAME) $< modRand.o modPar.o modSSA.o modODE.o modSolver.o dlsodar.a modHelper.o modSimulate.o

purge:  clean

clean:
	rm -f *~ *.o *.mod
	rm -f driver
help:
	@$(ECHO) "usage: make ?"


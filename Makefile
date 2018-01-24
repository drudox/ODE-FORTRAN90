## =====================================================================
## Makefile for the eqn 001 :
## f = 10*(t-1)*y with y0=exp(-5), t [0,2] 

#-----------------------------------------------------------------------
# declaration
#-----------------------------------------------------------------------

# compiler 
#
FC = gfortran

# flags for debugging or maximal performance

FCFLAGS = -g -fbounds-check
FCFLAGS = -O2 

# flags forall (e.g. look for system .mod files required in gfortran)

FCFLAGS += -I/usr/include

# executable to be built 
#PROGRAMS = eq1 eq Euler/Euler MultiStep/MultiStep RungeKutta/RungeKutta jacobian main2

PROGRAMS = main4

#eq1 eq Euler/Euler MultiStep/MultiStep RungeKutta/RungeKutta jacobian main2
# "make" builds all 

all: $(PROGRAMS)

#-----------------------------------------------------------------------
### GENERAL RULES 

main4.o : eq4.o 
main4.o : Euler/Euler.o 
main4.o : MultiStep/MultiStep.o 
main4.o : MultiStep/AdamsMethods/AdamsBashforth.o
main4.o : MultiStep/AdamsMethods/AdamsMoulton.o
main4.o : RungeKutta/RungeKutta.o
main4.o : Jacobian/jacobian.o 


main4: eq4.o
main4: Euler.o
main4: MultiStep.o
main4: RungeKutta.o
main4: jacobian.o 
main4: AdamsBashforth.o
main4: AdamsMoulton.o


# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD





FC = gfortran
FFLAGS = -g
SRC = mc_subroutines.f90 potential_LJ.f90 potential_LJ_NL.f90 potential_Sti.f90 potential_Sti_NL.f90 main.f90

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

MCS_AL.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o MCS_AL.exe


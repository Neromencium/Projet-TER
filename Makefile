#
#
F90 = gfortran -ffree-line-length-none
DEBUG_FLAG = -g -fcheck=all -fbounds-check -fbacktrace -C
OPTIM_FLAG = -O2
#USUAL_FLAG = -O3

PROG = run
SRC = precision.f90 type_def.f90 public.f90 \
condi_ini.f90 eos_dt.f90 ecriture.f90 \
condi_lim.f90 temps.f90 main.f90

usual :
	$(F90) $(USUAL_FLAG) $(SRC) -o $(PROG)

debug :
	$(F90) $(DEBUG_FLAG) $(SRC) -o $(PROG)

optim :
	$(F90) $(OPTIM_FLAG) $(SRC) -o $(PROG)

clean :
	@rm *.o *.mod *~ core a.out $(PROG)
	@echo "On a fait du nettoyage"

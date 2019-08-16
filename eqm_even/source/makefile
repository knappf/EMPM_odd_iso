# compiler 
XLF= ifort -O3 -parallel #-O3 -parallel
#XLF= ifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= eqm_even
#
all: input_sp.o anglib.o base.o choles.o ham.o admat_cp_par.o dmat_test.o eqm_even.o 
	$(XLF) input_sp.o anglib.o base.o choles.o ham.o admat_cp_par.o dmat_test.o eqm_even.o ${OPT} -o ../${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

base.o: base.f90
	${XLF} -c base.f90

choles.o: choles.f
	${XLF} -c choles.f

ham.o: ham.f
	${XLF} -c ham.f

admat_cp_par.o: admat_cp_par.f90
	${XLF} -c admat_cp_par.f90

dmat_test.o: dmat_test.f
	${XLF} -c dmat_test.f

eqm_even.o: eqm_even.f90
	${XLF} -c eqm_even.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}




	

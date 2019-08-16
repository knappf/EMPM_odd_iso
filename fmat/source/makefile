# compiler 
XLF= ifort -O3 -parallel 
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= fmat
#
all: anglib.o fmat_unfi.o 
	$(XLF) anglib.o fmat_unfi.o ${OPT} -o ../${EXF}

anglib.o: anglib.f90
	${XLF} -c anglib.f90

fmat_unfi.o: fmat_unfi.f90
	${XLF} -c fmat_unfi.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     




	

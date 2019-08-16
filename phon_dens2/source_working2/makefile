# compiler 
#XLF= ifort -g #-O3 -parallel #-O3 -parallel
XLF= mpiifort -O3 -mcmodel=large #-g -debug full #-O3 -parallel 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= phon_dens2_MPI
#
all: anglib.o ro4_MPI.o input_sp.o phon_dens_MPI.o
	$(XLF) anglib.o ro4_MPI.o input_sp.o phon_dens_MPI.o ${OPT} -o ../${EXF}

anglib.o: anglib.f90
	${XLF} -c anglib.f90

ro4_MPI.o: ro4_MPI.f90
	${XLF} -c ro4_MPI.f90

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

phon_dens_MPI.o: phon_dens_MPI.f90
	${XLF} -c phon_dens_MPI.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     


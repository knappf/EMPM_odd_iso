# compiler 
XLF= ifort -g #-O3 -parallel 
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= tda
#
all: input_sp.o cm.o tda_pn_coupl_ortog.o
	$(XLF) input_sp.o cm.o tda_pn_coupl_ortog.o ${OPT} -o ../${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

cm.o: cm.f90
	${XLF} -c cm.f90

tda_pn_coupl_ortog.o: tda_pn_coupl_ortog.f90
	${XLF} -c tda_pn_coupl_ortog.f90


clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     




	

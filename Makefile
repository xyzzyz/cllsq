F90 = ifort -free -u -debug -diag-enable warn -warn all
INCLUDE = -I/opt/intel/composer_xe_2011_sp1.8.273/mkl/include/intel64/lp64
LDFLAGS = -L/opt/intel/composer_xe_2011_sp1.8.273/mkl/lib/intel64
LINK = -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

llsq: llsq.o
	${F90} llsq.o -o llsq ${LDFLAGS} ${LINK}

llsq.o: llsq.f90 
	${F90} -c llsq.f90 -o llsq.o ${INCLUDE} 

clean:
	rm -f *.o *.mod llsq

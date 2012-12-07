# Makefile

FC = mpif90     # Fortran 90 compiler
F90 = ifort     # Fortran 90 compiler (for blas)

#OPT = -O0 -r8
# OPT = -O3 -CB -r8
OPT = -O3 -r8

MMOD = -module

MODULEDIR = ./Modules/

F90FLAGS = $(MMOD) $(MODULEDIR)


OBJ = $(patsubst %.f90,%.o, $(wildcard *.f90) ) 
OBJblas = ./BLAS/dgemm.o ./BLAS/lsame.o ./BLAS/xerbla.o
OBJ_MOD = $(patsubst %.f90,%.o,$(wildcard $(MODULEDIR)*.f90) )

all : spec.x

spec.x : make_module $(OBJ) $(OBJblas)
	$(FC) $(F90FLAGS) $(OPT) -o $@ $(OBJ) $(OBJ_MOD) $(OBJblas)

./BLAS/dgemm.o:	./BLAS/dgemm.f
			$(F90) -c ./BLAS/dgemm.f -o ./BLAS/dgemm.o

./BLAS/lsame.o:	./BLAS/lsame.f
			$(F90) -c ./BLAS/lsame.f -o ./BLAS/lsame.o

./BLAS/xerbla.o:	./BLAS/xerbla.f
			$(F90) -c ./BLAS/xerbla.f -o ./BLAS/xerbla.o

make_module : 
	$(MAKE) all -C $(MODULEDIR) "OPT=$(OPT)" "F90MOD=$(MMOD)"

%.o : %.f90 
	$(FC) $(F90FLAGS) $(OPT) -c $< 
clean: clean_modules
	rm -f *.o ; rm -f ./BLAS/*.o
cleanx:
	rm -f *.x
cleanall : clean cleanx

clean_modules:
	$(MAKE) clean -C $(MODULEDIR)


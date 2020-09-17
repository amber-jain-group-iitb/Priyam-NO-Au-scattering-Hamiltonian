all:
	gfortran -c project_new.f90
	gfortran -c mod_ham_modular.f90
	gfortran -g -ffpe-trap=zero,invalid,overflow,underflow project_new.o mod_ham_modular.o
	

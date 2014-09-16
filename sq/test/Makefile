LIB= ../src/libsq.a   
all:
	gcc -O3  -c fgamma_test.c -I../src
	gfortran -O3 -o fgamma_test fgamma_test.o  -lm ${LIB}  -llapack -lblas
	gcc -O3 -c int_test.c   -I../src
	gfortran -O3 -o int_test  int_test.o       -lm ${LIB}
	gcc -O3 -c basis_test.c -I../src
	gfortran -O3 -o basis_test  basis_test.o   -lm ${LIB}
	gcc -O3 -c lin_test.c -I../src
	gfortran -O3 -o lin_test    lin_test.o     -lm ${LIB}
	gcc -O3 -c matrix_test.c -I../src
	gfortran -O3 -o matrix_test matrix_test.o  -lm ${LIB}
	rm -f *.o

clean:
	rm -f fgamma_test int_test basis_test lin_test matrix_test


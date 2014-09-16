SRC = sq.c fgamma.c int.c util.c basis.c mol.c lin.c matrix.c rhf.c uhf.c option.c xsf.c conv.c grad.c optimize.c pop.c check.c mp2.c multipole.c quartet.c rpc.c mecp.c
OBJ = sq.o fgamma.o int.o util.o basis.o mol.o lin.o matrix.o rhf.o uhf.o option.o xsf.o conv.o grad.o optimize.o pop.o check.o mp2.o multipole.o quartet.o rpc.o mecp.o
LIB = ../../lapack_LINUX.a ../../blas_LINUX.a

#CC  = icc
#CCOPT = -O3 -funroll-loops -fno-alias -m64 -Wall
#LDOPT = -O3 -m64 -Wall

CC  = gcc
CCOPT = -O3 -m64 -Wall -msse2 -funroll-loops
LDOPT = -O3 -m64 -Wall

all: sq

%.o: %.c
	${CC} ${CCOPT} -c $<

sq: ${OBJ}
	${CC} ${LDOPT} -o sq $^ -lm -lgfortran ${LIB}
#	gfortran -o sq $^ -lm -static ${LIB}

clean:
	rm -f sq
	rm -f *.o

#libsq.so: ${SRC}
#	${CC} -Wall -shared -O3 -o libsq.so ${SRC}

#all: libsq.a sq

#libsq.a:${SRC}
#	${CC} ${CCOPT} -c *.c
#	ar cr libsq.a ${OBJ}
#	rm -f *.o	

#sq: sq.c libsq.a
#	${CC} ${CCOPT} -c sq.c
#	${CC} ${LDOPT} -o sq sq.o libsq.a -lm -lgfortran ${LIB} 
#	rm -f sq.o

#clean:
#	rm -f sq
#	rm -f libsq.a

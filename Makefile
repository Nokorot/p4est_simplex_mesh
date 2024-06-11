
CC=mpicc

CFLAGS=-lp4est -lsc -L../local/lib \
			 -I../local/include -I${MPI_INCLUDE} -I${MPI_LIB} \
			 -Wl,-rpath="${HOME}/Thisis/local/lib"

.PHONY: alla

all: simplex2 simplex3

simplex2: simplex2.c p4est_simplex_mesh.o
	${CC} ${CFLAGS} -o $@ $^

simplex3: simplex3.c p8est_simplex_mesh.o
	${CC} ${CFLAGS} -o $@ $^

p4est_simplex_mesh.o: p4est_simplex_mesh.c p4est_simplex_mesh.h
	${CC} ${CFLAGS} -c -o $@ $<

p8est_simplex_mesh.o: p4est_simplex_mesh.c p4est_simplex_mesh.h
	${CC} ${CFLAGS} -c -o $@ $<

	


CC=mpicc

CFLAGS=-ggdb -lp4est -lsc -L../local/lib \
			 -I../local/include -I${MPI_INCLUDE} -I${MPI_LIB} \
			 -Wl,-rpath="${HOME}/Thisis/local/lib" \
			 -DP4EST_SIMPLEX_DEBUG \
			 -DP4EST_ENABLE_DEBUG

.PHONY: alla

all: simplex2 simplex3

simplex2: simplex2.o p4est_simplex_mesh.o
	${CC} ${CFLAGS} -o $@ $^

simplex3: simplex3.o p8est_simplex_mesh.o
	${CC} ${CFLAGS} -o $@ $^



simplex2.o: simplex2.c
	${CC} ${CFLAGS} -c -o $@ $<

simplex3.o: simplex3.c simplex2.c
	${CC} ${CFLAGS} -c -o $@ $<

p4est_simplex_mesh.o: p4est_simplex_mesh.c p4est_simplex_mesh.h
	${CC} ${CFLAGS} -c -o $@ $<

p8est_simplex_mesh.o: p8est_simplex_mesh.c p8est_simplex_mesh.h p4est_simplex_mesh.c p4est_simplex_mesh.h
	${CC} ${CFLAGS} -c -o $@ $<



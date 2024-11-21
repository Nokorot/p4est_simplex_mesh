
PREFIX=${PWD}/../local
CC=mpicc

CFLAGS=-ggdb -lp4est -lsc -L"${PREFIX}/lib" \
			 -I"${PREFIX}/include" -I"${MPI_INCLUDE}" -I${MPI_LIB} \
			 -Wl,-rpath="${PREFIX}/lib" \
			 -DP4EST_SIMPLEX_DEBUG \
			 -DP4EST_ENABLE_DEBUG \
			 -DTIMINGS

# CFLAGS=-ggdb -lp4est -lsc -L../local/lib \
# 			 -I../local/include -I${MPI_INCLUDE} -I${MPI_LIB} \
# 			 -Wl,-rpath="${HOME}/Thisis/local/lib" \
# 			 -DP4EST_SIMPLEX_DEBUG \
# 			 -DP4EST_ENABLE_DEBUG

.PHONY: all, clean

all: simplex2 simplex3

clean:
	rm -f *.o simplex2 simplex3

simplex2: simplex2.o p4est_simplex_mesh.o statistics.o
	mkdir -p out
	${CC} ${CFLAGS} -o $@ $^

simplex3: simplex3.o p8est_simplex_mesh.o statistics.o
	mkdir -p out
	${CC} ${CFLAGS} -o $@ $^


simplex2.o: simplex2.c statistics.h
	${CC} ${CFLAGS} -c -o $@ $<

simplex3.o: simplex3.c simplex2.c statistics.h
	${CC} ${CFLAGS} -c -o $@ $<

p4est_simplex_mesh.o: p4est_simplex_mesh.c p4est_simplex_mesh.h statistics.h
	${CC} ${CFLAGS} -c -o $@ $<

p8est_simplex_mesh.o: p8est_simplex_mesh.c p8est_simplex_mesh.h p4est_simplex_mesh.c p4est_simplex_mesh.h statistics.h
	${CC} ${CFLAGS} -c -o $@ $<

statistics.o: statistics.c statistics.h
	${CC} ${CFLAGS} -c -o $@ $<

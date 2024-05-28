
PREFIX=../local


.PHONY: all
all: simplex2 simplex3

simplex2: simplex2.c p4est_simplex_to_gmsh.c
	gcc -o $@ $^ \
		-ggdb \
		-I$(PREFIX)/include \
		-I$(MPI_INCLUDE) \
		-L$(PREFIX)/lib \
		-L$(MPI_LIB) \
		-lp4est -lsc -lmpi -Wl,-rpath,$(PREFIX)/lib

simplex3: simplex2.c p4est_simplex_to_gmsh.c
	gcc -o $@ simplex3.c p8est_simplex_to_gmsh.c \
		-ggdb \
		-I$(PREFIX)/include \
		-I$(MPI_INCLUDE) \
		-L$(PREFIX)/lib \
		-L$(MPI_LIB) \
		-lp4est -lsc -lmpi -Wl,-rpath,$(PREFIX)/lib

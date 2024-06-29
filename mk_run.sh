#!/bin/sh

N=5

make -j && mpirun -n $N ./simplex2 -l1 -L2 \
  | tee out.log \
  && for i in {0..$N}; do \
        < out.log grep "\[[$i]\]"; echo ""; \
     done | tee out2.log | less

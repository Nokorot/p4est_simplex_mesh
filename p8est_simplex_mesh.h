#ifndef P8EST_SIMPLEX_MESH_H_
#define P8EST_SIMPLEX_MESH_H_

// #include <p4est_to_p8est.h>

#define p4est_simplex_mesh_t            p8est_simplex_mesh_t

#define p4est_simplex_mesh_destroy      p8est_simplex_mesh_destroy
#define p4est_new_simplex_mesh          p8est_new_simplex_mesh

#define p4est_simplex_mesh_write_gmsh_file \
        p8est_simplex_mesh_write_gmsh_file
#define p4est_simplex_mesh_write_vtk_file \
        p8est_simplex_mesh_write_vtk_file

#include "p4est_simplex_mesh.h"

#endif

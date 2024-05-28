#ifndef SIMPLEX_MESH_TO_GMSH_H
#define SIMPLEX_MESH_TO_GMSH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sc_containers.h"

// Define the p4est_simplex_mesh_t structure
typedef struct {
//     sc_array_t *element_nodes;  // Array of vertices, each vertex is 3 doubles (x, y, z)
//     sc_array_t *corner_nodes;   // Array of vertices, each vertex is 3 doubles (x, y, z)
//     sc_array_t *face_nodes;     // Array of vertices, each vertex is 3 doubles (x, y, z)
// #ifdef P4_TO_P8
//     sc_array_t *edge_nodes;     // Array of vertices, each vertex is 3 doubles (x, y, z)
// #endif

    sc_array_t *vertices;   // Array of vertices, each vertex is 3 doubles (x, y, z)
    sc_array_t *simplicies; // Array of simplices, each simplex is 4 p4est_locidx_t
} p4est_simplex_mesh_t;

// Function to write the mesh to a Gmsh file
void write_gmsh_file(const char *filename, p4est_simplex_mesh_t *mesh);

// Function to write the mesh to a VTK file
void write_vtk_file(const char *filename, p4est_simplex_mesh_t *mesh);

#endif // SIMPLEX_MESH_TO_GMSH_H

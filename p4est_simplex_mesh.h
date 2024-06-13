#ifndef P4EST_SIMPLEX_MESH_H_
#define P4EST_SIMPLEX_MESH_H_

#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_lnodes.h>
#include <p4est_ghost.h>
#include <p4est_geometry.h>
#else
#include <p8est_bits.h>
#include <p8est_lnodes.h>
#include <p8est_ghost.h>
#include <p8est_geometry.h>
#endif

// Define the p4est_simplex_mesh_t structure
typedef struct {
  
  union {
    struct {
      // lnode
      sc_MPI_Comm         mpicomm;
      p4est_locidx_t      num_local_nodes;
      p4est_locidx_t      owned_count;
      p4est_gloidx_t      global_offset;
      p4est_gloidx_t     *nonlocal_nodes;
      sc_array_t         *sharers;
      p4est_locidx_t     *global_owned_count;

      int                 degree, vnodes; // Not defined
      p4est_locidx_t      num_local_elements;
      p4est_lnodes_code_t *face_code; //  Not defined
      p4est_locidx_t     *element_nodes; // Not YET defined
      // end lnode
    };
    p4est_lnodes_t lnodes;
  };

    sc_array_t *vertices;   // Array of vertices, each vertex is 3 doubles (x, y, z)
    sc_array_t *simplicies; // Array of simplices, each simplex is 4 p4est_locidx_t
} 
p4est_simplex_mesh_t;

p4est_simplex_mesh_t *
p4est_new_simplex_mesh(
    p4est_t *p4est,
    p4est_geometry_t *geometry,
    p4est_ghost_t *ghost);

void
p4est_simplex_mesh_destroy(
    p4est_simplex_mesh_t *smesh);

void 
p4est_quad_corner_coords(
    p4est_quadrant_t *quad,
    int face,
    double abc[3]);

void 
p4est_quad_center_coords(
    p4est_quadrant_t *quad,
    double abc[3]);

void 
p4est_quad_face_center_coords(
    p4est_quadrant_t *quad,
    int face,
    double abc[3]);

#ifdef P4_TO_P8
void
p8est_quad_edge_center_coords(
    p4est_quadrant_t *quad,
    int edge,
    double abc[3]);
#endif

// Function to write the mesh to a Gmsh file
void p4est_simplex_mesh_write_gmsh_file(
    const char *filename, p4est_simplex_mesh_t *mesh);

// Function to write the mesh to a VTK file
void p4est_simplex_mesh_write_vtk_file(
    const char *filename, p4est_simplex_mesh_t *mesh);

#endif // P4EST_SIMPLEX_MESH_H_

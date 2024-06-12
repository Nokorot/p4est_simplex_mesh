
#ifdef P4_TO_P8
#define QUAD_COORDS(quad) \
            (quad).x / (double) P4EST_ROOT_LEN, \
            (quad).y / (double) P4EST_ROOT_LEN, \
            (quad).z / (double) P4EST_ROOT_LEN

#define PRT_QUAD(quad) printf("(%.4f %.4f %.4f)", QUAD_COORDS(quad))

#else
#define QUAD_COORDS(quad) \
            (quad).x / (double) P4EST_ROOT_LEN, \
            (quad).y / (double) P4EST_ROOT_LEN

#define PRT_QUAD(quad) printf("(%.4f %.4f)", QUAD_COORDS(quad))
#endif

#define PRTN_QUAD(quad) PRT_QUAD(quad); printf("\n");

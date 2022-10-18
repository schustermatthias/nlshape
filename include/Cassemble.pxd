# In this file all declarations of the C functions happen, in order to make them visible to cimport.
# If this .pxd is compiled with setup_linalg.py it results in a standalone package which can be cimported.
# In order import (to python) see linalg.pyx

# distutils: language = c++
# distutils: sources = RoadDensityPyMod.cpp
from libcpp cimport bool
from libcpp.string cimport string
# Assembly algorithm with BFS -----------------------------------------------------------------------------------------

cdef extern from "Cassemble.h":
    void par_assemble(const string compute, const string path_spAd, const string path_fd, const int K_Omega, const int K,
                  const long *ptrTriangles, const long *ptrLabelTriangles, const double *ptrVerts, const long * ptrLabelVerts,
                  const int nE,
                  const int nE_Omega, const int L, const int L_Omega, const double *Px, const int nPx, const double *dx,
                  const double *Py, const int nPy, const double *dy, const double sqdelta,
                  const int is_DiscontinuousGalerkin, const int is_NeumannBoundary, const string str_model_kernel,
                  const string str_model_f, const string str_integration_method_remote,
                  const string str_integration_method_close, const int is_PlacePointOnCap,
                  const int dim, const int outdim,
                  const long * ptrZetaIndicator_indices,
                  const long * ptrZetaIndicator_indptr, long nZeta,
                  const double * Pg, const int nPg, const double * dg, double maxDiameter,
                  double fractional_s,
                  int is_fullConnectedComponentSearch, int verbose, int is_ShapeDerivative, const double *state,
                  const double *adjoint)
    # Mass matrix evaluation ----------------------------------------------------------------------------------------------
    void par_evaluateMass(double *vd, const double *ud, long *Elements,
                          const long *ElementLabels, const double *Verts, const long * VertexLabels, int K_Omega, int J, int nP,
                          double *P, const double *dx, int dim, int outdim, int is_DiscontinuousGalerkin)

    # DEBUG Helpers and test functions
    int method_retriangulate(const double * x_center, const double * TE,
                           const double sqdelta, double * re_Triangle_list,
                           int is_placePointOnCap)
    int method_retriangulateInfty(const double * xCenter, const double * TE,
                                double sqdelta, double * reTriangleList,
                                int isPlacePointOnCap)
    void toRef(const double * E, const double * phys_x, double * ref_p)    # Pull point to Reference Element (performs 2x2 Solve)
    void toPhys(const double * E, const double * p, int dim, double * out_x)
    void solve2x2(const double * A, const double * b, double * x)
    void get_div(const double * E, double * div)
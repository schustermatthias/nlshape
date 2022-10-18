/*! \mainpage
 *
 * \section intro_sec Introduction
 *
 * This library provides a parallel assembly routine for a specific class of integral operators.
 * The focus lies on operators of the form
 * \f[
 *  \mathcal{L}(\mathbf{u})(\mathbf{x}) = p.v. \int_{B_{\delta}(\mathbf{x}) \cap \widetilde{\Omega}}(\mathbf{C}_\delta(\mathbf{x}, \mathbf{y})  \mathbf{u}(\mathbf{x}) - \mathbf{C}_\delta(\mathbf{y}, \mathbf{x})\mathbf{u}(\mathbf{y}))  d\mathbf{y}.
 *  \f]
 * where \f$ \mathbf{C}(\mathbf{x},\mathbf{y}) \f$ is a, possibly strongly singular, matrix-valued integral kernel with bounded interaction
 * radius. The domain \f$\Omega\f$ and the interaction set \f$\Omega_I\f$
 * can be subsets of \f$ R^d \f$ for \f$ d = 2,3\f$. The solution can be scalar or
 * vector valued, i.e. \f$ u(x) \in R^c\f$ for \f$c \geq 1\f$.
 *
 * This implementation addresses researchers who want to verify findings numerically. The computational
 * burden of the assembly of integral operators is so large that even computations in an
 * experimental setting require efficient and parallel implementations. However, we do not aim
 * to provide the fastest possible but rather a flexible implementation which runs at
 * convenient speed.
 *
 * \section Assembly
 *
 * The Python interface
 * calls the function par_assemble() which is a wrapper function. This function again calls the function par_system()
 * where the assembly of the nonlocal stiffness matrix happens, or the function par_forcing() for the assembly
 * of the right hand side.
 * The wrapper function par_assemble() translates the input data into
 * three objects which are defined in MeshStruct, QuadratureStruct, and ConfigurationStruct. It also performs
 * some basic input checks. The output of the integration is saved as
 * [Armadillo sparse matrix](http://arma.sourceforge.net/docs.html#SpMat).
 * The file is then copied again into the Python interface to make it available.
 *
 *\subsection Integration
 *
 * There are different integration routines available (e.g. integrate_retriangulate(), integrate_baryCenter(),
 * ...). The file integration.cpp also contains the implementations for the underlying truncation methods
 * (e.g. method_retriangulate(), method_baryCenter(), ...).
 *
 * \subsection Model
 *
 * The kernel and forcing functions, as well as the basis functions are defined in model.cpp.
 * Different options for kernels (e.g. kernel_constant(),
 * kernel_linearPrototypeMicroelastic, ...), and forcing terms (e.g. f_constant(), f_linear(),...)
 * are implemented. If you want to add a kernel or forcing function just add the corresponding
 * function to model.cpp obeying the common function signature. In order to make it accessible as
 * option you have change the function lookup_configuration() accordingly.
 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#ifndef CASSEMBLE_H
#define CASSEMBLE_H
#include <armadillo>
#include <MeshTypes.h>
#include "cstring"
#include "model.h"

//TODO stop it!!
using namespace std;
//Put your stuff into a namespace

////// Retriangulation Routine ---------------------------------------------------------------------------------------------

int method_retriangulate(const double * xCenter, const double * TE,
                         double sqdelta, double * reTriangleList,
                         int isPlacePointOnCap);
int method_retriangulateInfty(const double * xCenter, const double * TE,
                              double sqdelta, double * reTriangleList,
                              int isPlacePointOnCap);
void toRef(const double * E, const double * phys_x, double * ref_p);     // Pull point to Reference Element (performs 2x2 Solve)
void toPhys(const double * E, const double * p, double * out_x);         // Push point to Physical Element
void toPhys(const double * E, const double * p, int dim, double * out_x);
void solve2x2(const double * A, const double * b, double * x);
void get_div(const double * E, double * div);
double test_integrate_shape(const string compute, const string path_spAd, const string path_fd, const int K_Omega, const int K,
                            const long *ptrTriangles, const long *ptrLabelTriangles, const double *ptrVerts, const long * ptrLabelVerts, const int nE,
                            const int nE_Omega, const int nV, const int nV_Omega, const double *Px, const int nPx, const double *dx,
                            const double *Py, const int nPy, const double *dy, const double sqdelta,
                            const int is_DiscontinuousGalerkin, const int is_NeumannBoundary, const string str_model_kernel,
                            const string str_model_f, const string str_integration_method_remote,
                            const string str_integration_method_close,
                            const int is_PlacePointOnCap,
                            const int dim, const int outdim,
                            const long * ptrZetaIndicator_indices,
                            const long * ptrZetaIndicator_indptr,
                            const long nZeta,
                            const double * Pg, const int degree, const double * dg, double maxDiameter, double fractional_s,
                            int is_fullConnectedComponentSearch, int verbose,
                            int is_ShapeDerivative, const double *state, const double *adjoint, int triangle_i, int triangle_j);

////// Assembly algorithm with BFS -----------------------------------------------------------------------------------------
/*!
 * @brief Parallel assembly of nonlocal operator using a finite element approach. This function is a wrapper
 * for the functions par_system() and par_forcing(). It allows to call the C++ functions from Cython without
 * writing wrapper classes or structs.
 *
 * @details Any 2-dimensional array is handed over by a pointer only. The expected shape of the input arrays is given
 * in the parameter list where the underlying data is expected to be in C-contiguous
 * (i.e. row-major) order. We denote the dimension of the domain by d. See par_system() and par_forcing() for
 * mor details on the assembly.
 *
 * @param compute   string "system", "forcing", "systemforcing". This string determines whether the stiffnes matrix,
 * load or both should be computed.
 * @param path_spAd Save path for sparse matrix (armadillo binary)
 * @param path_fd   Save path for forcing (armadillo binary)
 * @param K_Omega   Number of rows of the stiffness matrix A. Example: If you want to solve a scalar equation using
 * continuous Galerkin basis functions then K_Omega is equal to the number of basisfunctions (nV_Omega) which are not part of the
 * Dirichlet boundary. If your problem is not scalar as for example in peridynamics then K_Omega = nV_Omega*outdim.
 * @param K Number of Columns of the stiffness matrix A
 * @param ptrTriangles <B>(nE, d+1)</B> Pointer to elements. Label 1 for elements in Omega, Label 2 for elements in OmegaI.
 * @param ptrLabelTriangles <B>(nE,)</B> Pointer to element labels
 * @param ptrVerts <B>(L, d)</B> Pointer to vertices
 * @param ptrLabelVerts Pointer to vertex labels
 * @param nE         Number of elements
 * @param nE_Omega   Number of elements in Omega
 * @param nV         Number of vertices
 * @param nV_Omega   Number of vertices in Omega
 * @param Px        <B>(nPx, d)</B> Pointer to quadrature points for the outer integral
 * @param nPx       Number of quadrature points in the outer integral
 * @param dx        <B>(nPx,)</B> Pointer to quadrature weights of the outer integral
 * @param Py        Number of quadrature points in the inner integral
 * @param nPy       <B>(nPy, d)</B> Pointer to quadrature points for the inner integral
 * @param dy        <B>(nPx,)</B> Pointer to quadrature weights of the inner integral
 * @param sqdelta   Squared delta
 * @param is_DiscontinuousGalerkin Switch for discontinuous Galerkin
 * @param is_NeumannBoundary Switch of Neumann Boundary Conditions
 * @param str_model_kernel  Name of kernel
 * @param str_model_f   Name of right hand side
 * @param str_integration_method_remote Name of integration method
 * @param is_PlacePointOnCap Switch for withcaps parameter in retriangulation
 * @param dim       Dimension of the domain
 * @param outdim    Dimension of the solution space (e.g. 1 for scalar problems, dim for linear elasticity)
 * @param ptrZeta   Pointer to overlap counter of decomposed Mesh (optional)
 * @param nZeta     Number of rows of Zeta
 * @param Pg        <B>(nPg, dim^2)</B> Quadrature points for tensor Gauss quadrature (optional, needed for singular kernels).
 * @param nPg       <B>(tensorGaussDegree^dim,)</B> Number of quadrature points.
 * @param dg        <B>(nPg,)</B> Weights for tensor Gauss quadrature.
 * @param maxDiameter Maximal diameter of finite elements (optional). Might increase speed of retriangulation if provided.
 * @param fractional_s  Degree of fractional kernel (d=2 only). Required for the correct choice of the integration routine.
 *                      Default is s=-1.0, which corresponds to a kernel with no singularity.
 * @param verbose    switch for verbose mode.
 * */
void par_assemble(string compute, string path_spAd, string path_fd, int K_Omega, int K,
                  const long *ptrTriangles, const long *ptrLabelTriangles, const double *ptrVerts, const long * ptrLabelVerts,
                  int nE,
                  int nE_Omega, int nV, int nV_Omega, const double *Px, int nPx, const double *dx,
                  const double *Py, int nPy, const double *dy, double sqdelta,
                  //const long *ptrNeighbours,
                  //int nNeighbours,
                  int is_DiscontinuousGalerkin, int is_NeumannBoundary, string str_model_kernel,
                  string str_model_f, string str_integration_method_remote,
                  string str_integration_method_close,
                  int is_PlacePointOnCap,
                  int dim, int outdim,
                  const long * ptrZetaIndicator_indices,
                  const long * ptrZetaIndicator_indptr, long nZeta = 0,
                  const double * Pg = nullptr, int degree = 0, const double * dg = nullptr, double maxDiameter = 0.0,
                  double fractional_s=-1.0, int is_fullConnectedComponentSearch=0, int verbose=0,
                  int is_ShapeDerivative=0, const double *state = nullptr, const double *adjoint = nullptr);

/**
 * @brief Parallel assembly of nonlocal operator using a finite element approach.
 * Kernel functions can be defined in *model* see model_kernel() for more information.
 *
 * This function assembles the stiffness
 * matrix corresponding to the operator
 *
 *  * \f[
 *  \mathcal{L}(\mathbf{u})(\mathbf{x}) = p.v. \int_{B_{\delta}(\mathbf{x}) \cap \widetilde{\Omega}}(\mathbf{C}_\delta(\mathbf{x}, \mathbf{y})  \mathbf{u}(\mathbf{x}) - \mathbf{C}_\delta(\mathbf{y}, \mathbf{x})\mathbf{u}(\mathbf{y}))  d\mathbf{y}.
 *  \f]
 *
 *  It traverses all elements aT with element Label != 0 and adds up their contribution to the global stiffness matrix. The integration starts with the domain aT x aT and then proceeds with the neighbouring
 *  elements of aT in the mesh until the interaction domain is exceeded. For each pair aT, bT an integration routine
 *  (options are e.g. integrate_retriangulate(), integrate_baryCenter(), integrate_baryCenterRT(), ...)
 *  is called. If it all computed integrals are 0 the elements are considered as non-interacting. The integrals,
 *  and the interaction sets
 *  again depend on the truncation routines, which are called inside integrate(). This approach allows to
 *  define interaction sets directly in the truncation routines. The code then automatically
 *  finds the interacting elements bT without traversing all of them. This approach requires setting up the dual
 *  graph of the mesh which is to be found in mesh.neighbours.
 *
 * @param Ad, Pointer to map<unsigned long, double> which is used to store matrix values. This map is shared between all threads and the value updates
 * happen in critical sections. This makes sense because in case of nonlocal operators Ad can turn out to be very large.
 * Storing Ad separately in separate threads and finally merging the different maps results in a significant overhead.
 * @param mesh The mesh is of type MeshType and contains all information about the finite element descretization.
 * See MeshStruct for more information.
 * @param quadRule Quadrature rules for inner, and outer elements as well as for the singular kernels.
 * @param conf General configuration, namely kernel, and forcing functions, as well as integration method.
 */
void par_system(MeshType &mesh, QuadratureType &quadRule, ConfigurationType &conf,
                const double *state = nullptr, const double *adjoint = nullptr);

/**
 * @brief Parallel assembly of forcing term. Forcing functions can be defined in *model* see model_f() for
 * more information.
 *
 * @param mesh The mesh is of type MeshType and contains all information about the finite element descretization.
 * See MeshStruct for more information.
 * @param quadRule Quadrature rules for inner, and outer elements as well as for the singular kernels.
 * @param conf General configuration, namely kernel, and forcing functions, as well as integration method.
 */
void par_forcing(MeshType &mesh, QuadratureType &quadRule, ConfigurationType &conf);

// Mass matrix evaluation ----------------------------------------------------------------------------------------------
/*!
 * @brief Evaluate the mass matrix v = Mu.
 *
 * @param vd        Pointer to the first entry of the output vector.
 * @param ud        Pointer to the first entry of the input vector.
 * @param Elements  List of elements of a finite element triangulation (CSR-format, row major order).
 * @param ElementLabels List of element Labels.
 * @param Verts     List of vertices (row major order).
 * @param VertexLabels List of vertex Labels.
 * @param K_Omega   Number of rows and columns in M. Example: If you use continuous Galerkin basis functions and
 * want to solve a scalar problem K_Omega = J.
 * @param J         Number of elements in the triangulation.
 * @param nP        Number of quadrature points in the outer integral
 * @param P         <B>(nPx, d)</B> Pointer to quadrature points.
 * @param dx        <B>(nPx,)</B> Pointer to quadrature weights.
 * @param dim       Dimension of the domain Omega (2 or 3).
 */
void par_evaluateMass(double *vd, const double *ud, long *Elements,
                      const long *ElementLabels, const double *Verts, const long * VertexLabels, int K_Omega, int J, int nP,
                      double *P, const double *dx, int dim, int outdim, int is_DiscontinuousGalerkin);
//[DEBUG]
#endif /* Cassemble.h */
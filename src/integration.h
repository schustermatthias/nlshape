//
// Created by klar on 18.09.20.
//

#ifndef NONLOCAL_ASSEMBLY_INTEGRATION_H
#define NONLOCAL_ASSEMBLY_INTEGRATION_H
// ___ INTEGRATION DECLARATION _________________________________________________________________________________________
// This headers have been commented out due to compatibility issues with a system.
// Install the required C++ library, uncomment this and the function method_retriangulateInfty in
// integration.cpp to get access to the Linfinity-retriangulation

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

// Library for integration of L-infty Ball
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>         Triangulation;
typedef Triangulation::Point             Point;
typedef Triangulation::Finite_faces_iterator Finite_face_iterator;

// Integration Routine #################################################################################################

int integrate(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                  const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                  double *termLocalPrime, double *termNonlocPrime);
extern int (*integrate_remote)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                         const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                         double *termLocalPrime, double *termNonlocPrime);
extern int (*integrate_close)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                        const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                        double *termLocalPrime, double *termNonlocPrime);
extern int (*method)(const double * xCenter, const ElementType & T, const MeshType & mesh, double * reTriangleList,
              int isPlacePointOnCap);

int integrate_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                    const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                    double *termLocalPrime, double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values);
extern int (*integrate_remote_shape)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                               const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                               double *termLocalPrime, double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values);
extern int (*integrate_close_shape)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                              const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                              double *termLocalPrime, double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values);

int ERROR_wrongAccess(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
              const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
              double *termLocalPrime, double *termNonlocPrime);
int ERROR_wrongAccess_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                            const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                            double *termLocalPrime, double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values);
// Integration Methods #################################################################################################
// Methods -------------------------------------------------------------------------------------------------------------
/**
 * @brief This integration routines uses method_retriangulate() to truncate the *inner domain* bT. See integrate()
 * for general information about the integration routines.
 *
 *  termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 *  termNonloc = int_aT phiA(x) int_bT phiB(y) ker'(y,x) dy dx.
 *  termLocalPrime = int_aT  int_bT phiA(y) phiB(y) ker'(x,y) dy dx,\n
 *  termNonlocPrime = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int integrate_retriangulate(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                             const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                             double *termNonloc, double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief This integration routines uses method_retriangulate() to truncate the *inner domain* bT. See integrate()
 * for general information about the integration routines. The approx balls are not symmetrified.
 *
 *  termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 *  termNonloc = int_aT phiA(x) int_bT phiB(y) ker'(y,x) dy dx.
 *  termLocalPrime = int_aT  int_bT phiA(y) phiB(y) ker'(x,y) dy dx,\n
 *  termNonlocPrime = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */

int integrate_retriangulate_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                                  const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                                  double *termLocal, double *termNonloc, double *termLocalPrime,
                                  double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values);
/**
 * @brief This integration routines uses method_retriangulate() to truncate the *inner domain* bT. See integrate()
 * for general information about the integration routines. The approx balls are not symmetrified.
 *
 *  termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 *  termNonloc = int_aT phiA(x) int_bT phiB(y) ker'(y,x) dy dx.
 *  termLocalPrime = int_aT  int_bT phiA(y) phiB(y) ker'(x,y) dy dx,\n
 *  termNonlocPrime = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */

int integrate_retriangulate_unysmm(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                                   const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                                   double *termNonloc, double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief This integration routines uses method_exact() to truncate the *inner domain* bT. See integrate()
 * for general information about the integration routines.
 *
 *  termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 *  termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int integrate_exact(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                     const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                     double *termNonloc, double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief . This integration routines uses method_baryCenter() to truncate the *inner domain* bT. See integrate()
 * for general information about the integration routines.
 *
 * termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx\n,
 * termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int integrate_baryCenter(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                          const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                          double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief This integration routines uses method_retriangulate() to truncate the *outer domain* bT. See integrate()
 * for general information about the integration routines.
 *
 * termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 * termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int integrate_baryCenterRT(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                            const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                            double *termLocal, double *termNonloc, double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief This integration routines uses method_subSuperSetBalls() to truncate the *inner domain* bT. See integrate()
 * for general information about the integration routines.
 *
 * termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 * termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int
integrate_subSuperSetBalls(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                           const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                           double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief This integration routine performs no truncation of any domain. It can be applied to integrate in cases where
 * no truncation is needed. See integrate()
 * for general information about the integration routines.
 *
 * termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 * termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int integrate_fullyContained(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                              const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                              double *termLocal, double *termNonloc, double *termLocalPrime, double *termNonlocPrime);
/**
 * @brief This integration routine handles the weak singularity close to the origin.
 * Due to the transformation of the set aT x bT truncations are not possible. However,
 * the truncation is necessary only close to the boundary of the interaction set. There, this function
 * is used in combination with integration routines which truncate the domain.
 *
 * termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx,\n
 * termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.
 *
 * Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
 * matrix.
 * @param aT    Triangle of the outer integral.
 * @param bT    Triangle of the inner integral.
 * @param quadRule Quadrature rule.
 * @param mesh  Mesh.
 * @param conf  Confuration.
 * @param is_firstbfslayer (Unused)
 * @param termLocal This term contains the local part of the integral
 * @param termNonloc This term contains the nonlocal part of the integral
 * @param termLocalPrime This term contains the local part of the integral
 * @param termNonlocPrime This term contains the nonlocal part of the integral
 */
int integrate_weakSingular(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                            const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                            double * termLocal, double * termNonloc,
                            double *termLocalPrime, double *termNonlocPrime);

int integrate_fractional(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                           const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                           double * termLocal, double * termNonloc,
                           double *termLocalPrime, double *termNonlocPrime);
// Helpers -------------------------------------------------------------------------------------------------------------
/**
 * @brief This truncation method checks whether the distance of the point point x_center (of aT) to
 * the bary center of bT is smaller than delta. The function writes the list of triangles
 * which are obtained from the trunaction into reTriangle_list. In this case the list
 * is either untouched or contains bT itself (there is no retriangulation).
 *
 * @param x_center Physical quadrature point. This point is obtained by mapping a quadrature point
 * of the reference element to the triangle aT of the outer domain.
 * @param T Triangle of the inner integral.
 * @param mesh Mesh
 * @param reTriangle_list List of triangles which are obtained from the truncation.
 * @param is_placePointOnCap (Unused)
 * @return 0 if there is no interaction, -1 otherwise.
 */

int integrate_fractional_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                         const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                         double * termLocal, double * termNonloc,
                         double *termLocalPrime, double *termNonlocPrime, double * state_nodal_values,
                         double * adjoint_nodal_values);
// Helpers -------------------------------------------------------------------------------------------------------------
/**
 * @brief This truncation method checks whether the distance of the point point x_center (of aT) to
 * the bary center of bT is smaller than delta. The function writes the list of triangles
 * which are obtained from the trunaction into reTriangle_list. In this case the list
 * is either untouched or contains bT itself (there is no retriangulation).
 *
 * @param x_center Physical quadrature point. This point is obtained by mapping a quadrature point
 * of the reference element to the triangle aT of the outer domain.
 * @param T Triangle of the inner integral.
 * @param mesh Mesh
 * @param reTriangle_list List of triangles which are obtained from the truncation.
 * @param is_placePointOnCap (Unused)
 * @return 0 if there is no interaction, -1 otherwise.
*/

int method_baryCenter(const double * x_center, const ElementType & T, const MeshType & mesh, double * reTriangle_list, int is_placePointOnCap);
/**
 * @brief This truncation method retriangulates the triangle bT depending
 * on the distance of xCenter to bT w.r.t. the L2-norm ball. The function writes the list of triangles
 * which are obtained from the trunaction into reTriangle_list.
 *
 * @param x_center Physical quadrature point. This point is obtained by mapping a quadrature point
 * of the reference element to the triangle aT of the outer domain.
 * @param T Triangle of the inner integral.
 * @param mesh Mesh
 * @param reTriangle_list List of triangles which are obtained from the retriangulation.
 * @param is_placePointOnCap Switch for with caps integration.
 * @return 0 if there is no interaction, -1 otherwise.
 */
int method_retriangulate(const double * xCenter, const ElementType & T, const MeshType & mesh, double * reTriangleList, int isPlacePointOnCap);
/**
 * @brief This truncation method retriangulates the triangle bT depending
 * on the distance of xCenter to bT w.r.t. the L-infinity-norm ball. The function writes the list of triangles
 * which are obtained from the trunaction into reTriangle_list.
 *
 * @param x_center Physical quadrature point. This point is obtained by mapping a quadrature point
 * of the reference element to the triangle aT of the outer domain.
 * @param T Triangle of the inner integral.
 * @param mesh Mesh
 * @param reTriangle_list List of triangles which are obtained from the retriangulation.
 * @param is_placePointOnCap Switch for with caps integration.
 * @return 0 if there is no interaction, -1 otherwise.
 */
int method_retriangulateInfty(const double * xCenter, const ElementType & T, const MeshType & mesh, double * reTriangleList,
                              int isPlacePointOnCap);
/**
 * @brief This truncation method retriangulates the triangle bT depending
 * on the distance of xCenter to bT w.r.t. an exact quadrature rule for caps. The function writes the list of triangles
 * which are obtained from the trunaction into reTriangle_list, and additionally probides one point in the cap center
 * with the corresponding weight.
 *
 * @param x_center Physical quadrature point. This point is obtained by mapping a quadrature point
 * of the reference element to the triangle aT of the outer domain.
 * @param T Triangle of the inner integral.
 * @param mesh Mesh
 * @param reTriangle_list List of triangles which are obtained from the retriangulation.
 * @param capsList List of cap center points
 * @param capsWeights Corresponding weights of the quadrature
 * @param prtnCaps Number of caps.
 * @return
 */
int method_exact(const double * xCenter, const ElementType & T,
                 const MeshType & mesh, double * reTriangleList, double * capsList, double * capsWeights,
                 int * prtnCaps);
/**
 * @brief This truncation method returns the number of vertices of triangle T which interact
 * with x_center, i.e. have a l2 distance smaller than delta.
 *
 * @param x_center Physical quadrature point. This point is obtained by mapping a quadrature point
 * of the reference element to the triangle aT of the outer domain.
 * @param T Triangle of the inner integral.
 * @param mesh Mesh
 * @return 0 if there is no interaction. 1,2 or 3 otherwise.
 */
int method_subSuperSetBalls(const double * x_center, const ElementType & T, const MeshType & mesh);

int placePointOnCap(const double * y_predecessor, const double * y_current,
                    const double * x_center, double sqdelta, const double * TE,
                    const double * nu_a, const double * nu_b, const double * nu_c,
                    double orientation, int Rdx, double * R);
double placePointCapCenter(const double * y_predecessor, const double * y_current,
                           const double * x_center, const double sqdelta, const double * TE,
                           const double * nu_a, const double * nu_b, const double * nu_c,
                           const double orientation, double * capsList);
bool inTriangle(const double * y_new, const double * p, const double * q, const double * r,
                const double *  nu_a, const double * nu_b, const double * nu_c);

bool isFullyContained(const ElementType &aT, const ElementType &bT, const MeshType &mesh);

// Peridynamic Helper functions
void setupElement(const MeshType &mesh, const long * Vdx_new, ElementType &T);
int join(const ElementType &aT, const ElementType &bT, const MeshType &mesh,
         ElementType &aTsorted, ElementType &bTsorted, int * argSortA, int * argSortB);
double traffoCommonVertex0(double * alpha);
double traffoCommonVertex1(double * alpha);

double traffoCommonEdge0( double * alpha);
double traffoCommonEdge1( double * alpha);
double traffoCommonEdge2( double * alpha);
double traffoCommonEdge3( double * alpha);
double traffoCommonEdge4( double * alpha);

double traffoIdentical0( double * alpha);
double traffoIdentical1( double * alpha);
double traffoIdentical2( double * alpha);
double traffoIdentical3( double * alpha);
double traffoIdentical4( double * alpha);
double traffoIdentical5( double * alpha);

void scale(double * alpha);
void mirror(double * alpha);

#endif //NONLOCAL_ASSEMBLY_INTEGRATION_H

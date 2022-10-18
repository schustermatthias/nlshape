//
// Created by klar on 18.09.20.
//

#ifndef NONLOCAL_ASSEMBLY_MODEL_H
#define NONLOCAL_ASSEMBLY_MODEL_H

#include <cmath>
#include "MeshTypes.h"

// Pointer -------------------------------------------------------------------------------------------------------------
//// Kernel pointer and implementations
//TODO Use lambdas instead (check how to obtain inlining! Function pointer can never be inlined...)
//TODO change to template interface such that Kernel class with operator() which is a TwoPointFunction depending on Points
//TODO kernel is a class containing kernelFunction (TwoPointFunction) and truncation (TwoPointFunction, delta (TwoPointFunction)).
// The kernel knows its own delta, and this might be a function!
extern void (*model_kernel)(const double * x, long labelx, const double * y, long labely, const MeshType &mesh, double * kernel_val);
//TODO Default initialization of operator()...
void ERROR_wrongAccess(const double * x, long labelx, const double * y, long labely, const MeshType &mesh, double * kernel_val);

/**
 * @brief Constant kernel in 2D case. The constant is chosen such that the operator is equivalent to the laplacian for
 * polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_constant(const double * x, long labelx, const double * y, long labely, const MeshType &mesh,
                     double * kernel_val);
/**
 * @brief Constant kernel in 2D case with Linf truncation. The constant is chosen such that the operator is equivalent
 * to the laplacian for polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_constantLinf2D(const double *x, const long labelx, const double *y, const long labely, const MeshType &mesh,
                     double *kernel_val);

/**
 * @brief Constant kernel in 2D case which truncates itself.
 * This means that the kernel returns 0, if \f$ |x-y| > \delta \f$. It is used to simulate a specific ball approximation.
 * The constant is chosen such that the operator is equivalent to the laplacian for
 * polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_constantTruncated(const double *x, const long labelx, const double *y, const long labely, const MeshType &mesh,
                              double *kernel_val);
/**
  * @brief This kernel is defined by \f$ \delta^2 - \|z\|^2 \f$. It is radial but not
  * singular or constant. The second moment is given by \f$ \pi \delta ^6 / 12 \f$.
  * @param x
  * @param labelx
  * @param y
  * @param labely
  * @param sqdelta
  * @param kernel_val
  */
void kernel_parabola(const double *x, const long labelx, const double *y, const long labely, const MeshType &mesh,
                     double *kernel_val);
/**
 * @brief Constant kernel in 1D case. The constant is chosen such that the operator is equivalent to the laplacian for
 * polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_constant1D(const double *x, const long labelx, const double *y, const long labely, const MeshType &mesh,
                       double *kernel_val);
/**
 * @brief Antisymmetric kernel in the 1D case.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_antisymmetric1D(const double *x, const long labelx, const double *y, const long labely, const MeshType &mesh,
                       double *kernel_val);
/**
 * @brief Constant kernel in 3D case. The constant is chosen such that the operator is equivalent to the laplacian for
 * polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_constant3D(const double * x, long labelx, const double * y, long labely, const MeshType &mesh,
                       double * kernel_val);
/**
 * @brief Kernel depending on the triangle labels. Can be used to model nonlocal to nonlocal coupling.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_labeled(const double * x, long labelx, const double * y, long labely, const MeshType &mesh,
                    double * kernel_val);

void kernel_integrable_sym(const double * x, long labelx, const double * y, long labely, const MeshType &mesh,
                           double * kernel_val);
                    	
void kernel_integrable_unsym(const double * x, long labelx, const double * y, long labely, const MeshType &mesh,
                             double * kernel_val);

void kernel_notch(const double * x, const long labelx, const double * y, const long labely, const MeshType &mesh,
                  double * kernel_val);
void kernel_labeledNotch(const double * x, const long labelx, const double * y, const long labely, const MeshType &mesh,
                         double * kernel_val);
void kernel_labeledValve(const double * x, const long labelx, const double * y, const long labely, const MeshType &mesh,
                         double * kernel_val);
/**
 * @brief Kernel for peridynamics diffusion model. The scalar valued weakly singular kernel reads as
 *
 *  \f[
 * \gamma(x,y) = \frac{1}{\| x - y \|}  \frac{3}{\pi \delta ^3}.
 *  \f]
 * The constant is chosen such that the operator is equivalent to the laplacian for
 * polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_linearPrototypeMicroelastic(const double * x, long labelx, const double * y, long labely,
                                           const MeshType &mesh, double * kernel_val);

/**
 * @brief Kernel for fractional Laplacian. The scalar valued weakly singular kernel reads as
 *
 *  \f[
 * \gamma(x,y) = \frac{1}{\| x - y \|}  \frac{3}{\pi \delta ^3}.
 *  \f]
 * The constant is chosen such that the operator is equivalent to the laplacian for
 * polynomials of degree less or equal to 2.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernel_fractional(const double * x, const long labelx, const double * y, const long labely,
                         const MeshType &mesh, double * kernel_val);
/**
 * @brief Kernel for the peridynamics model. The matrix valued weakly singular kernel reads as
 *
 *  \f[
 * \gamma(x,y) = (x-y) \otimes (x-y) \frac{1}{\| x - y \|^3}  \frac{12}{\pi \delta ^3}.
 *  \f]
 *
 * The constant is chosen such that the operator is equivalent to linear elasticity for
 * polynomials of degree less or equal to 2 and lame paramter
 * \f$\mu = \lambda = \frac{\pi}{4}\f$.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */

void kernel_fractional_shape(const double * x, const long labelx, const double * y, const long labely,
                             const MeshType &mesh, double * kernel_val);

void kernelField_linearPrototypeMicroelastic(const double * x, long labelx, const double * y, long labely,
                                                const MeshType &mesh, double * kernel_val);
/**
 * @brief Kernel for the peridynamics model. The matrix valued weakly singular kernel reads as
 *
 *  \f[
 * \gamma(x,y) = (x-y) \otimes (x-y) \frac{1}{\| x - y \|^3}  \frac{12}{\pi \delta ^3}.
 *  \f]
 *
 * The constant is chosen such that the operator is equivalent to linear elasticity for
 * polynomials of degree less or equal to 2 and lame paramter
 * \f$\mu = \lambda = \frac{\pi}{4}\f$.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernelField_linearPrototypeMicroelastic3D(const double * x, long labelx, const double * y, long labely,
                                               const MeshType &mesh, double * kernel_val);
/**
 * @brief Constant matrix valued kernel. For testing. The constant is the same as in the scalar case.
 *
 * @param x Physical point of the outer integration region.
 * @param labelx Label of the outer triangle.
 * @param y Physical point of the inner integration region.
 * @param labely Label of inner triangle.
 * @param mesh Mesh
 * @param kernel_val Value of the the kernel. Pointer to double in case of diffusion. Pointer to a array
 * of shape d x d in case of peridynamics.
 */
void kernelField_constant(const double * x, long labelx, const double * y, long labely,
                             const MeshType &mesh, double * kernel_val);


extern void (*model_f)(const double * x, double * forcing_out);
void ERROR_wrongAccess(const double * x, double * forcing_out);
/**
 * @brief Constant forcing function f = 1.
 * @param x
 * @param forcing_out
 */
void f_constant(const double * x, double * forcing_out);
/**
 * @brief Linear vector valued forcing function \f$f = (1 + 2 x_1, x_2 ) \frac{\pi}{2} \f$. This
 * function is used for checking rates. The corresponding solution reads as
 * \f[
 *      u(x) = (x_2^2, x_1^2 x_2) \frac{2}{5}.
 * \f]
 * @param x
 * @param forcing_out
 */
void fField_linear(const double * x, double * forcing_out);
/**
 * @brief Constant vector valued forcing function.
 * @param x
 * @param forcing_out
 */
void fField_linear3D(const double * x, double * forcing_out);
/**
 * @brief Constant vector valued forcing function.
 * @param x
 * @param forcing_out
 */
void fField_constantRight(const double * x, double * forcing_out);
/**
 * @brief Constant vector valued forcing function.
 * @param x
 * @param forcing_out
 */
void fField_constantDown(const double * x, double * forcing_out);
/**
 * @brief Constant vector valued forcing function. Used for rates checks.
 * @param x
 * @param forcing_out
 */
void fField_constantBoth(const double * x, double * forcing_out);
/**
 * @brief Linear scalar valued forcing function \f$f = -2 (x_1 + 1)\f$. Used for computation of rates.
 * @param x
 * @param forcing_out
 */

void f_linear(const double * x, double * forcing_out);
/**
 * @brief Linear scalar valued forcing function \f$f = -2 (x_1+ 1)\f$. Used for computation of rates.
 * @param x
 * @param forcing_out
 */
void f_gaussian(const double * x, double * forcing_out);
/**
 * @brief Right hand side \f$f = - 4 (\| x\|^2 - 1)e^{\|x\|^2}\f$.
 * @param x
 * @param forcing_out
 */
void f_jump(const double * x, double * forcing_out);
/**
 * @brief Right hand side \f$f = - 2 \chi_{R^+}(x_1 - x_2)\f$.
 * @param x
 * @param forcing_out
 */
void f_tensorsin(const double * x, double * forcing_out);
/**
 * @brief Right hand side \f$f = 32 \pi^2 sin(4 \pi x)sin(4 \pi y )\f$.
 * @param x
 * @param forcing_out
 */
void f_linear3D(const double * x, double * forcing_out);
void f_linear1D(const double * x, double * forcing_out);

#endif //NONLOCAL_ASSEMBLY_MODEL_H

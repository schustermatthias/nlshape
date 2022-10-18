//
// Created by klar on 18.09.20.
//

#ifndef NONLOCAL_ASSEMBLY_MATHHELPERS_H
#define NONLOCAL_ASSEMBLY_MATHHELPERS_H

#include "MeshTypes.h"

// ___ MATH HELERPS DECLARATION ________________________________________________________________________________________
const double EPSILON=1e-15;
const double EPSILON_CHKS=1e-8;

// Miscellaneous helpers ###############################################################################################
double evaluateZeta(const long * ptr_indices,
                    const long * ptr_indptr,
                    long nColumns,
                    long aTdx, long bTdx);
void solve2x2(const double *, const double *, double *);                 // Solve 2x2 System with LU
void rightNormal(const double * y0, const double * y1, double orientation, double * normal);
int faculty(int n);
/**
 * @brief Rescales a cube \f$ [-1, 1]^4 \f$ to a cube
 * \f$ [0, 1]^4 \f$. The determinant of this
 * transformation is given by the constant SCALEDET.
 *
 * @param alpha, a single input point of dimension 4.
 */
void scale(double * alpha);
/**
 * @brief Maps a tuple of two triangles \f$\left\lbrace z | z_2 \in
 * [0,1], z_1 \ in [0,z_2] \right\rbrace\f$
 * to a tuple of two standard simplices
 * \f$\left\lbrace z | z_1 \in [0,1], z_2 \ in [0,1-z_1] \right\rbrace\f$.
 * The determinant of this transformation is 1.
 *
 * @param alpha, a single input point of dimension 4.
 */
void mirror(double * alpha);

// Matrix operations ###################################################################################################
// Double
double absDet(const double * E);                                         // Compute determinant
double absDet(const double * E, int dim);
double signDet(const double * E);
double signDet(const double * E, const MeshType & mesh);
void baryCenter(const double * E, double * bary);                        // Bary Center
void baryCenter(int dim, const double * E, double * bary);
void toRef(const double * E, const double * phys_x, double * ref_p);     // Pull point to Reference Element (performs 2x2 Solve)
void toPhys(const double * E, const double * p, double * out_x);         // Push point to Physical Element
void toPhys(const double * E, const double * p, int dim, double * out_x);
void get_div(const double * E, double * div);

// Vector operations ###################################################################################################

// Double
double vec_sqL2dist(const double * x, const double * y, int len);      // L2 Distance
double vec_LInfdist(const double * x, const double * y, const int len);
double vec_dot(const double * x, const double * y, int len);           // Scalar Product
double vec_sum(const double *x, int len);                               // Vector Sum
int doubleVec_any(const double * vec, int len);                        // Any
void doubleVec_tozero(double *, int);               // Reset to zero
void doubleVec_subtract(const double * vec1, const double * vec2, double * out, int len);
void doubleVec_midpoint(const double * vec1, const double * vec2, double * midpoint, int len);
void doubleVec_scale(double lambda, const double * vec, double * out, int len);
void doubleVec_add(const double * vec1, const double * vec2, double * out, int len);
void doubleVec_copyTo(const double * input, double * output, int len);
// Long
int longVec_all(const long *, int);                // All
int longVec_any(const long *, int);                // Any

// Int
void intVec_tozero(int *, int);                    // Reset to zero

// Scalar operations ###################################################################################################
double absolute(double);                                  // Get absolute value
bool double_eq(double x, double y, double eps=EPSILON);     // Compare to double values
double scal_sqL2dist(double x, double y);           // L2 Distance



#endif //NONLOCAL_ASSEMBLY_MATHHELPERS_H

/**
    Integration and retriangulation routines.
  All integration routines share a common signature.
The central difference between different routines is the way they handle the truncation of the domain
 of integration. In almost all cases, the truncation of the inner triangle bT
is performed for each quadrature point in the outer
triangle aT.
In addition for singular kernels special care is needed. Due to the necessity of
transformations of the domain of integration those integration routines are tied to a kernel.

The integration routine changes the data in *termLocal and *termNonloc to. **The arrays have to be zero-initialized.**

    * termLocal = int_aT phiA(x) phiB(x) int_bT ker(x,y) dy dx\n,
    * termNonloc = int_aT phiA(x) int_bT phiB(y) ker(y,x) dy dx.

Please note that the nonlocal term has to be subtracted, while the local term has to be added to the stiffness
matrix.

@param aT    Triangle of the outer integral.
@param bT    Triangle of the inner integral.
@param quadRule Quadrature rule.
@param mesh  Mesh.
@param conf  Confuration.
@param is_firstbfslayer Switch to tell whether the integration is happening in the first layer of the breadth first
search. This variable is true only if the kernel is singular. In that case the integrals between aT and its immediate
neighbours have to be handled with special care.
@param termLocal This term contains the local part of the integral
@param termNonloc This term contains the nonlocal part of the integral

@file integration.cpp
@author Manuel Klar
@version 0.1 25/08/20
**/

#ifndef NONLOCAL_ASSEMBLY_INTEGRATION_CPP
#define NONLOCAL_ASSEMBLY_INTEGRATION_CPP

#include <list>
#include "MeshTypes.h"
#include "mathhelpers.h"
#include "model.h"
#include "integration.h"

const double SCALEDET = 0.0625;
//int (*integrate)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
//                         const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
//                         double *termLocalPrime, double *termNonlocPrime);
int (*integrate_remote)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                  const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                  double *termLocalPrime, double *termNonlocPrime);
int (*integrate_close)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                  const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                  double *termLocalPrime, double *termNonlocPrime);
int (*method)(const double * xCenter, const ElementType & T, const MeshType & mesh, double * reTriangleList,
                int isPlacePointOnCap);

int (*integrate_remote_shape)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                        const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                        double *termLocalPrime, double *termNonlocPrime, double state_nodal_values[6], double adjoint_nodal_values[6]);
int (*integrate_close_shape)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                       const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                       double *termLocalPrime, double *termNonlocPrime, double state_nodal_values[6], double adjoint_nodal_values[6]);


int ERROR_wrongAccess(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                      const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                      double *termLocalPrime, double *termNonlocPrime){
    cout << "ERROR in integration.cpp/ERROR_wrongAccess(): You chose no integration method, but the assembly routine tried to call it." << endl;
    cout << "The names of possible choices are given in lookup_configuration() in Cassemble.cpp." << endl;
    abort();
    return -1;
};
int ERROR_wrongAccess_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                      const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                      double *termLocalPrime, double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values){
    cout << "ERROR in integration.cpp/ERROR_wrongAccess(): You chose no integration method, but the assembly routine tried to call it." << endl;
    cout << "The names of possible choices are given in lookup_configuration() in Cassemble.cpp." << endl;
    abort();
    return -1;
};
// ___ INTEGRATION IMPLEMENTATION ______________________________________________________________________________________
// TODO Return int = 1, if interaction, 0 otherwise. This helps to set other values here if required
// TODO integrate should take elements, two-points-functions, and quadrule and return
int integrate(const ElementType &aT, const ElementType &bT,
                                                         const QuadratureType &quadRule,
                                                         const MeshType &mesh, const ConfigurationType &conf,
                                                         bool is_firstbfslayer, double *termLocal,
                                                         double *termNonloc, double *termLocalPrime, double *termNonlocPrime){
    if (is_firstbfslayer) {
        return integrate_close(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal,
                                                          termNonloc, termLocalPrime, termNonlocPrime);
    } else {
        return integrate_remote(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal, termNonloc,
                                termLocalPrime, termNonlocPrime);
    }
}

int integrate_shape(const ElementType &aT, const ElementType &bT,
              const QuadratureType &quadRule,
              const MeshType &mesh, const ConfigurationType &conf,
              bool is_firstbfslayer, double *termLocal,
              double *termNonloc, double *termLocalPrime, double *termNonlocPrime,
              double * state_nodal_values, double * adjoint_nodal_values){
    if (is_firstbfslayer) {
        return integrate_close_shape(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal,
                               termNonloc, termLocalPrime, termNonlocPrime, state_nodal_values, adjoint_nodal_values);
    } else {
        return integrate_remote_shape(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal, termNonloc,
                                      termLocalPrime, termNonlocPrime, state_nodal_values, adjoint_nodal_values);
    }
}

// Integration Methods #################################################################################################
int integrate_weakSingular(const ElementType &aT, const ElementType &bT,
                           const QuadratureType &quadRule,
                           const MeshType &mesh, const ConfigurationType &conf,
                           bool is_firstbfslayer,
                           double * termLocal, double * termNonloc,
                           double *termLocalPrime, double *termNonlocPrime){

    ElementStruct aTsorted, bTsorted;
    int argSortA[3], argSortB[3];
    int nTraffos;

    int nEqual = join(aT, bT, mesh, aTsorted, bTsorted, argSortA, argSortB);
    std::list<double (*)(double *)> traffoList;

    double factors;
    //const int dim = mesh.dim;
    //double alpha[4], alphaCanceled[4]; // traffodet
    const arma::Mat<double> * traffodetCanceled;
    const arma::Cube<double> * alpha;
    const arma::Cube<double> * alphaCanceled;

    //double traffodetCanceled;//, x[2], y[2],
    double x_canceled[2], y_canceled[2];//, kernel_val=0.;

    double kernel_val[mesh.outdim*mesh.outdim];
    double psix[3], psiy[3];
    //cout << "Tensor Gauss" << endl;
    if (nEqual == 1){
        traffoList = traffoCommonVertex;
        alpha = &quadRule.alpha_CommonVertex;
        alphaCanceled = &quadRule.alphaCanceled_CommonVertex;
        traffodetCanceled = &quadRule.traffodetWeakCanceled_CommonVertex;
        nTraffos = 2;
    } else if (nEqual == 2){
        traffoList = traffoCommonEdge;
        alpha = &quadRule.alpha_CommonEdge;
        alphaCanceled = &quadRule.alphaCanceled_CommonEdge;
        traffodetCanceled = &quadRule.traffodetWeakCanceled_CommonEdge;
        nTraffos = 5;
    } else if (nEqual == 3) {
        traffoList = traffoIdentical;
        alpha = &quadRule.alpha_Identical;
        alphaCanceled = &quadRule.alphaCanceled_Identical;
        traffodetCanceled = &quadRule.traffodetWeakCanceled_Identical;
        nTraffos = 6;
    } else {
        printf("Error in integrate_weakSingular: nEqual == %i", nEqual );
        printf( "Error in integrate_weakSingular: This should not have happened." );
        abort();
    }

    //TODO Consider untying retriangulation and integral evaluation..(?)
    //The code below is much simpler that the retriangulation code..
    for(int traffoCounter=0; traffoCounter < nTraffos; traffoCounter++) {
        //if (aT.Tdx == 0) cout << "traffo " << traffoCounter << endl;
        for (int k = 0; k < quadRule.nPg; k++) {
            double traffoWeakCanceled =   (*traffodetCanceled)(k, traffoCounter);

            toPhys(aTsorted.E, &(*alphaCanceled)(0, k, traffoCounter), 2, x_canceled);
            toPhys(bTsorted.E, &(*alphaCanceled)(2, k, traffoCounter), 2, y_canceled);

            // Eval Kernel(x-y)
            model_kernel(x_canceled, aTsorted.label, y_canceled, bTsorted.label, mesh, kernel_val);


            model_basisFunction(&(*alpha)(0, k, traffoCounter), 2, psix);
            model_basisFunction(&(*alpha)(2, k, traffoCounter), 2, psiy);

            for (int a = 0; a < mesh.dVertex; a++) {
                for (int aOut = 0; aOut < mesh.outdim; aOut++) {
                    for (int b = 0; b < mesh.dVertex; b++) {
                        for (int bOut = 0; bOut < mesh.outdim; bOut++) {

                            factors =
                                    kernel_val[mesh.outdim * aOut + bOut] * traffoWeakCanceled *
                                    quadRule.dg[k] * aTsorted.absDet * bTsorted.absDet;

                            termLocal[mesh.outdim * mesh.dVertex * (argSortA[a]*mesh.outdim + aOut) +
                                      argSortA[b]*mesh.outdim + bOut] += factors * psix[a] * psix[b];

                            termLocalPrime[mesh.outdim * mesh.dVertex * (argSortB[a]*mesh.outdim + aOut) +
                                           argSortB[b]*mesh.outdim + bOut] += factors * psiy[a] * psiy[b];

                            termNonloc[mesh.outdim * mesh.dVertex * (argSortA[a]*mesh.outdim + aOut)
                                       + argSortB[b]*mesh.outdim + bOut] += factors * psix[a] * psiy[b];

                            termNonlocPrime[mesh.outdim * mesh.dVertex * (argSortA[a]*mesh.outdim + aOut)
                                            + argSortB[b]*mesh.outdim + bOut] += factors * psix[a] * psiy[b];
                            //TODO Check whether this is mesh.outdim * mesh.dVertex * (argSortA[b]*mesh.outdim + aOut)
                            //              + argSortB[a- shift]*mesh.outdim + bOut
                            //termNonloc[mesh.dVertex * a + b] += 2*factors * psix[a] * psiy[b];
                            //cout << termLocal[mesh.outdim*mesh.dVertex * a + b] << endl;
                        }
                    }
                }
                //cout << endl;
            }
            //cout << "Thank you!" << endl;
            //abort();
        }
        //traffoCounter++;
    }
    return 1;
}

int integrate_fractional(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                           const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                           double * termLocal, double * termNonloc, double *termLocalPrime, double *termNonlocPrime){
        //printf("Ok, is_firstbfslayer!");
        //abort();
        ElementStruct aTsorted, bTsorted;
        int argSortA[3], argSortB[3];
        const int nEqual = join(aT, bT, mesh, aTsorted, bTsorted, argSortA, argSortB);
        double factors;
        double x_canceled[2], y_canceled[2];
        const arma::Mat<double> * traffodetCanceled;
        const arma::Cube<double> * alphaCanceled;
        int nTraffos;
        double kernel_val[mesh.outdim*mesh.outdim];
        const int nLocalVerts=2*mesh.dVertex-nEqual;
        double psix[3], psiy[3], psixy[nLocalVerts];

        //cout << "Tensor Gauss" << endl;
        if (nEqual == 1){
            alphaCanceled = &quadRule.alphaCanceled_CommonVertex;
            traffodetCanceled = &quadRule.traffodetFractionalCanceled_CommonVertex;
            nTraffos = 2;
        } else if (nEqual == 2){
            alphaCanceled = &quadRule.alphaCanceled_CommonEdge;
            traffodetCanceled = &quadRule.traffodetFractionalCanceled_CommonEdge;
            nTraffos = 5;
        } else if (nEqual == 3) {
            alphaCanceled = &quadRule.alphaCanceled_Identical;
            traffodetCanceled = &quadRule.traffodetFractionalCanceled_Identical;
            nTraffos = 6;
        } else {
            printf("Error in integrate_fractional: nEqual == %i", nEqual );
            printf( "Error in integrate_fractional: This should not have happened." );
            abort();
        }

        for(int traffoCounter=0; traffoCounter < nTraffos; traffoCounter++) {
            //if (aT.Tdx == 0) cout << "Traffo Number " << traffoCounter << ": \n";
            for (int k = 0; k < quadRule.nPg; k++) {
                // Depending on the defree of the singularity we have to cancel out some parts of the
                // determinant. This happens here, based on the correct power of the variable xi = alpha[0]
                double traffoFractionalCanceled =  (*traffodetCanceled)(k, traffoCounter);
                //if (aT.Tdx == 0) cout << "      traffodetCanceled : " <<  ",\t" << traffoFractionalCanceled << endl;

                toPhys(aTsorted.E, &(*alphaCanceled)(0, k, traffoCounter), 2, x_canceled);
                toPhys(bTsorted.E, &(*alphaCanceled)(2, k, traffoCounter), 2, y_canceled);

                // Eval Kernel(x-y)
                model_kernel(x_canceled, aTsorted.label, y_canceled, bTsorted.label, mesh, kernel_val);

                model_basisFunction(&(*alphaCanceled)(0, k, traffoCounter), 2, psix);
                model_basisFunction(&(*alphaCanceled)(2, k, traffoCounter), 2, psiy);

                int shift = mesh.dVertex-nEqual;

                for(int a = 0; a < nEqual; a++){
                    psixy[a] = psix[a] - psiy[a];
                }
                for(int a = nEqual; a < mesh.dVertex; a++){
                    psixy[a] = psix[a];
                    psixy[a + shift] = - psiy[a];
                }

                for (int a = 0; a < nLocalVerts; a++){
                    for (int aOut = 0; aOut < mesh.outdim; aOut++) {
                        //cout << "----------------------\n" << a << endl;
                        //cout << "a" << a << endl;
                        for (int b = 0; b < nLocalVerts; b++) {

                            for (int bOut = 0; bOut < mesh.outdim; bOut++) {
                                //cout <<  "b" <<b << endl;
                                factors = kernel_val[mesh.outdim * aOut + bOut] * traffoFractionalCanceled *
                                          psixy[a] * psixy[b] *
                                          quadRule.dg[k] * aTsorted.absDet * bTsorted.absDet;
                                if (a < mesh.dVertex) {
                                    if (b < mesh.dVertex) {
                                        termLocal[mesh.outdim * mesh.dVertex * (argSortA[a]*mesh.outdim + aOut) +
                                                  argSortA[b]*mesh.outdim + bOut] += factors;

                                    } else {
                                        termNonloc[mesh.outdim * mesh.dVertex * (argSortA[a]*mesh.outdim + aOut)
                                                   + argSortB[b- shift]*mesh.outdim + bOut] -= factors;
                                    }
                                } else {
                                    if (b < mesh.dVertex) {
                                        termNonlocPrime[mesh.outdim * mesh.dVertex * (argSortA[b]*mesh.outdim + aOut)
                                                        + argSortB[a- shift]*mesh.outdim + bOut] -= factors;

                                    } else {
                                        termLocalPrime[mesh.outdim * mesh.dVertex * (argSortB[a- shift]*mesh.outdim + aOut) +
                                                       argSortB[b- shift]*mesh.outdim + bOut] += factors;

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return 1;
}

int integrate_fractional_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                               const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                               double * termLocal, double * termNonloc, double *termLocalPrime, double *termNonlocPrime,
                               double * state_nodal_values, double * adjoint_nodal_values){
    //printf("Ok, is_firstbfslayer!");
    //abort();
    ElementStruct aTsorted, bTsorted;
    int argSortA[3], argSortB[3];
    const int nEqual = join(aT, bT, mesh, aTsorted, bTsorted, argSortA, argSortB);
    double factors;
    double x_canceled[2], y_canceled[2];
    const arma::Mat<double> * traffodetCanceled;
    const arma::Cube<double> * alphaCanceled;
    int nTraffos;
    double kernel_val[mesh.outdim*mesh.outdim];
    const int nLocalVerts=2*mesh.dVertex-nEqual;
    double psix[3], psiy[3], psixy[nLocalVerts];
    double state_values[2] = {0.0, 0.0};
    double adjoint_values[2] = {0.0, 0.0};

    //cout << "Tensor Gauss" << endl;
    if (nEqual == 1){
        alphaCanceled = &quadRule.alphaCanceled_CommonVertex;
        traffodetCanceled = &quadRule.traffodetFractionalCanceled_CommonVertex;
        nTraffos = 2;
    } else if (nEqual == 2){
        alphaCanceled = &quadRule.alphaCanceled_CommonEdge;
        traffodetCanceled = &quadRule.traffodetFractionalCanceled_CommonEdge;
        nTraffos = 5;
    } else if (nEqual == 3) {
        alphaCanceled = &quadRule.alphaCanceled_Identical;
        traffodetCanceled = &quadRule.traffodetFractionalCanceled_Identical;
        nTraffos = 6;
    } else {
        printf("Error in integrate_fractional: nEqual == %i", nEqual );
        printf( "Error in integrate_fractional: This should not have happened." );
        abort();
    }

    for(int traffoCounter=0; traffoCounter < nTraffos; traffoCounter++) {
        //if (aT.Tdx == 0) cout << "Traffo Number " << traffoCounter << ": \n";
        for (int k = 0; k < quadRule.nPg; k++) {
            // Depending on the defree of the singularity we have to cancel out some parts of the
            // determinant. This happens here, based on the correct power of the variable xi = alpha[0]
            double traffoFractionalCanceled =  (*traffodetCanceled)(k, traffoCounter);
            //if (aT.Tdx == 0) cout << "      traffodetCanceled : " <<  ",\t" << traffoFractionalCanceled << endl;

            toPhys(aTsorted.E, &(*alphaCanceled)(0, k, traffoCounter), 2, x_canceled);
            toPhys(bTsorted.E, &(*alphaCanceled)(2, k, traffoCounter), 2, y_canceled);

            // Eval Kernel(x-y)
            model_kernel(x_canceled, aTsorted.label, y_canceled, bTsorted.label, mesh, kernel_val);

            model_basisFunction(&(*alphaCanceled)(0, k, traffoCounter), 2, psix);
            model_basisFunction(&(*alphaCanceled)(2, k, traffoCounter), 2, psiy);

            // int shift = mesh.dVertex-nEqual;

            for(int a = 0; a < nEqual; a++){
                psixy[a] = psix[a] - psiy[a];
            }
            for(int a = nEqual; a < mesh.dVertex; a++){
                psixy[a] = psix[a];
                //psixy[a + shift] = - psiy[a];
            }
            state_values[0] = 0.0;
            state_values[1] = 0.0;
            adjoint_values[0] = 0.0;
            adjoint_values[1] = 0.0;

            for (int a = 0; a < mesh.dVertex; a++){
                state_values[0] += state_nodal_values[argSortA[a]] * psix[a];
                adjoint_values[0] += adjoint_nodal_values[argSortA[a]] * psix[a];
                state_values[1] += state_nodal_values[argSortB[a] + mesh.dVertex] * psiy[a];
                adjoint_values[1] += adjoint_nodal_values[argSortB[a] + mesh.dVertex] * psiy[a];
            }

            factors = kernel_val[0] * traffoFractionalCanceled * (state_values[0] - state_values[1]) *
                      (adjoint_values[0] - adjoint_values[1]) * quadRule.dg[k] * aTsorted.absDet * bTsorted.absDet;
            termLocal[0] += factors;

            for (int a = 0; a < mesh.dVertex; a++){
                termNonloc[argSortA[a]] += factors * psixy[a]
                        * (x_canceled[0] - y_canceled[0])/(vec_sqL2dist(x_canceled, y_canceled, 2));
                termNonlocPrime[argSortA[a]] += factors * psixy[a]
                        * (x_canceled[1] - y_canceled[1])/(vec_sqL2dist(x_canceled, y_canceled, 2));
            }
        }
    }
    return 1;
}

int integrate_fullyContained(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                              const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                              double *termLocal, double *termNonloc, double *termLocalPrime, double *termNonlocPrime){

    const int dim = mesh.dim;
    double x[2], y[2];//, kernel_val=0.;
    double kernel_val[mesh.outdim*mesh.outdim];
    double kernelPrime_val[mesh.outdim*mesh.outdim];

    for (int k = 0; k < quadRule.nPx; k++) {
        toPhys(aT.E, &(quadRule.Px[dim * k]), dim, x);
        for (int i = 0; i < quadRule.nPy; i++) {
            toPhys(bT.E, &(quadRule.Py[dim * i]), dim, y);

            // Eval Kernel(x-y)
            model_kernel(x, aT.label, y, bT.label, mesh, kernel_val);
            model_kernel(y, bT.label, x, aT.label, mesh, kernelPrime_val);


            for (int a = 0; a < mesh.dVertex; a++) {
                for (int aOut = 0; aOut < mesh.outdim; aOut++) {
                    for (int b = 0; b < mesh.dVertex; b++) {
                        for (int bOut = 0; bOut < mesh.outdim; bOut++) {
                            //outTest[mesh.outdim*mesh.dVertex * a + b] = kernel_val[mesh.outdim * (a%mesh.outdim) + b%mesh.outdim]
                            //        + psitest[a/mesh.outdim]*psitest[b/mesh.outdim];
                            //cout << outTest[mesh.outdim*mesh.dVertex * a + b] << ",   ";
                            double factors = quadRule.dx[k] * quadRule.dy[i] * aT.absDet * bT.absDet;
                            //factors = kernel_val * traffodet * scaledet * quadule.dg[k] * aTsorted.absDet * bTsorted.absDet;
                            termLocal[mesh.outdim * mesh.dVertex * (a*mesh.outdim + aOut) + b*mesh.outdim + bOut]
                                += factors * quadRule.psix(a, k) * quadRule.psix(b, k)
                                        * kernel_val[mesh.outdim * aOut + bOut];
                            termLocalPrime[mesh.outdim * mesh.dVertex * (a*mesh.outdim + aOut) + b*mesh.outdim + bOut]
                                += factors * quadRule.psiy(a, i) * quadRule.psiy(b, i)
                                        * kernelPrime_val[mesh.outdim * aOut + bOut];
                            //termLocal[mesh.dVertex * a + b] += 2*factors* quadRule.psix[a] * quadRule.psix[b];
                            // [10] siehe [9]
                            termNonloc[mesh.outdim * mesh.dVertex * (a*mesh.outdim + aOut) + b*mesh.outdim + bOut]
                                += factors * quadRule.psix(a, k) * quadRule.psiy(b, i)
                                        * kernelPrime_val[mesh.outdim * aOut + bOut];
                            termNonlocPrime[mesh.outdim * mesh.dVertex * (a*mesh.outdim + aOut) + b*mesh.outdim + bOut]
                                += factors * quadRule.psix(a, k) * quadRule.psiy(b, i)
                                        * kernel_val[mesh.outdim * aOut + bOut];
                            //TODO Check corresponding Todo in weak-singular integration.
                        }
                    }
                }
                //cout << endl;
            }
            //cout << "Thank you!" << endl;
            //abort()
        }
    }
    return 1;
}

int integrate_retriangulate(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                             const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                             double *termNonloc, double *termLocalPrime, double *termNonlocPrime) {
    // One would assume, the code became faster by the check below - but until delta = .3 it still didn't
    // Code can be more efficient, if we use very simple rules for elements which fully lie in the interaction set.
    //if ((mesh.maxDiameter > EPSILON) && (mesh.delta - 2*mesh.maxDiameter > 0) && isFullyContained(aT, bT, mesh)){
    //    return integrate_fullyContained(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal, termNonloc,
    //                             termLocalPrime, termNonlocPrime);
    //}

    const int dim = mesh.dim;
    int k = 0, a = 0, b = 0;
    double x[dim];
    // [x 11] [mesh.outdim*mesh.outdim]
    //double innerLocal = 0;
    double innerLocal[mesh.outdim*mesh.outdim];
    //double innerLocalPrime[mesh.outdim*mesh.outdim];

    // [x 12] [mesh.outdim*mesh.outdim*mesh.dVertex]
    // double innerNonloc[mesh.dVertex];
    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];
    double innerNonlocPrime[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0, rTdx = 0, Rdx = 0;
    // [x 13] kernel_val [mesh.outdim*mesh.outdim]
    double kernel_val[mesh.outdim*mesh.outdim], kernelPrime_val[mesh.outdim*mesh.outdim];;

    double rTdet = 0;
    double physical_quad[dim];
    double reference_quad[dim];
    double psi_value[mesh.dVertex];
    //double psi_value_test[3] = {20., 30., 50.};
    double reTriangle_list[36 * mesh.dVertex * dim];
    doubleVec_tozero(reTriangle_list, 36 * mesh.dVertex * dim);
    int doesInteract=0;

    //[DEBUG]
    //printf("\nouterInt_full----------------------------------------\n");
    for (k = 0; k < quadRule.nPx; k++) {
        //printf("k %i, quadRule.nPx %i\n", k, quadRule.nPx);
        toPhys(aT.E, &(quadRule.Px[dim * k]), mesh.dim, x);
        //printf("\nInner Integral, Iterate %i\n", k);
        //printf("\Physical x [%17.16e, %17.16e]\n",  x[0], x[1]);
        //innerInt_retriangulate(x, aT, bT, quadRule, sqdelta, &innerLocal, innerNonloc);

        //is_placePointOnCap = true;
        Rdx = method(x, bT, mesh, reTriangle_list, conf.is_placePointOnCap); // innerInt_retriangulate
        doesInteract += Rdx;
        //Rdx = baryCenterMethod(x, bT, mesh, reTriangle_list, is_placePointOnCap);
        //Rdx = quadRule.interactionMethod(x, bT, mesh, reTriangle_list);

        //[DEBUG]
        //printf("Retriangulation Rdx %i\n", Rdx);
        //for (i=0;i<Rdx;i++){
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i], reTriangle_list[2 * 3 * i+1]);
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i+2], reTriangle_list[2 * 3 * i+3]);
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i+4], reTriangle_list[2 * 3 * i+5]);
        //printf("absDet %17.16e\n", absDet(&reTriangle_list[2 * 3 * i]));
        //}
        // [x 14] doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);
        //innerLocal = 0.0;
        doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);

        // [x 15] doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim*mesh.dVertex);
        // doubleVec_tozero(innerNonloc, mesh.dVertex);
        doubleVec_tozero(innerNonloc, mesh.outdim*mesh.outdim*mesh.dVertex);
        doubleVec_tozero(innerNonlocPrime, mesh.outdim*mesh.outdim*mesh.dVertex);

        if (Rdx == 0) {
        } else {
            //printf("\nInner Integral\n");
            for (rTdx = 0; rTdx < Rdx; rTdx++) {
                //printf("rTdx %i \n",rTdx);
                for (i = 0; i < quadRule.nPy; i++) {
                    // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                    toPhys(&reTriangle_list[dim * mesh.dVertex * rTdx], &(quadRule.Py[dim * i]), physical_quad);
                    // Determinant of Triangle of retriangulation
                    rTdet = absDet(&reTriangle_list[dim * mesh.dVertex * rTdx]);
                    // inner Local integral with ker
                    model_kernel(x, aT.label, physical_quad, bT.label, mesh, kernel_val);
                    // [x 16]
                    // INNER LOCAL ORDER [(0,0), (0,1), (1,0), (1,1)] = KERNEL ORDER
                    for (int o=0; o<mesh.outdim*mesh.outdim; o++){
                              innerLocal[o] += kernel_val[o] * quadRule.dy[i] * rTdet; // Local Term
                    }
                    //innerLocal += kernel_val * quadRule.dy[i] * rTdet; // Local Term

                    // Pull resulting physical point ry to the (underlying!) reference Triangle aT.
                    toRef(bT.E, physical_quad, reference_quad);
                    // Evaluate ker on physical quad (note this is ker')
                    model_kernel(physical_quad, bT.label, x, aT.label, mesh, kernelPrime_val);
                    // Evaluate basis function on resulting reference quadrature point
                    model_basisFunction(reference_quad, mesh.dim, psi_value);

                    // [17]
                    // INNER NON-LOCAL ORDER
                    // [(b 0, ker (0,0)), (b 0, ker (0,1)), (b 0, ker (1,0)), (b 0, ker (1,1)),
                    //  (b 1, ker (0,0)), (b 1, ker (0,1)), (b 1, ker (1,0)), (b 1, ker (1,1)),
                    //  (b 2, ker (0,0)), (b 2, ker (0,1)), (b 2, ker (1,0)), (b 2, ker (1,1))]
                    //  = (PSI ORDER) * (KERNEL ORDER)
                    //TODO
                    // Function is much too complex
                    // A different expression than kernel * u * v in the integration requires a new integration routine!
                    for (b = 0; b < mesh.dVertex*mesh.outdim*mesh.outdim; b++) {
                    // for (b = 0; b < mesh.dVertex; b++) {
                        // [x 18]
                        innerNonloc[b] +=
                                psi_value[b/(mesh.outdim*mesh.outdim)] *
                                kernel_val[b%(mesh.outdim*mesh.outdim)] *
                                quadRule.dy[i] * rTdet; // Nonlocal Term
                        //innerNonloc[b] += psi_value[b] * kernel_val * quadRule.dy[i] * rTdet; // Nonlocal Term
                        innerNonlocPrime[b] +=
                                psi_value[b/(mesh.outdim*mesh.outdim)] *
                                kernelPrime_val[b%(mesh.outdim*mesh.outdim)] *
                                quadRule.dy[i] * rTdet; // Nonlocal Term

                    }
                    for (a = 0; a < mesh.dVertex; a++) {
                        for (b = 0; b < mesh.dVertex; b++) {
                            for (int o = 0; o < mesh.outdim * mesh.outdim; o++) {
                                termLocalPrime[a * mesh.outdim * mesh.outdim * mesh.dVertex +
                                               (o / mesh.outdim) * mesh.outdim * mesh.dVertex +
                                               b * mesh.outdim +
                                               o % mesh.outdim]
                                                +=  psi_value[a] *
                                                    psi_value[b] *
                                                    kernelPrime_val[o] *
                                                    quadRule.dy[i] *
                                                    quadRule.dx[k] * rTdet * aT.absDet;
                            }
                        }
                    }

                    //[DEBUG]
                    //printf("i %i \n",i);
                    //printf("GAM %17.16e\n", ker * dy[i] * rTdet);
                    //printf("Basis0 %17.16e\n", psi_value[0]);
                    //printf("Basis1 %17.16e\n", psi_value[1]);
                    //printf("Basis2 %17.16e\n", psi_value[2]);
                }
                //printf("Chris: v0 %17.16e\nv1 %17.16e\nv2 %17.16e\n", innerNonloc[0], innerNonloc[1], innerNonloc[2]);
                //printf("Chris: v %17.16e\n", innerLocal);
            }
        }

        // TERM LOCAL & TERM NON-LOCAL ORDER
        // Note: This order is not trivially obtained from innerNonloc, as b switches in between.
        // However it mimics the matrix which results from the multiplication.
        //                                      Kernel switches back here. v
        // [(a 0, b 0, ker (0,0)), (a 0, b 0, ker (0,1)), (a 0, b 1, ker (0,0)), (a 0, b 1, ker (0,1)), (a 0, b 2, ker (0,0)), (a 0, b 2, ker (0,1)),
        //  (a 0, b 0, ker (1,0)), (a 0, b 0, ker (1,1)), (a 0, b 0, ker (1,0)), (a 0, b 1, ker (1,1)), (a 0, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 1, b 0, ker (0,0)), (a 1, b 0, ker (0,1)), (a 1, b 1, ker (0,0)), (a 1, b 1, ker (0,1)), (a 1, b 2, ker (0,0)), (a 1, b 2, ker (0,1)),
        //  (a 1, b 0, ker (1,0)), (a 1, b 0, ker (1,1)), (a 1, b 0, ker (1,0)), (a 1, b 1, ker (1,1)), (a 1, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 2, b 0, ker (0,0)), (a 2, b 0, ker (0,1)), (a 2, b 1, ker (0,0)), (a 2, b 1, ker (0,1)), (a 2, b 2, ker (0,0)), (a 2, b 2, ker (0,1)),
        //  (a 2, b 0, ker (1,0)), (a 2, b 0, ker (1,1)), (a 2, b 0, ker (1,0)), (a 2, b 1, ker (1,1)), (a 2, b 2, ker (1,0)), (a 2, b 2, ker (1,1))]

        //  = (PSI ORDER) * (PSI ORDER) * (INNER LOCAL ORDER)
        //  = (PSI ORDER) *' (INNER NON-LOCAL ORDER)

        //printf("Local %17.16e\n", innerLocal);
        //printf("Nonloc [%17.16e, %17.16e, %17.16e, %17.16e] \n", innerNonloc[0], innerNonloc[1], innerNonloc[2], innerNonloc[3]);

        // [x 19] for (a = 0; a < mesh.dVertex*mesh.outdim; a++) {
        for (a = 0; a < mesh.dVertex * mesh.outdim; a++) {
            // [x 20] for (b = 0; b < mesh.dVertex*mesh.outdim; b++) {
            for (b = 0; b < mesh.dVertex * mesh.outdim; b++) {
                // [x 21] termLocal[mesh.dVertex * mesh.outputdim * a + b] +=
                termLocal[mesh.dVertex * mesh.outdim * a + b] +=
                        aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) * quadRule.dx[k] *
                        innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]; //innerLocal
                // psi_value_test[a/mesh.outdim]*psi_value_test[b/mesh.outdim]+innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)];
                //printf("a %6.4e, b %6.4e, innerLocal %6.4e \n", psi_value_test[a/mesh.outdim], psi_value_test[b/mesh.outdim], innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]);
                // [x 22] 2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) * quadRule.dx[k] * ...

                //printf("quadRule.psix(%i,%i) %17.16e\nquadRule.psix(%i,%i) %17.16e \n", a,k, quadRule.psix(a,k), b,k,quadRule.psix(b,k));
                // [x 24] termNonloc[mesh.dVertex * mesh.outputdim * a + b] +=
                termNonloc[mesh.dVertex * mesh.outdim * a + b] +=
                        aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.dx[k] *
                        innerNonloc[(a%mesh.outdim)*mesh.outdim +
                                    mesh.outdim*mesh.outdim*(b/mesh.outdim) +
                                    (b%mesh.outdim)];
                //printf("a %6.4e, innerNonloc %6.4e \n", psi_value_test[a/mesh.outdim],
                // innerNonloc[(a%mesh.outdim)*mesh.outdim + mesh.outdim*mesh.outdim*(b/mesh.outdim) + (b%mesh.outdim)]);
                // [x 25] 2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.dx[k] *
                //2 * aT.absDet * quadRule.psix(a, k) * quadRule.dx[k] * innerNonloc[b]; //innerNonloc
                termNonlocPrime[mesh.dVertex * mesh.outdim * a + b] +=
                        aT.absDet * quadRule.psix(a / mesh.outdim, k) * quadRule.dx[k] *
                        innerNonlocPrime[(a % mesh.outdim) * mesh.outdim +
                                         mesh.outdim * mesh.outdim * (b / mesh.outdim) +
                                         (b % mesh.outdim)];

            }
        }
    }
    return (doesInteract || conf.is_fullConnectedComponentSearch);
}

int integrate_retriangulate_shape(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                                  const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer,
                                  double *termLocal, double *termNonloc, double *termLocalPrime,
                                  double *termNonlocPrime, double * state_nodal_values, double * adjoint_nodal_values) {
    // See integrate_retriangulate for additional comments.
    const int dim = mesh.dim;
    int k = 0, a = 0;
    double x[dim];
    double state_values[2] = {0.0, 0.0};
    double adjoint_values[2] = {0.0, 0.0};

    double value = 0.0;

    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];
    double innerNonlocPrime[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0, rTdx = 0, Rdx = 0;
    double kernel_val[1];
    double kernelPrime_val[1];

    double rTdet = 0;
    double physical_quad[dim];
    double reference_quad[dim];
    double psi_value[mesh.dVertex];
    double reTriangle_list[36 * mesh.dVertex * dim];
    doubleVec_tozero(reTriangle_list, 36 * mesh.dVertex * dim);
    int doesInteract=0;

    for (k = 0; k < quadRule.nPx; k++) {
        toPhys(aT.E, &(quadRule.Px[dim * k]), mesh.dim, x);

        Rdx = method(x, bT, mesh, reTriangle_list, conf.is_placePointOnCap);

        doesInteract += Rdx;

        state_values[0] = 0.0;
        adjoint_values[0] = 0.0;
        for (a = 0; a < mesh.dVertex; a++) {
            state_values[0] += state_nodal_values[a] * quadRule.psix(a, k);
            adjoint_values[0] += adjoint_nodal_values[a] * quadRule.psix(a, k);
        }
        // printf("%f, %f\n", state_values[0], adjoint_values[0]);
        doubleVec_tozero(innerNonloc, mesh.outdim*mesh.outdim*mesh.dVertex);
        doubleVec_tozero(innerNonlocPrime, mesh.outdim*mesh.outdim*mesh.dVertex);

        if (Rdx == 0) {
        } else {
            for (rTdx = 0; rTdx < Rdx; rTdx++) {
                for (i = 0; i < quadRule.nPy; i++) {
                    // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                    toPhys(&reTriangle_list[dim * mesh.dVertex * rTdx], &(quadRule.Py[dim * i]), physical_quad);
                    // Determinant of Triangle of retriangulation
                    rTdet = absDet(&reTriangle_list[dim * mesh.dVertex * rTdx]);

                    model_kernel(x, aT.label, physical_quad, bT.label, mesh, kernel_val);

                    // Pull resulting physical point ry to the (underlying!) reference Triangle aT.
                    toRef(bT.E, physical_quad, reference_quad);
                    // Evaluate ker on physical quad (note this is ker')
                    model_kernel(physical_quad, bT.label, x, aT.label, mesh, kernelPrime_val);
                    // Evaluate basis function on resulting reference quadrature point
                    model_basisFunction(reference_quad, mesh.dim, psi_value);

                    state_values[1] = 0.0;
                    adjoint_values[1] = 0.0;
                    for (a = 0; a < mesh.dVertex; a++) {
                        state_values[1] += state_nodal_values[a + mesh.dVertex] * psi_value[a];
                        adjoint_values[1] += adjoint_nodal_values[a + mesh.dVertex] * psi_value[a];
                    }

                    value = (adjoint_values[0] - adjoint_values[1])
                            *(state_values[0]*kernel_val[0] - state_values[1]*kernelPrime_val[0])*quadRule.dy[i]
                            *quadRule.dx[k]*rTdet*aT.absDet;
                    termLocal[0] += value;
                    if (conf.is_singularKernel){
                        for(a = 0; a < mesh.dVertex; a++) {
                            termNonloc[a] += value * (x[0] - physical_quad[0])/(vec_sqL2dist(x, physical_quad, 2))
                                    * quadRule.psix(a, k);
                            termNonlocPrime[a] += value *
                                    (x[1] - physical_quad[1])/(vec_sqL2dist(x, physical_quad, 2)) * quadRule.psix(a, k);
                        }
                    }
                }
            }
        }
    }

    return (doesInteract || conf.is_fullConnectedComponentSearch);
}

int integrate_retriangulate_unysmm(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                             const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                             double *termNonloc, double *termLocalPrime, double *termNonlocPrime) {


    //if ((mesh.maxDiameter > EPSILON) && (mesh.delta - 2*mesh.maxDiameter > 0) && isFullyContained(aT, bT, mesh)){
    //    integrate_fullyContained(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal, termNonloc);
    //    return;
    //}

    const int dim = mesh.dim;
    int k = 0, a = 0, b = 0;
    double x[dim];
    // [x 11] [mesh.outdim*mesh.outdim]
    //double innerLocal = 0;
    double innerLocal[mesh.outdim*mesh.outdim];

    // [x 12] [mesh.outdim*mesh.outdim*mesh.dVertex]
    // double innerNonloc[mesh.dVertex];
    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0, rTdx = 0, Rdx = 0;
    // [x 13] kernel_val [mesh.outdim*mesh.outdim]
    double kernel_val[mesh.outdim*mesh.outdim];

    double rTdet = 0;
    double physical_quad[dim];
    double reference_quad[dim];
    double psi_value[mesh.dVertex];
    //double psi_value_test[3] = {20., 30., 50.};
    double reTriangle_list[36 * mesh.dVertex * dim];
    doubleVec_tozero(reTriangle_list, 36 * mesh.dVertex * dim);
    int doesInteract=0;
    //[DEBUG]
    //printf("\nouterInt_full----------------------------------------\n");
    for (k = 0; k < quadRule.nPx; k++) {
        //printf("k %i, quadRule.nPx %i\n", k, quadRule.nPx);
        toPhys(aT.E, &(quadRule.Px[dim * k]), mesh.dim, x);
        //printf("\nInner Integral, Iterate %i\n", k);
        //printf("\Physical x [%17.16e, %17.16e]\n",  x[0], x[1]);
        //innerInt_retriangulate(x, aT, bT, quadRule, sqdelta, &innerLocal, innerNonloc);

        //is_placePointOnCap = true;
        Rdx = method(x, bT, mesh, reTriangle_list, conf.is_placePointOnCap); // innerInt_retriangulate
        doesInteract += Rdx;
        //Rdx = baryCenterMethod(x, bT, mesh, reTriangle_list, is_placePointOnCap);
        //Rdx = quadRule.interactionMethod(x, bT, mesh, reTriangle_list);

        //[DEBUG]
        //printf("Retriangulation Rdx %i\n", Rdx);
        //for (i=0;i<Rdx;i++){
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i], reTriangle_list[2 * 3 * i+1]);
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i+2], reTriangle_list[2 * 3 * i+3]);
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i+4], reTriangle_list[2 * 3 * i+5]);
        //printf("absDet %17.16e\n", absDet(&reTriangle_list[2 * 3 * i]));
        //}
        // [x 14] doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);
        //innerLocal = 0.0;
        doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);

        // [x 15] doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim*mesh.dVertex);
        // doubleVec_tozero(innerNonloc, mesh.dVertex);
        doubleVec_tozero(innerNonloc, mesh.outdim*mesh.outdim*mesh.dVertex);

        if (Rdx == 0) {
        } else {
            //printf("\nInner Integral\n");
            for (rTdx = 0; rTdx < Rdx; rTdx++) {
                //printf("rTdx %i \n",rTdx);
                for (i = 0; i < quadRule.nPy; i++) {
                    // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                    toPhys(&reTriangle_list[dim * mesh.dVertex * rTdx], &(quadRule.Py[dim * i]), physical_quad);
                    // Determinant of Triangle of retriangulation
                    rTdet = absDet(&reTriangle_list[dim * mesh.dVertex * rTdx]);
                    // inner Local integral with ker
                    model_kernel(x, aT.label, physical_quad, bT.label, mesh, kernel_val);
                    // [x 16]
                    // INNER LOCAL ORDER [(0,0), (0,1), (1,0), (1,1)] = KERNEL ORDER
                    for (int o=0; o<mesh.outdim*mesh.outdim; o++){
                        innerLocal[o] += kernel_val[o] * quadRule.dy[i] * rTdet; // Local Term
                    }
                    //innerLocal += kernel_val * quadRule.dy[i] * rTdet; // Local Term

                    // Pull resulting physical point ry to the (underlying!) reference Triangle aT.
                    toRef(bT.E, physical_quad, reference_quad);
                    // Evaluate ker on physical quad (note this is ker')
                    model_kernel(physical_quad, bT.label, x, aT.label, mesh, kernel_val);
                    // Evaluate basis function on resulting reference quadrature point
                    model_basisFunction(reference_quad, mesh.dim, psi_value);

                    // [17]
                    // INNER NON-LOCAL ORDER
                    // [(b 0, ker (0,0)), (b 0, ker (0,1)), (b 0, ker (1,0)), (b 0, ker (1,1)),
                    //  (b 1, ker (0,0)), (b 1, ker (0,1)), (b 1, ker (1,0)), (b 1, ker (1,1)),
                    //  (b 2, ker (0,0)), (b 2, ker (0,1)), (b 2, ker (1,0)), (b 2, ker (1,1))]
                    //  = (PSI ORDER) * (KERNEL ORDER)

                    for (b = 0; b < mesh.dVertex*mesh.outdim*mesh.outdim; b++) {
                        // for (b = 0; b < mesh.dVertex; b++) {
                        // [x 18]
                        innerNonloc[b] +=
                                psi_value[b/(mesh.outdim*mesh.outdim)] *
                                kernel_val[b%(mesh.outdim*mesh.outdim)] *
                                quadRule.dy[i] * rTdet; // Nonlocal Term
                        //innerNonloc[b] += psi_value[b] * kernel_val * quadRule.dy[i] * rTdet; // Nonlocal Term
                    }
                    //[DEBUG]
                    //printf("i %i \n",i);
                    //printf("GAM %17.16e\n", ker * dy[i] * rTdet);
                    //printf("Basis0 %17.16e\n", psi_value[0]);
                    //printf("Basis1 %17.16e\n", psi_value[1]);
                    //printf("Basis2 %17.16e\n", psi_value[2]);
                }
                //printf("Chris: v0 %17.16e\nv1 %17.16e\nv2 %17.16e\n", innerNonloc[0], innerNonloc[1], innerNonloc[2]);
                //printf("Chris: v %17.16e\n", innerLocal);
            }
        }

        // TERM LOCAL & TERM NON-LOCAL ORDER
        // Note: This order is not trivially obtained from innerNonloc, as b switches in between.
        // However it mimics the matrix which results from the multiplication.
        //                                      Kernel switches back here. v
        // [(a 0, b 0, ker (0,0)), (a 0, b 0, ker (0,1)), (a 0, b 1, ker (0,0)), (a 0, b 1, ker (0,1)), (a 0, b 2, ker (0,0)), (a 0, b 2, ker (0,1)),
        //  (a 0, b 0, ker (1,0)), (a 0, b 0, ker (1,1)), (a 0, b 0, ker (1,0)), (a 0, b 1, ker (1,1)), (a 0, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 1, b 0, ker (0,0)), (a 1, b 0, ker (0,1)), (a 1, b 1, ker (0,0)), (a 1, b 1, ker (0,1)), (a 1, b 2, ker (0,0)), (a 1, b 2, ker (0,1)),
        //  (a 1, b 0, ker (1,0)), (a 1, b 0, ker (1,1)), (a 1, b 0, ker (1,0)), (a 1, b 1, ker (1,1)), (a 1, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 2, b 0, ker (0,0)), (a 2, b 0, ker (0,1)), (a 2, b 1, ker (0,0)), (a 2, b 1, ker (0,1)), (a 2, b 2, ker (0,0)), (a 2, b 2, ker (0,1)),
        //  (a 2, b 0, ker (1,0)), (a 2, b 0, ker (1,1)), (a 2, b 0, ker (1,0)), (a 2, b 1, ker (1,1)), (a 2, b 2, ker (1,0)), (a 2, b 2, ker (1,1))]

        //  = (PSI ORDER) * (PSI ORDER) * (INNER LOCAL ORDER)
        //  = (PSI ORDER) *' (INNER NON-LOCAL ORDER)

        //printf("Local %17.16e\n", innerLocal);
        //printf("Nonloc [%17.16e, %17.16e, %17.16e, %17.16e] \n", innerNonloc[0], innerNonloc[1], innerNonloc[2], innerNonloc[3]);

        // [x 19] for (a = 0; a < mesh.dVertex*mesh.outdim; a++) {
        for (a = 0; a < mesh.dVertex * mesh.outdim; a++) {
            // [x 20] for (b = 0; b < mesh.dVertex*mesh.outdim; b++) {
            for (b = 0; b < mesh.dVertex * mesh.outdim; b++) {
                // [x 21] termLocal[mesh.dVertex * mesh.outputdim * a + b] +=
                termLocal[mesh.dVertex * mesh.outdim * a + b] +=
                        2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) * quadRule.dx[k] *
                        innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]; //innerLocal
                // psi_value_test[a/mesh.outdim]*psi_value_test[b/mesh.outdim]+innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)];
                //printf("a %6.4e, b %6.4e, innerLocal %6.4e \n", psi_value_test[a/mesh.outdim], psi_value_test[b/mesh.outdim], innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]);
                // [x 22] 2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) * quadRule.dx[k] * ...

                //printf("quadRule.psix(%i,%i) %17.16e\nquadRule.psix(%i,%i) %17.16e \n", a,k, quadRule.psix(a,k), b,k,quadRule.psix(b,k));
                // [x 24] termNonloc[mesh.dVertex * mesh.outputdim * a + b] +=
                termNonloc[mesh.dVertex * mesh.outdim * a + b] +=
                        2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.dx[k] *
                        innerNonloc[(a%mesh.outdim)*mesh.outdim +
                                    mesh.outdim*mesh.outdim*(b/mesh.outdim) +
                                    (b%mesh.outdim)];
                //printf("a %6.4e, innerNonloc %6.4e \n", psi_value_test[a/mesh.outdim],
                // innerNonloc[(a%mesh.outdim)*mesh.outdim + mesh.outdim*mesh.outdim*(b/mesh.outdim) + (b%mesh.outdim)]);
                // [x 25] 2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.dx[k] *
                //2 * aT.absDet * quadRule.psix(a, k) * quadRule.dx[k] * innerNonloc[b]; //innerNonloc
            }
        }
    }
    return (doesInteract || conf.is_fullConnectedComponentSearch);
}

int integrate_exact(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                     const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                     double *termNonloc, double *termLocalPrime, double *termNonlocPrime) {
    // One would assume, the code became faster by the check below - but until delta = .3 it still didn't
    // Code can be more efficient, if we use very simple rules for elements which fully lie in the interaction set.
    //if ((mesh.maxDiameter > EPSILON) && (mesh.delta - 2*mesh.maxDiameter > 0) && isFullyContained(aT, bT, mesh)){
    //    return integrate_fullyContained(aT, bT, quadRule, mesh, conf, is_firstbfslayer, termLocal, termNonloc,
    //                             termLocalPrime, termNonlocPrime);
    //}

    const int dim = mesh.dim;
    int k = 0, a = 0, b = 0;
    double x[dim];
    // [x 11] [mesh.outdim*mesh.outdim]
    //double innerLocal = 0;
    double innerLocal[mesh.outdim*mesh.outdim];

    // [x 12] [mesh.outdim*mesh.outdim*mesh.dVertex]
    // double innerNonloc[mesh.dVertex];
    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0, rTdx = 0, Rdx = 0;
    // [x 13] kernel_val [mesh.outdim*mesh.outdim]
    double kernel_val[mesh.outdim*mesh.outdim];

    double rTdet = 0;
    double physical_quad[dim];
    double reference_quad[dim];
    double psi_value[mesh.dVertex];
    //double psi_value_test[3] = {20., 30., 50.};
    double reTriangle_list[36 * mesh.dVertex * dim];
    doubleVec_tozero(reTriangle_list, 36 * mesh.dVertex * dim);
    double capsList[3*2];
    doubleVec_tozero(capsList, 3*2);
    double capsWeights[3];
    doubleVec_tozero(capsWeights, 3);
    int nCaps = 0;
    int doesInteract=0;
    //[DEBUG]
    //printf("\nouterInt_full----------------------------------------\n");
    for (k = 0; k < quadRule.nPx; k++) {
        //printf("k %i, quadRule.nPx %i\n", k, quadRule.nPx);
        toPhys(aT.E, &(quadRule.Px[dim * k]), mesh.dim, x);
        //printf("\nInner Integral, Iterate %i\n", k);
        //printf("\Physical x [%17.16e, %17.16e]\n",  x[0], x[1]);
        //innerInt_retriangulate(x, aT, bT, quadRule, sqdelta, &innerLocal, innerNonloc);

        //is_placePointOnCap = true;
        Rdx = method_exact(x, bT, mesh, reTriangle_list, capsList, capsWeights, &nCaps); // innerInt_retriangulate
        doesInteract += Rdx;
        //[DEBUG]
        //printf("Retriangulation Rdx %i\n", Rdx);
        //for (i=0;i<Rdx;i++){
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i], reTriangle_list[2 * 3 * i+1]);
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i+2], reTriangle_list[2 * 3 * i+3]);
        //printf("[%17.16e, %17.16e]\n", reTriangle_list[2 * 3 * i+4], reTriangle_list[2 * 3 * i+5]);
        //printf("absDet %17.16e\n", absDet(&reTriangle_list[2 * 3 * i]));
        //}
        // [x 14] doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);
        //innerLocal = 0.0;
        doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);

        // [x 15] doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim*mesh.dVertex);
        // doubleVec_tozero(innerNonloc, mesh.dVertex);
        doubleVec_tozero(innerNonloc, mesh.outdim*mesh.outdim*mesh.dVertex);

        if (Rdx == 0) {
        } else {
            //printf("\nInner Integral\n");
            for (rTdx = 0; rTdx < Rdx; rTdx++) {
                //printf("rTdx %i \n",rTdx);
                for (i = 0; i < quadRule.nPy; i++) {
                    // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                    toPhys(&reTriangle_list[dim * mesh.dVertex * rTdx], &(quadRule.Py[dim * i]), physical_quad);
                    // Determinant of Triangle of retriangulation
                    rTdet = absDet(&reTriangle_list[dim * mesh.dVertex * rTdx]);
                    // inner Local integral with ker
                    model_kernel(x, aT.label, physical_quad, bT.label, mesh, kernel_val);
                    // [x 16]
                    // INNER LOCAL ORDER [(0,0), (0,1), (1,0), (1,1)] = KERNEL ORDER
                    for (int o=0; o<mesh.outdim*mesh.outdim; o++){
                        innerLocal[o] += kernel_val[o] * quadRule.dy[i] * rTdet; // Local Term
                    }
                    //innerLocal += kernel_val * quadRule.dy[i] * rTdet; // Local Term

                    // Pull resulting physical point ry to the (underlying!) reference Triangle aT.
                    toRef(bT.E, physical_quad, reference_quad);
                    // Evaluate ker on physical quad (note this is ker')
                    model_kernel(physical_quad, bT.label, x, aT.label, mesh, kernel_val);
                    // Evaluate basis function on resulting reference quadrature point
                    model_basisFunction(reference_quad, mesh.dim, psi_value);

                    // [17]
                    // INNER NON-LOCAL ORDER
                    // [(b 0, ker (0,0)), (b 0, ker (0,1)), (b 0, ker (1,0)), (b 0, ker (1,1)),
                    //  (b 1, ker (0,0)), (b 1, ker (0,1)), (b 1, ker (1,0)), (b 1, ker (1,1)),
                    //  (b 2, ker (0,0)), (b 2, ker (0,1)), (b 2, ker (1,0)), (b 2, ker (1,1))]
                    //  = (PSI ORDER) * (KERNEL ORDER)

                    for (b = 0; b < mesh.dVertex*mesh.outdim*mesh.outdim; b++) {
                        // for (b = 0; b < mesh.dVertex; b++) {
                        // [x 18]
                        innerNonloc[b] +=
                                psi_value[b/(mesh.outdim*mesh.outdim)] *
                                kernel_val[b%(mesh.outdim*mesh.outdim)] *
                                quadRule.dy[i] * rTdet; // Nonlocal Term
                        //innerNonloc[b] += psi_value[b] * kernel_val * quadRule.dy[i] * rTdet; // Nonlocal Term
                    }
                    //[DEBUG]
                    //printf("i %i \n",i);
                    //printf("GAM %17.16e\n", ker * dy[i] * rTdet);
                    //printf("Basis0 %17.16e\n", psi_value[0]);
                    //printf("Basis1 %17.16e\n", psi_value[1]);
                    //printf("Basis2 %17.16e\n", psi_value[2]);
                }
                //printf("Chris: v0 %17.16e\nv1 %17.16e\nv2 %17.16e\n", innerNonloc[0], innerNonloc[1], innerNonloc[2]);
                //printf("Chris: v %17.16e\n", innerLocal);
            }
        }

        // integration for the caps
        for (i = 0; i < nCaps; i++) {
            //cout << "hallo" << endl;
            //abort();
            // inner Local integral with ker
            model_kernel(x, aT.label, &capsList[2*i], bT.label, mesh, kernel_val);
            // INNER LOCAL ORDER [(0,0), (0,1), (1,0), (1,1)] = KERNEL ORDER
            for (int o=0; o<mesh.outdim*mesh.outdim; o++){
                innerLocal[o] += kernel_val[o] * capsWeights[i]; // Local Term
            }
            //innerLocal += kernel_val * quadRule.dy[i] * rTdet; // Local Term

            // Pull resulting physical point ry to the (underlying!) reference Triangle aT.
            toRef(bT.E, &capsList[2*i], reference_quad);
            // Evaluate ker on physical quad (note this is ker')
            model_kernel(&capsList[2*i], bT.label, x, aT.label, mesh, kernel_val);
            // Evaluate basis function on resulting reference quadrature point
            model_basisFunction(reference_quad, mesh.dim, psi_value);

            for (b = 0; b < mesh.dVertex*mesh.outdim*mesh.outdim; b++) {
                // for (b = 0; b < mesh.dVertex; b++) {
                // [x 18]
                innerNonloc[b] +=
                        psi_value[b/(mesh.outdim*mesh.outdim)] *
                        kernel_val[b%(mesh.outdim*mesh.outdim)] *
                        capsWeights[i]; // Nonlocal Term
                //innerNonloc[b] += psi_value[b] * kernel_val * quadRule.dy[i] * rTdet; // Nonlocal Term
            }

        }







        // TERM LOCAL & TERM NON-LOCAL ORDER
        // Note: This order is not trivially obtained from innerNonloc, as b switches in between.
        // However it mimics the matrix which results from the multiplication.
        //                                      Kernel switches back here. v
        // [(a 0, b 0, ker (0,0)), (a 0, b 0, ker (0,1)), (a 0, b 1, ker (0,0)), (a 0, b 1, ker (0,1)), (a 0, b 2, ker (0,0)), (a 0, b 2, ker (0,1)),
        //  (a 0, b 0, ker (1,0)), (a 0, b 0, ker (1,1)), (a 0, b 0, ker (1,0)), (a 0, b 1, ker (1,1)), (a 0, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 1, b 0, ker (0,0)), (a 1, b 0, ker (0,1)), (a 1, b 1, ker (0,0)), (a 1, b 1, ker (0,1)), (a 1, b 2, ker (0,0)), (a 1, b 2, ker (0,1)),
        //  (a 1, b 0, ker (1,0)), (a 1, b 0, ker (1,1)), (a 1, b 0, ker (1,0)), (a 1, b 1, ker (1,1)), (a 1, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 2, b 0, ker (0,0)), (a 2, b 0, ker (0,1)), (a 2, b 1, ker (0,0)), (a 2, b 1, ker (0,1)), (a 2, b 2, ker (0,0)), (a 2, b 2, ker (0,1)),
        //  (a 2, b 0, ker (1,0)), (a 2, b 0, ker (1,1)), (a 2, b 0, ker (1,0)), (a 2, b 1, ker (1,1)), (a 2, b 2, ker (1,0)), (a 2, b 2, ker (1,1))]

        //  = (PSI ORDER) * (PSI ORDER) * (INNER LOCAL ORDER)
        //  = (PSI ORDER) *' (INNER NON-LOCAL ORDER)

        //printf("Local %17.16e\n", innerLocal);
        //printf("Nonloc [%17.16e, %17.16e, %17.16e, %17.16e] \n", innerNonloc[0], innerNonloc[1], innerNonloc[2], innerNonloc[3]);

        // [x 19] for (a = 0; a < mesh.dVertex*mesh.outdim; a++) {
        for (a = 0; a < mesh.dVertex * mesh.outdim; a++) {
            // [x 20] for (b = 0; b < mesh.dVertex*mesh.outdim; b++) {
            for (b = 0; b < mesh.dVertex * mesh.outdim; b++) {
                // [x 21] termLocal[mesh.dVertex * mesh.outputdim * a + b] +=
                termLocal[mesh.dVertex * mesh.outdim * a + b] +=
                        2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) * quadRule.dx[k] *
                        innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]; //innerLocal
                // psi_value_test[a/mesh.outdim]*psi_value_test[b/mesh.outdim]+innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)];
                //printf("a %6.4e, b %6.4e, innerLocal %6.4e \n", psi_value_test[a/mesh.outdim], psi_value_test[b/mesh.outdim], innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]);
                // [x 22] 2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) * quadRule.dx[k] * ...

                //printf("quadRule.psix(%i,%i) %17.16e\nquadRule.psix(%i,%i) %17.16e \n", a,k, quadRule.psix(a,k), b,k,quadRule.psix(b,k));
                // [x 24] termNonloc[mesh.dVertex * mesh.outputdim * a + b] +=
                termNonloc[mesh.dVertex * mesh.outdim * a + b] +=
                        2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.dx[k] *
                        innerNonloc[(a%mesh.outdim)*mesh.outdim +
                                    mesh.outdim*mesh.outdim*(b/mesh.outdim) +
                                    (b%mesh.outdim)];
                //printf("a %6.4e, innerNonloc %6.4e \n", psi_value_test[a/mesh.outdim],
                // innerNonloc[(a%mesh.outdim)*mesh.outdim + mesh.outdim*mesh.outdim*(b/mesh.outdim) + (b%mesh.outdim)]);
                // [x 25] 2 * aT.absDet * quadRule.psix(a/mesh.outdim, k) * quadRule.dx[k] *
                //2 * aT.absDet * quadRule.psix(a, k) * quadRule.dx[k] * innerNonloc[b]; //innerNonloc
            }
        }
    }
    return (doesInteract || conf.is_fullConnectedComponentSearch);
}


int
integrate_baryCenter(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                     const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                     double *termLocalPrime, double *termNonlocPrime) {
    const int dim = mesh.dim;
    int k = 0, a = 0, b = 0;
    double x[dim];
    // [x 26]
    double innerLocal[mesh.outdim*mesh.outdim];

    // [x 27]
    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0;
    double kernel_val[mesh.outdim*mesh.outdim];
    double rTdet = 0.0;
    double physical_quad[dim];
    double reTriangle_list[36 * mesh.dVertex * dim];
    doubleVec_tozero(reTriangle_list, 36 * mesh.dVertex * dim);
    int doesInteract=0, Rdx=0;

    //[DEBUG]
    //printf("\nouterInt_full----------------------------------------\n");
    for (k = 0; k < quadRule.nPx; k++) {
        toPhys(aT.E, &(quadRule.Px[dim * k]), mesh.dim, x);

        doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);
        // [x 28]
        doubleVec_tozero(innerNonloc, mesh.outdim*mesh.outdim*mesh.dVertex);
        Rdx = method_baryCenter(x, bT, mesh, reTriangle_list, false);
        doesInteract+=Rdx;
        if (Rdx) {
            for (i = 0; i < quadRule.nPy; i++) {
                // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                toPhys(bT.E, &(quadRule.Py[dim * i]), mesh.dim, physical_quad);
                // Determinant of Triangle of retriangulation
                rTdet = absDet(bT.E, mesh.dim);
                // inner Local integral with ker
                model_kernel(x, aT.label, physical_quad, bT.label, mesh, kernel_val);
                // [x 29]
                for (int o=0; o<mesh.outdim*mesh.outdim; o++) {
                    innerLocal[o] += rTdet * kernel_val[o] * quadRule.dy[i]; // Local Term
                }
                // Evaluate ker on physical quad (note this is ker')
                model_kernel(physical_quad, bT.label, x, aT.label, mesh, kernel_val);
                // Evaluate basis function on resulting reference quadrature point

                // INNER NON-LOCAL ORDER
                // [(b 0, ker (0,0)), (b 0, ker (0,1)), (b 0, ker (1,0)), (b 0, ker (1,1)),
                //  (b 1, ker (0,0)), (b 1, ker (0,1)), (b 1, ker (1,0)), (b 1, ker (1,1)),
                //  (b 2, ker (0,0)), (b 2, ker (0,1)), (b 2, ker (1,0)), (b 2, ker (1,1))]
                //  = (PSI ORDER) * (KERNEL ORDER)

                for (b = 0; b < mesh.dVertex*mesh.outdim*mesh.outdim; b++) {
                    // [x 31]
                    //innerNonloc[b] += quadRule.psiy(b, i) * kernel_val * quadRule.dy[i] * rTdet; // Nonlocal Term
                    innerNonloc[b] +=
                            quadRule.psiy(b/(mesh.outdim*mesh.outdim), i) *
                            kernel_val[b%(mesh.outdim*mesh.outdim)] *
                            quadRule.dy[i] * rTdet; // Nonlocal Term
                }
            }
        }
        //printf("Local %17.16e\n", innerLocal);
        //printf("Nonloc [%17.16e, %17.16e, %17.16e, %17.16e] \n", innerNonloc[0], innerNonloc[1], innerNonloc[2], innerNonloc[3]);

        // TERM LOCAL & TERM NON-LOCAL ORDER
        // Note: This order is not trivially obtained from innerNonloc, as b switches in between.
        // However it mimics the matrix which results from the multiplication.
        //                                      Kernel switches back here. v
        // [(a 0, b 0, ker (0,0)), (a 0, b 0, ker (0,1)), (a 0, b 1, ker (0,0)), (a 0, b 1, ker (0,1)), (a 0, b 2, ker (0,0)), (a 0, b 2, ker (0,1)),
        //  (a 0, b 0, ker (1,0)), (a 0, b 0, ker (1,1)), (a 0, b 0, ker (1,0)), (a 0, b 1, ker (1,1)), (a 0, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 1, b 0, ker (0,0)), (a 1, b 0, ker (0,1)), (a 1, b 1, ker (0,0)), (a 1, b 1, ker (0,1)), (a 1, b 2, ker (0,0)), (a 1, b 2, ker (0,1)),
        //  (a 1, b 0, ker (1,0)), (a 1, b 0, ker (1,1)), (a 1, b 0, ker (1,0)), (a 1, b 1, ker (1,1)), (a 1, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 2, b 0, ker (0,0)), (a 2, b 0, ker (0,1)), (a 2, b 1, ker (0,0)), (a 2, b 1, ker (0,1)), (a 2, b 2, ker (0,0)), (a 2, b 2, ker (0,1)),
        //  (a 2, b 0, ker (1,0)), (a 2, b 0, ker (1,1)), (a 2, b 0, ker (1,0)), (a 2, b 1, ker (1,1)), (a 2, b 2, ker (1,0)), (a 2, b 2, ker (1,1))]

        //  = (PSI ORDER) * (PSI ORDER) * (INNER LOCAL ORDER)
        //  = (PSI ORDER) *' (INNER NON-LOCAL ORDER)

        // [x 32]
        for (a = 0; a < mesh.dVertex*mesh.outdim; a++) {
            // [x 33]
            for (b = 0; b < mesh.dVertex*mesh.outdim; b++) {
                // [x 33]
                termLocal[mesh.dVertex * mesh.outdim * a + b] +=
                    2 * aT.absDet *
                    quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) *
                    quadRule.dx[k] *
                    innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]; //innerLocal
                    // [x 34]
                    // [x 35]
                //2 * aT.absDet * quadRule.psix(a, k) * quadRule.psix(b, k) * quadRule.dx[k] * innerLocal;

                // [x 36]
                termNonloc[mesh.dVertex * mesh.outdim * a + b] +=
                // [x 37]
                //2 * aT.absDet * quadRule.psix(a, k) * quadRule.dx[k] * innerNonloc[b];
                    2 * aT.absDet *
                    quadRule.psix(a/mesh.outdim, k) *
                    quadRule.dx[k] *
                    innerNonloc[(a%mesh.outdim)*mesh.outdim +
                                mesh.outdim*mesh.outdim*(b/mesh.outdim) +
                                (b%mesh.outdim)];//innerNonloc
            }
        }
    }
    return (doesInteract  ||  conf.is_fullConnectedComponentSearch);
}


int
integrate_subSuperSetBalls(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                     const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                     double *termLocalPrime, double *termNonlocPrime) {

    double averageWeights[mesh.dVertex];
    doubleVec_tozero(averageWeights, mesh.dVertex);
    //int doesInteract;

    if (conf.integration_method_remote == "superSetBall"){
        //doesInteract = 1;
        for (auto ptr = averageWeights; ptr < averageWeights + mesh.dVertex; ptr++){
            *ptr = 1.0;
        }
    } else if (conf.integration_method_remote == "subSetBall"){
        //doesInteract = 3;
        averageWeights[2]=1.0;
    } else if (conf.integration_method_remote == "averageBall") {
        for (auto ptr = averageWeights; ptr < averageWeights + mesh.dVertex; ptr++){
            *ptr = 0.5;
        }
        averageWeights[2]=1.0;
    } else {
            cout << "Error in integrate_subSuperSetBalls: No such integration_method: " <<
                 conf.integration_method_remote << endl;
            abort();
    }

    const int dim = mesh.dim;
    int k = 0, a = 0, b = 0;
    double x[dim];
    // [x 26]
    double innerLocal[mesh.outdim*mesh.outdim];

    // [x 27]
    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0;
    double kernel_val[mesh.outdim*mesh.outdim];
    double rTdet = 0.0;
    double physical_quad[dim];
    int doesInteract=0;

    //[DEBUG]
    //printf("\nouterInt_full----------------------------------------\n");
    for (k = 0; k < quadRule.nPx; k++) {
        toPhys(aT.E, &(quadRule.Px[dim * k]), mesh.dim, x);

        doubleVec_tozero(innerLocal, mesh.outdim*mesh.outdim);
        // [x 28]
        doubleVec_tozero(innerNonloc, mesh.outdim*mesh.outdim*mesh.dVertex);
        int nContained = method_subSuperSetBalls(x, bT, mesh);
        doesInteract+=nContained;
        if (nContained >= 1) {
            //averageWeight = isAverage ? 0.5 * static_cast<double>(interaction) : 1.0;

            for (i = 0; i < quadRule.nPy; i++) {
                // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                toPhys(bT.E, &(quadRule.Py[dim * i]), mesh.dim, physical_quad);
                // Determinant of Triangle of retriangulation
                rTdet = absDet(bT.E, mesh.dim);
                // inner Local integral with ker
                model_kernel(x, aT.label, physical_quad, bT.label, mesh, kernel_val);
                // [x 29]
                for (int o=0; o<mesh.outdim*mesh.outdim; o++) {
                    innerLocal[o] += averageWeights[nContained-1] * rTdet * kernel_val[o] * quadRule.dy[i]; // Local Term
                }
                // Evaluate ker on physical quad (note this is ker')
                model_kernel(physical_quad, bT.label, x, aT.label, mesh, kernel_val);
                // Evaluate basis function on resulting reference quadrature point

                // INNER NON-LOCAL ORDER
                // [(b 0, ker (0,0)), (b 0, ker (0,1)), (b 0, ker (1,0)), (b 0, ker (1,1)),
                //  (b 1, ker (0,0)), (b 1, ker (0,1)), (b 1, ker (1,0)), (b 1, ker (1,1)),
                //  (b 2, ker (0,0)), (b 2, ker (0,1)), (b 2, ker (1,0)), (b 2, ker (1,1))]
                //  = (PSI ORDER) * (KERNEL ORDER)

                for (b = 0; b < mesh.dVertex*mesh.outdim*mesh.outdim; b++) {
                    // [x 31]
                    //innerNonloc[b] += quadRule.psiy(b, i) * kernel_val * quadRule.dy[i] * rTdet; // Nonlocal Term
                    innerNonloc[b] += averageWeights[nContained-1] *
                            quadRule.psiy(b/(mesh.outdim*mesh.outdim), i) *
                            kernel_val[b%(mesh.outdim*mesh.outdim)] *
                            quadRule.dy[i] * rTdet; // Nonlocal Term
                }
            }
        }
        //printf("Local %17.16e\n", innerLocal);
        //printf("Nonloc [%17.16e, %17.16e, %17.16e, %17.16e] \n", innerNonloc[0], innerNonloc[1], innerNonloc[2], innerNonloc[3]);

        // TERM LOCAL & TERM NON-LOCAL ORDER
        // Note: This order is not trivially obtained from innerNonloc, as b switches in between.
        // However it mimics the matrix which results from the multiplication.
        //                                      Kernel switches back here. v
        // [(a 0, b 0, ker (0,0)), (a 0, b 0, ker (0,1)), (a 0, b 1, ker (0,0)), (a 0, b 1, ker (0,1)), (a 0, b 2, ker (0,0)), (a 0, b 2, ker (0,1)),
        //  (a 0, b 0, ker (1,0)), (a 0, b 0, ker (1,1)), (a 0, b 0, ker (1,0)), (a 0, b 1, ker (1,1)), (a 0, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 1, b 0, ker (0,0)), (a 1, b 0, ker (0,1)), (a 1, b 1, ker (0,0)), (a 1, b 1, ker (0,1)), (a 1, b 2, ker (0,0)), (a 1, b 2, ker (0,1)),
        //  (a 1, b 0, ker (1,0)), (a 1, b 0, ker (1,1)), (a 1, b 0, ker (1,0)), (a 1, b 1, ker (1,1)), (a 1, b 2, ker (1,0)), (a 0, b 2, ker (1,1)),
        //  (a 2, b 0, ker (0,0)), (a 2, b 0, ker (0,1)), (a 2, b 1, ker (0,0)), (a 2, b 1, ker (0,1)), (a 2, b 2, ker (0,0)), (a 2, b 2, ker (0,1)),
        //  (a 2, b 0, ker (1,0)), (a 2, b 0, ker (1,1)), (a 2, b 0, ker (1,0)), (a 2, b 1, ker (1,1)), (a 2, b 2, ker (1,0)), (a 2, b 2, ker (1,1))]

        //  = (PSI ORDER) * (PSI ORDER) * (INNER LOCAL ORDER)
        //  = (PSI ORDER) *' (INNER NON-LOCAL ORDER)

        // [x 32]
        for (a = 0; a < mesh.dVertex*mesh.outdim; a++) {
            // [x 33]
            for (b = 0; b < mesh.dVertex*mesh.outdim; b++) {
                // [x 33]
                termLocal[mesh.dVertex * mesh.outdim * a + b] +=
                        2 * aT.absDet *
                        quadRule.psix(a/mesh.outdim, k) * quadRule.psix(b/mesh.outdim, k) *
                        quadRule.dx[k] *
                        innerLocal[mesh.outdim*(a%mesh.outdim) + (b%mesh.outdim)]; //innerLocal
                // [x 34]
                // [x 35]
                //2 * aT.absDet * quadRule.psix(a, k) * quadRule.psix(b, k) * quadRule.dx[k] * innerLocal;

                // [x 36]
                termNonloc[mesh.dVertex * mesh.outdim * a + b] +=
                        // [x 37]
                        //2 * aT.absDet * quadRule.psix(a, k) * quadRule.dx[k] * innerNonloc[b];
                        2 * aT.absDet *
                        quadRule.psix(a/mesh.outdim, k) *
                        quadRule.dx[k] *
                        innerNonloc[(a%mesh.outdim)*mesh.outdim +
                                    mesh.outdim*mesh.outdim*(b/mesh.outdim) +
                                    (b%mesh.outdim)];//innerNonloc
            }
        }
    }
    return (doesInteract  || conf.is_fullConnectedComponentSearch);
}

int integrate_baryCenterRT(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule,
                            const MeshType &mesh, const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal,
                            double *termNonloc, double *termLocalPrime, double *termNonlocPrime) {
    const int dim = mesh.dim;
    int k = 0, a = 0, b = 0;
    double physical_quad[dim], reference_quad[dim], psix[mesh.dVertex];
    double bTbaryC[dim];
    // [x 37]
    double innerLocal[mesh.outdim*mesh.outdim];// = 0;
    // [x 38]
    double innerNonloc[mesh.outdim*mesh.outdim*mesh.dVertex];

    int i = 0, Rdx, rTdx;
    // [x 39]
    double kernel_val[mesh.outdim*mesh.outdim];
    double rTdet = 0, bTdet = 0;
    double y[dim];
    double reTriangle_list[36 * mesh.dVertex * dim];
    doubleVec_tozero(reTriangle_list, 36 * mesh.dVertex * dim);

    //[DEBUG]
    //printf("\nouterInt_full----------------------------------------\n");
    baryCenter(bT.E, &bTbaryC[0]);
    Rdx = method_retriangulate(bTbaryC, aT, mesh, reTriangle_list, conf.is_placePointOnCap);

    if (!Rdx) {
        return conf.is_fullConnectedComponentSearch;
    } else {
        // Determinant of Triangle of retriangulation
        bTdet = absDet(bT.E);

        for (rTdx = 0; rTdx < Rdx; rTdx++) {
            rTdet = absDet(&reTriangle_list[dim * mesh.dVertex * rTdx]);

            for (k = 0; k < quadRule.nPx; k++) {
                toPhys(&reTriangle_list[dim * mesh.dVertex * rTdx], &(quadRule.Px[dim * k]), mesh.dim,
                       physical_quad);

                // Compute Integral over Triangle bT
                doubleVec_tozero(innerLocal, mesh.outdim * mesh.outdim);
                doubleVec_tozero(innerNonloc, mesh.dVertex);

                for (i = 0; i < quadRule.nPy; i++) {
                    // Push quadrature point P[i] to physical triangle reTriangle_list[rTdx] (of the retriangulation!)
                    toPhys(bT.E, &(quadRule.Py[dim * i]), mesh.dim, y);
                    // inner Local integral with ker
                    // Local Term
                    model_kernel(physical_quad, aT.label, y, bT.label, mesh, kernel_val);
                    // [x 40]
                    for (int o = 0; o < mesh.outdim * mesh.outdim; o++) {
                        // innerLocal[o] += kernel_val[o] * quadRule.dy[i] * rTdet; // Local Term
                        innerLocal[o] += kernel_val[o] * quadRule.dy[i] * bTdet;
                    }
                    // Evaluate kernel on physical quad (note this is kernel')
                    model_kernel(y, bT.label, physical_quad, aT.label, mesh, kernel_val);
                    // Evaluate basis function on resulting reference quadrature point
                    // [x 41]
                    for (b = 0; b < mesh.dVertex * mesh.outdim * mesh.outdim; b++) {
                        // [42]
                        innerNonloc[b] += quadRule.psiy(b / (mesh.outdim * mesh.outdim), i) *
                                          // kernel_val[b%(mesh.outputdim*mesh.outputdim)] * quadRule.dy[i] * rTdet; // Nonlocal Term
                                          kernel_val[b % (mesh.outdim * mesh.outdim)] *
                                          quadRule.dy[i] * bTdet; // Nonlocal Term

                        //innerNonloc[b] +=
                        //        quadRule.psiy(b, i) * kernel_val * quadRule.dy[i] * bTdet; // Nonlocal Term
                    }
                }

                toRef(aT.E, physical_quad, reference_quad);
                model_basisFunction(reference_quad, psix);

                // [x 43]
                for (a = 0; a < mesh.dVertex * mesh.outdim; a++) {
                    // [x 44]
                    for (b = 0; b < mesh.dVertex * mesh.outdim; b++) {
                        // [x 45]
                        termLocal[mesh.dVertex * mesh.outdim * a + b] +=
                                //2 * rTdet * psix[a] * psix[b] * quadRule.dx[k] * innerLocal; //innerLocal
                                2 * rTdet * psix[a / mesh.outdim] * psix[b / mesh.outdim] *
                                quadRule.dx[k] * innerLocal[mesh.outdim * (a % mesh.outdim) + (b % mesh.outdim)];

                        // [x 47]
                        termNonloc[mesh.dVertex * mesh.outdim * a + b] +=
                                // [48]
                                2 * rTdet * psix[a / mesh.outdim] * quadRule.dx[k] *
                                innerNonloc[mesh.outdim * (a % mesh.outdim) + b % mesh.outdim]; //innerNonloc
                        //2 * rTdet * psix[a] * quadRule.dx[k] * innerNonloc[b]; //innerNonloc
                    }
                }
            }
        }
        return (Rdx  || conf.is_fullConnectedComponentSearch);
    }
}
// Helpers -------------------------------------------------------------------------------------------------------------
// Sub Set Method ------------------------------------------------------------------------------------------------------
// Super, Sub and Average Set Method -----------------------------------------------------------------------------------
int method_subSuperSetBalls(const double * x, const ElementType & T, const MeshType & mesh){
    int nContained = 0;
    for (int i=0; i<mesh.dVertex; i++){
        nContained += (vec_sqL2dist(x, &T.E[mesh.dim * i], mesh.dim) < mesh.sqdelta);
    }
    return nContained;
}

// Bary Center Method --------------------------------------------------------------------------------------------------
int method_baryCenter(const double * x_center, const ElementType & T, const MeshType & mesh, double * reTriangle_list, int is_placePointOnCap){
    int i,k;
    double distance;
    arma::vec baryC(mesh.dim);
    //void baryCenter(const int dim, const double * E, double * bary);
    baryCenter(mesh.dim, T.E, &baryC[0]);
    distance = vec_sqL2dist(x_center, &baryC[0], mesh.dim);

    if (distance > mesh.sqdelta){
        return 0;
    } else {
        for (i=0; i<mesh.dim; i++){
            for(k=0; k<mesh.dVertex; k++) {
                reTriangle_list[2 * k + i] = T.E[2 * k + i];
            }
        }
        return -1;
    }
}

// Retriangulation Method ----------------------------------------------------------------------------------------------
bool inTriangle(const double * y_new, const double * p, const double * q, const double * r,
                const double *  nu_a, const double * nu_b, const double * nu_c){
    bool a, b, c;
    double vec[2];

    doubleVec_subtract(y_new, p, vec, 2);
    a = vec_dot(nu_a, vec , 2) >= 0;

    doubleVec_subtract(y_new, q, vec, 2);
    b = vec_dot(nu_b, vec , 2) >= 0;

    doubleVec_subtract(y_new, r, vec, 2);
    c = vec_dot(nu_c, vec , 2) >= 0;

    return a && b && c;
}

double placePointCapCenter(const double * y_predecessor, const double * y_current,
                           const double * x_center, const double sqdelta, const double * TE,
                           const double * nu_a, const double * nu_b, const double * nu_c,
                           const double orientation, double * capsList){
    // Place a point on the cap.
    //y_predecessor = &R[2*(Rdx-1)];
    double capCentroid[2], s_midpoint[2], s_projectionDirection[2];
    double scalingFactor;
    double alpha = 0;
    double p1[2];
    double p2[2];

    doubleVec_scale(-1, y_predecessor, p1, 2);
    doubleVec_scale(-1, y_current, p2, 2);
    doubleVec_add(x_center, p1, p1, 2);
    doubleVec_add(x_center, p2, p2, 2);
    alpha = 0.5 * acos(vec_dot(p1, p2,2) / sqdelta );


    doubleVec_midpoint(y_predecessor, y_current, s_midpoint, 2);
    // Note, this yields the left normal from y_predecessor to y0
    rightNormal(y_current, y_predecessor, orientation, s_projectionDirection);

    double E[6];
    double scalingJohn = sqrt( sqdelta / vec_dot(s_projectionDirection, s_projectionDirection, 2));
    doubleVec_copyTo(y_predecessor, &E[0], 2);
    doubleVec_copyTo(y_current, &E[2], 2);
    doubleVec_copyTo(s_projectionDirection, &E[4], 2);
    doubleVec_scale(scalingJohn, &E[4], &E[4], 2);
    doubleVec_add(x_center, &E[4], &E[4], 2);

    scalingFactor = (4*sqrt(sqdelta)*pow(sin(alpha),3)) / (3*(2*alpha-sin(2*alpha))) * sqrt( 1. / vec_dot(s_projectionDirection, s_projectionDirection, 2)); // here change!!!
    doubleVec_scale(scalingFactor, s_projectionDirection, s_projectionDirection, 2);
    doubleVec_add(x_center, s_projectionDirection, capCentroid, 2);

    // a = y_predecessor, b = y_current, c = scalingJohn * s_projectionDirection
    //printf("P1^TP2 %f\n", vec_dot(p1, p2,2) / sqdelta);
    //printf("alpha %f\n", alpha);
    //printf("scaling cap, scaling john %f, %f\n", (4*sqrt(sqdelta)*pow(sin(alpha),3)) / (3*(2*alpha-sin(2*alpha))),sqrt( sqdelta));


    if ( inTriangle(capCentroid, &TE[0], &TE[2], &TE[4], nu_a, nu_b, nu_c)){
        // Append capCentroid (Point on the cap)
        doubleVec_copyTo(capCentroid, capsList, 2);
        //printf("Approx Caps %f, Exact Caps: %f\n", area, (sqdelta/2) * (2*alpha - sin(2*alpha)));
        return (sqdelta/2) * (2*alpha - sin(2*alpha));
    } else {
        return 0.0;
    }
}



int placePointOnCap(const double * y_predecessor, const double * y_current,
                    const double * x_center, const double sqdelta, const double * TE,
                    const double * nu_a, const double * nu_b, const double * nu_c,
                    const double orientation, const int Rdx, double * R){
    // Place a point on the cap.
    //y_predecessor = &R[2*(Rdx-1)];
    double y_new[2], s_midpoint[2], s_projectionDirection[2];
    double scalingFactor;

    doubleVec_midpoint(y_predecessor, y_current, s_midpoint, 2);
    // Note, this yields the left normal from y_predecessor to y0
    rightNormal(y_current, y_predecessor, orientation, s_projectionDirection);
    // Simple way
    scalingFactor = sqrt( sqdelta / vec_dot(s_projectionDirection, s_projectionDirection, 2));
    doubleVec_scale(scalingFactor, s_projectionDirection, s_projectionDirection, 2);
    doubleVec_add(x_center, s_projectionDirection, y_new, 2);

    if ( inTriangle(y_new, &TE[0], &TE[2], &TE[4], nu_a, nu_b, nu_c)){
        // Append y_new (Point on the cap)
        doubleVec_copyTo(y_new, &R[2*Rdx], 2);
        return 1;
    } else {
        return 0;
    }
}

bool isFullyContained(const ElementType & aT, const ElementType & bT, const MeshType & mesh){
    //cout << "Max Diameter: " << mesh.maxDiameter << endl;
    //abort();

    double abary[mesh.dim], bbary[mesh.dim];

    baryCenter(mesh.dim, aT.E, abary);
    baryCenter(mesh.dim, bT.E, bbary);

    double l2dist = vec_sqL2dist(abary, bbary, mesh.dim);

    return ((mesh.delta - 2*mesh.maxDiameter) > 0) && ( l2dist < pow(mesh.delta - 2*mesh.maxDiameter,2) );
}

int method_retriangulate(const double * xCenter, const ElementType & T,
                         const MeshType & mesh, double * reTriangleList,
                         int isPlacePointOnCap){
    // C Variables and Arrays.
    int i=0, k=0, edgdx0=0, edgdx1=0, Rdx=0;
    double v=0, lam1=0, lam2=0, term1=0, term2=0;
    double nu_a[2], nu_b[2], nu_c[2]; // Normals
    arma::vec p(2);
    arma::vec q(2);
    arma::vec a(2);
    arma::vec b(2);
    arma::vec y1(2);
    arma::vec y2(2);
    arma::vec vec_x_center(xCenter, 2);
    double orientation;

    bool is_onEdge=false, is_firstPointLiesOnVertex=true;
    // The upper bound for the number of required points is 9
    // Hence 9*2 is an upper bound to encode all resulting triangles
    // Hence we can hardcode how much space needs to bee allocated
    // (This upper bound is thight! Check Christian Vollmann's thesis for more information.)

    double R[9*2]; // Vector containing all intersection points.
    doubleVec_tozero(R, 9*2);

    // Compute Normals of the Triangle
    orientation = -signDet(T.E);
    rightNormal(&T.E[0], &T.E[2], orientation, nu_a);
    rightNormal(&T.E[2], &T.E[4], orientation, nu_b);
    rightNormal(&T.E[4], &T.E[0], orientation, nu_c);

    for (k=0; k<3; k++){
        edgdx0 = k;
        edgdx1 = (k+1) % 3;

        doubleVec_copyTo(&T.E[2*edgdx0], &p[0], 2);
        doubleVec_copyTo(&T.E[2*edgdx1], &q[0], 2);

        a = q - vec_x_center;
        b = p - q;

        if (vec_sqL2dist(&p[0], xCenter, 2) <= mesh.sqdelta){
            doubleVec_copyTo(&p[0], &R[2*Rdx], 2);
            is_onEdge = false; // This point does not lie on the edge.
            Rdx += 1;
        }
        // PQ-Formula to solve quadratic problem
        v = pow( dot(a, b), 2) - (dot(a, a) - mesh.sqdelta) * dot(b, b);
        // If there is no sol to the quadratic problem, there is nothing to do.
        if (v >= 0){
            term1 = - dot(a, b) / dot(b, b);
            term2 = sqrt(v) / dot(b, b);

            // Vieta's Formula for computing the roots
            if (term1 > 0){
                lam1 = term1 + term2;
                lam2 = 1/lam1*(dot(a, a) - mesh.sqdelta) / dot(b, b);
            } else {
                lam2 = term1 - term2;
                lam1 = 1/lam2*(dot(a, a) - mesh.sqdelta) / dot(b, b);
            }
            y1 = lam1*b + q;
            y2 = lam2*b + q;

            // Check whether the first lambda "lies on the Triangle".
            if ((0 <= lam1) && (lam1 <= 1)){
                is_firstPointLiesOnVertex = is_firstPointLiesOnVertex && (bool)Rdx;
                // Check whether the predecessor lied on the edge
                if (is_onEdge && isPlacePointOnCap){
                    Rdx += placePointOnCap(&R[2*(Rdx-1)], &y1[0], xCenter, mesh.sqdelta, T.E, nu_a, nu_b, nu_c, orientation, Rdx, R);
                }
                // Append y1
                doubleVec_copyTo(&y1[0], &R[2*Rdx], 2);
                is_onEdge = true; // This point lies on the edge.
                Rdx += 1;
            }
            // Check whether the second lambda "lies on the Triangle".
            if ((0 <= lam2) && (lam2 <= 1) && (scal_sqL2dist(lam1, lam2) > 0)){
                is_firstPointLiesOnVertex = is_firstPointLiesOnVertex && (bool)Rdx;

                // Check whether the predecessor lied on the edge
                if (is_onEdge && isPlacePointOnCap){
                    Rdx += placePointOnCap(&R[2*(Rdx-1)], &y2[0], xCenter, mesh.sqdelta, T.E, nu_a, nu_b, nu_c, orientation, Rdx, R);
                }
                // Append y2
                doubleVec_copyTo(&y2[0], &R[2*Rdx], 2);
                is_onEdge = true; // This point lies on the edge.
                Rdx += 1;
            }
        }
    }
    //[DEBUG]
    //(len(RD)>1) cares for the case that either the first and the last point lie on an endge
    // and there is no other point at all.
    //shift=1;
    if (is_onEdge && (!is_firstPointLiesOnVertex && Rdx > 1) && isPlacePointOnCap){
        Rdx += placePointOnCap(&R[2*(Rdx-1)], &R[0], xCenter, mesh.sqdelta, T.E, nu_a, nu_b, nu_c, orientation, Rdx, R);
    }

    // Construct List of Triangles from intersection points
    if (Rdx < 3){
        // In this case the content of the array out_RE will not be touched.
        return 0;
    } else {

        for (k=0; k < (Rdx - 2); k++){
            for (i=0; i<2; i++){
                // i is the index which runs first, then h (which does not exist here), then k
                // hence if we increase i, the *-index (of the pointer) inreases in the same way.
                // if we increase k, there is quite a 'jump'
                reTriangleList[2 * (3 * k + 0) + i] = R[i];
                reTriangleList[2 * (3 * k + 1) + i] = R[2 * (k + 1) + i];
                reTriangleList[2 * (3 * k + 2) + i] = R[2 * (k + 2) + i];
            }
        }
        // Excessing the bound out_Rdx will not lead to an error but simply to corrupted data!

        return Rdx - 2; // So that, it acutally contains the number of triangles in the retriangulation
    }
}
/*int method_retriangulateInfty(const double * xCenter, const double * TE,
                              double sqdelta, double * reTriangleList,
                              int isPlacePointOnCap) {
    cout << "ERROR: method_retriangulateInfty not implemented. Uncomment code in integration.cpp" << endl;
    abort();
    return 0;
}
*/
int method_retriangulateInfty(const double * xCenter, const double * TE,
                          double sqdelta, double * reTriangleList,
                          int isPlacePointOnCap){
    vector<Point> points;
    double delta = sqrt(sqdelta);
    double nu_a[2], nu_b[2], nu_c[2]; // Normals
    arma::vec p(2);
    arma::vec q(2);
    arma::vec a(2);
    arma::vec b(2);
    arma::vec y1(2);
    arma::vec y2(2);
    arma::vec vec_x_center(xCenter, 2);
    arma::vec xPoint(2);

    for (int k=0; k<3; k++) {
        int edgdx0 = k;
        int edgdx1 = (k + 1) % 3;

        doubleVec_copyTo(&TE[2 * edgdx0], &p[0], 2);
        doubleVec_copyTo(&TE[2 * edgdx1], &q[0], 2);

        if (vec_LInfdist(&TE[2 * edgdx0], xCenter, 2) <= delta){
            points.emplace_back(Point(TE[2 * edgdx0], TE[2 * edgdx0+1] ));
        }
        a = q - vec_x_center;
        b = p - q;
        int lambdaType[4];
        intVec_tozero(lambdaType, 4);
        double lambda[4];
        doubleVec_tozero(lambda, 4);

        if(abs(b[0]) > 0.){
            lambda[0] = (xCenter[0] + delta - q[0]) / b[0];
            lambda[1] = (xCenter[0] - delta - q[0]) / b[0];
        }
        if(abs(b[1]) > 0.){
            lambda[2] = (xCenter[1] + delta - q[1]) / b[1];
            lambda[3] = (xCenter[1] - delta - q[1]) / b[1];
        }

        for (double i : lambda){
            //cout << "i:" << i << ", lam " << lambda[i] << endl;
            if ((i > 0.) and (i < 1.)) {
                xPoint = q + b*i;
                if (vec_LInfdist(xPoint.memptr(), xCenter, 2) <= delta*(1+EPSILON)){
                    points.emplace_back(Point(xPoint[0], xPoint[1]));
                }
            }
        }

    }

    double lInfVert[2];
    int orientation = -signDet(TE);
    rightNormal(&TE[0], &TE[2], orientation, nu_a);
    rightNormal(&TE[2], &TE[4], orientation, nu_b);
    rightNormal(&TE[4], &TE[0], orientation, nu_c);

    lInfVert[0] = xCenter[0] + delta;
    lInfVert[1] = xCenter[1] + delta;
    if (inTriangle(lInfVert, &TE[0], &TE[2], &TE[4], nu_a, nu_b, nu_c)){
        points.emplace_back(Point(lInfVert[0], lInfVert[1]));
    }

    lInfVert[0] = xCenter[0] + delta;
    lInfVert[1] = xCenter[1] - delta;
    if (inTriangle(lInfVert, &TE[0], &TE[2], &TE[4], nu_a, nu_b, nu_c)){
        points.emplace_back(Point(lInfVert[0], lInfVert[1]));
    }

    lInfVert[0] = xCenter[0] - delta;
    lInfVert[1] = xCenter[1] - delta;
    if (inTriangle(lInfVert, &TE[0], &TE[2], &TE[4], nu_a, nu_b, nu_c)){
        points.emplace_back(Point(lInfVert[0], lInfVert[1]));
    }

    lInfVert[0] = xCenter[0] - delta;
    lInfVert[1] = xCenter[1] + delta;
    if (inTriangle(lInfVert, &TE[0], &TE[2], &TE[4], nu_a, nu_b, nu_c)){
        points.emplace_back(Point(lInfVert[0], lInfVert[1]));
    }

    //cout << "Set up points" << endl;
    Triangulation T;
    T.insert(points.begin(), points.end());
    int Rdx=0;
    //cout << "Set up Triangulation" << endl;
    for(Finite_face_iterator it = T.finite_faces_begin();
        it != T.finite_faces_end();
        ++it){
        for (int k=0; k<3; k++){
            Point p = it->vertex(k)->point();
            //std::cout << "x " << p.x() << ", ";
            //std::cout << "y " << p.y() << std::endl;
            reTriangleList[2 * (3 * Rdx + k) + 0] = p.x();
            reTriangleList[2 * (3 * Rdx + k) + 1] = p.y();
        }
        Rdx++;
    }
    //cout << "Success!" << endl;
    return Rdx;
}

int method_retriangulateInfty(const double * xCenter, const ElementType & T,
                              const MeshType & mesh, double * reTriangleList,
                              int isPlacePointOnCap){
    return method_retriangulateInfty(xCenter, T.E, mesh.sqdelta, reTriangleList, isPlacePointOnCap);
}

// Signature for debugging purpose only (Appears in Cassemble.h)
int method_exact(const double * xCenter, const double * TE,
                 const double sqdelta, double * reTriangleList, double * capsList, double * capsWeights,
                 int * prtnCaps){
    *prtnCaps = 0;
    // C Variables and Arrays.
    int i=0, k=0, edgdx0=0, edgdx1=0, Rdx=0;
    double v=0, lam1=0, lam2=0, term1=0, term2=0;
    double nu_a[2], nu_b[2], nu_c[2]; // Normals
    arma::vec p(2);
    arma::vec q(2);
    arma::vec a(2);
    arma::vec b(2);
    arma::vec y1(2);
    arma::vec y2(2);
    arma::vec vec_x_center(xCenter, 2);
    double orientation;
    bool is_onEdge=false, is_firstPointLiesOnVertex=true;
    // The upper bound for the number of required points is 9
    // Hence 9*2 is an upper bound to encode all resulting triangles
    // Hence we can hardcode how much space needs to bee allocated
    // (This upper bound is thight! Check Christian Vollmann's thesis for more information.)

    double R[9*2]; // Vector containing all intersection points.
    doubleVec_tozero(R, 9*2);

    // Compute Normals of the Triangle
    orientation = -signDet(TE);
    rightNormal(&TE[0], &TE[2], orientation, nu_a);
    rightNormal(&TE[2], &TE[4], orientation, nu_b);
    rightNormal(&TE[4], &TE[0], orientation, nu_c);

    for (k=0; k<3; k++){
        edgdx0 = k;
        edgdx1 = (k+1) % 3;

        doubleVec_copyTo(&TE[2*edgdx0], &p[0], 2);
        doubleVec_copyTo(&TE[2*edgdx1], &q[0], 2);

        a = q - vec_x_center;
        b = p - q;

        if (vec_sqL2dist(&p[0], xCenter, 2) <= sqdelta){
            doubleVec_copyTo(&p[0], &R[2*Rdx], 2);
            is_onEdge = false; // This point does not lie on the edge.
            Rdx += 1;
        }
        // PQ-Formula to solve quadratic problem
        v = pow( dot(a, b), 2) - (dot(a, a) - sqdelta) * dot(b, b);
        // If there is no sol to the quadratic problem, there is nothing to do.
        if (v >= 0){
            term1 = - dot(a, b) / dot(b, b);
            term2 = sqrt(v) / dot(b, b);

            // Vieta's Formula for computing the roots
            if (term1 > 0){
                lam1 = term1 + term2;
                lam2 = 1/lam1*(dot(a, a) - sqdelta) / dot(b, b);
            } else {
                lam2 = term1 - term2;
                lam1 = 1/lam2*(dot(a, a) - sqdelta) / dot(b, b);
            }
            y1 = lam1*b + q;
            y2 = lam2*b + q;

            // Check whether the first lambda "lies on the Triangle".
            if ((0 <= lam1) && (lam1 <= 1)){
                is_firstPointLiesOnVertex = is_firstPointLiesOnVertex && (bool)Rdx;

                // HERE CHANGE FOR EXACT CAPS --------------------------------------------------
                // Check whether the predecessor lied on the edge
                if (is_onEdge && (*prtnCaps < 3)){
                    //Rdx += placePointOnCap(&R[2*(Rdx-1)], &y1[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, Rdx, R);
                    capsWeights[*prtnCaps] = placePointCapCenter(&R[2*(Rdx-1)], &y1[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, &(capsList[(*prtnCaps)*2]));
                    *prtnCaps += !double_eq(capsWeights[*prtnCaps], 0.);
                }

                // Append y1
                doubleVec_copyTo(&y1[0], &R[2*Rdx], 2);
                is_onEdge = true; // This point lies on the edge.
                Rdx += 1;
            }
            // Check whether the second lambda "lies on the Triangle".
            if ((0 <= lam2) && (lam2 <= 1) && (scal_sqL2dist(lam1, lam2) > 0)){
                is_firstPointLiesOnVertex = is_firstPointLiesOnVertex && (bool)Rdx;

                // HERE CHANGE FOR EXACT CAPS --------------------------------------------------
                // Check whether the predecessor lied on the edge
                if (is_onEdge && (*prtnCaps < 3)){
                    //Rdx += placePointOnCap(&R[2*(Rdx-1)], &y1[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, Rdx, R);
                    capsWeights[*prtnCaps] = placePointCapCenter(&R[2*(Rdx-1)], &y2[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, &(capsList[(*prtnCaps)*2]));
                    *prtnCaps += !double_eq(capsWeights[*prtnCaps], 0.);
                }

                // Append y2
                doubleVec_copyTo(&y2[0], &R[2*Rdx], 2);
                is_onEdge = true; // This point lies on the edge.
                Rdx += 1;
            }
        }
    }
    //[DEBUG]
    //(len(RD)>1) cares for the case that either the first and the last point lie on an endge
    // and there is no other point at all.
    //shift=1;

    // HERE CHANGE FOR EXACT CAPS --------------------------------------------------
    if (is_onEdge && (*prtnCaps < 3)){
        //Rdx += placePointOnCap(&R[2*(Rdx-1)], &y1[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, Rdx, R);
        capsWeights[*prtnCaps] = placePointCapCenter(&R[2*(Rdx-1)], &R[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, &(capsList[(*prtnCaps)*2]));
        *prtnCaps += !double_eq(capsWeights[*prtnCaps], 0.);
    }

    // Construct List of Triangles from intersection points
    if (Rdx < 3){
        // In this case the content of the array out_RE will not be touched.
        return 0;
    } else {

        for (k=0; k < (Rdx - 2); k++){
            for (i=0; i<2; i++){
                // i is the index which runs first, then h (which does not exist here), then k
                // hence if we increase i, the *-index (of the pointer) inreases in the same way.
                // if we increase k, there is quite a 'jump'
                reTriangleList[2 * (3 * k + 0) + i] = R[i];
                reTriangleList[2 * (3 * k + 1) + i] = R[2 * (k + 1) + i];
                reTriangleList[2 * (3 * k + 2) + i] = R[2 * (k + 2) + i];
            }
        }
        // Excessing the bound out_Rdx will not lead to an error but simply to corrupted data!

        return Rdx - 2; // So that, it acutally contains the number of triangles in the retriangulation
    }
}


// Actual function signature (Appears in integration.h)
int method_exact(const double * xCenter, const ElementType & T,
                 const MeshType & mesh, double * reTriangleList, double * capsList, double * capsWeights,
                 int * prtnCaps){
    return method_exact(xCenter, T.E, mesh.sqdelta, reTriangleList, capsList, capsWeights, prtnCaps);
}


int method_retriangulate(const double * xCenter, const double * TE,
                         double sqdelta, double * reTriangleList,
                         int isPlacePointOnCap){
    // C Variables and Arrays.
    int i=0, k=0, edgdx0=0, edgdx1=0, Rdx=0;
    double v=0, lam1=0, lam2=0, term1=0, term2=0;
    double nu_a[2], nu_b[2], nu_c[2]; // Normals
    arma::vec p(2);
    arma::vec q(2);
    arma::vec a(2);
    arma::vec b(2);
    arma::vec y1(2);
    arma::vec y2(2);
    arma::vec vec_x_center(xCenter, 2);
    double orientation;

    bool is_onEdge=false, is_firstPointLiesOnVertex=true;
    // The upper bound for the number of required points is 9
    // Hence 9*2 is an upper bound to encode all resulting triangles
    // Hence we can hardcode how much space needs to bee allocated
    // (This upper bound is thight! Check Christian Vollmann's thesis for more information.)

    double R[9*2]; // Vector containing all intersection points.
    doubleVec_tozero(R, 9*2);

    // Compute Normals of the Triangle
    orientation = -signDet(TE);
    rightNormal(&TE[0], &TE[2], orientation, nu_a);
    rightNormal(&TE[2], &TE[4], orientation, nu_b);
    rightNormal(&TE[4], &TE[0], orientation, nu_c);

    for (k=0; k<3; k++){
        edgdx0 = k;
        edgdx1 = (k+1) % 3;

        doubleVec_copyTo(&TE[2*edgdx0], &p[0], 2);
        doubleVec_copyTo(&TE[2*edgdx1], &q[0], 2);

        a = q - vec_x_center;
        b = p - q;

        if (vec_sqL2dist(&p[0], xCenter, 2) <= sqdelta){
            doubleVec_copyTo(&p[0], &R[2*Rdx], 2);
            is_onEdge = false; // This point does not lie on the edge.
            Rdx += 1;
        }
        // PQ-Formula to solve quadratic problem
        v = pow( dot(a, b), 2) - (dot(a, a) - sqdelta) * dot(b, b);
        // If there is no sol to the quadratic problem, there is nothing to do.
        if (v >= 0){
            term1 = - dot(a, b) / dot(b, b);
            term2 = sqrt(v) / dot(b, b);

            // Vieta's Formula for computing the roots
            if (term1 > 0){
                lam1 = term1 + term2;
                lam2 = 1/lam1*(dot(a, a) - sqdelta) / dot(b, b);
            } else {
                lam2 = term1 - term2;
                lam1 = 1/lam2*(dot(a, a) - sqdelta) / dot(b, b);
            }
            y1 = lam1*b + q;
            y2 = lam2*b + q;

            // Check whether the first lambda "lies on the Triangle".
            if ((0 <= lam1) && (lam1 <= 1)){
                is_firstPointLiesOnVertex = is_firstPointLiesOnVertex && (bool)Rdx;
                // Check whether the predecessor lied on the edge
                if (is_onEdge && isPlacePointOnCap){
                    Rdx += placePointOnCap(&R[2*(Rdx-1)], &y1[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, Rdx, R);
                }
                // Append y1
                doubleVec_copyTo(&y1[0], &R[2*Rdx], 2);
                is_onEdge = true; // This point lies on the edge.
                Rdx += 1;
            }
            // Check whether the second lambda "lies on the Triangle".
            if ((0 <= lam2) && (lam2 <= 1) && (scal_sqL2dist(lam1, lam2) > 0)){
                is_firstPointLiesOnVertex = is_firstPointLiesOnVertex && (bool)Rdx;

                // Check whether the predecessor lied on the edge
                if (is_onEdge && isPlacePointOnCap){
                    Rdx += placePointOnCap(&R[2*(Rdx-1)], &y2[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, Rdx, R);
                }
                // Append y2
                doubleVec_copyTo(&y2[0], &R[2*Rdx], 2);
                is_onEdge = true; // This point lies on the edge.
                Rdx += 1;
            }
        }
    }
    //[DEBUG]
    //(len(RD)>1) cares for the case that either the first and the last point lie on an endge
    // and there is no other point at all.
    //shift=1;
    if (is_onEdge && (!is_firstPointLiesOnVertex && Rdx > 1) && isPlacePointOnCap){
        Rdx += placePointOnCap(&R[2*(Rdx-1)], &R[0], xCenter, sqdelta, TE, nu_a, nu_b, nu_c, orientation, Rdx, R);
    }

    // Construct List of Triangles from intersection points
    if (Rdx < 3){
        // In this case the content of the array out_RE will not be touched.
        return 0;
    } else {

        for (k=0; k < (Rdx - 2); k++){
            for (i=0; i<2; i++){
                // i is the index which runs first, then h (which does not exist here), then k
                // hence if we increase i, the *-index (of the pointer) inreases in the same way.
                // if we increase k, there is quite a 'jump'
                reTriangleList[2 * (3 * k + 0) + i] = R[i];
                reTriangleList[2 * (3 * k + 1) + i] = R[2 * (k + 1) + i];
                reTriangleList[2 * (3 * k + 2) + i] = R[2 * (k + 2) + i];
            }
        }
        // Excessing the bound out_Rdx will not lead to an error but simply to corrupted data!

        return Rdx - 2; // So that, it acutally contains the number of triangles in the retriangulation
    }
}


// Helpers Peridynamics ------------------------------------------------------------------------------------------------
void setupElement(const MeshType &mesh, const long * Vdx_new, ElementType &T){
    T.matE = arma::vec(mesh.dim*(mesh.dim+1));
    for (int k=0; k<mesh.dVertex; k++) {
        //Vdx = mesh.Triangles(k, Tdx);
        for (int j = 0; j < mesh.dim; j++) {
            T.matE[mesh.dim * k + j] = mesh.Verts(j, Vdx_new[k]);
            //printf ("aT %3.2f ", T.matE[mesh.dim * k + j]);
        }
    }

    // Initialize Structs
    T.E = T.matE.memptr();
    T.absDet = absDet(T.E, mesh.dim);
    T.signDet = static_cast<int>(signDet(T.E, mesh));
    T.dim = mesh.dim;
}
int join(const ElementType &aT, const ElementType &bT, const MeshType &mesh,
         ElementType &aTsorted, ElementType &bTsorted, int * argSortA, int * argSortB){
    //cout << "Welcome to join()" << endl;
    //TODO Try usage if STL Algorithms
    int nEqual = 0;
    int AinB[3], BinA[3];
    const long * aVdx = &(mesh.Triangles(0, aT.Tdx));
    const long * bVdx = &(mesh.Triangles(0, bT.Tdx));
    long aVdxsorted[3], bVdxsorted[3];

    intVec_tozero(AinB, 3);
    intVec_tozero(BinA, 3);

    for (int a=0; a<mesh.dVertex; a++){
        for (int b=0; b<mesh.dVertex; b++) {
            if (aVdx[a] == bVdx[b]) {
                AinB[a] += 1;
                BinA[b] += 1;
                aVdxsorted[nEqual] = aVdx[a];
                argSortA[nEqual] = a;
                bVdxsorted[nEqual] = bVdx[b];
                argSortB[nEqual] = b;
                nEqual += 1;
            }
        }
    }
    int ia = 0, ib = 0;
    for (int i=0; i<3; i++){
        if (!AinB[i]){
            aVdxsorted[nEqual + ia] = aVdx[i];
            argSortA[nEqual + ia] = i;
            ia++;
        }
        if (!BinA[i]){
            bVdxsorted[nEqual + ib] = bVdx[i];
            argSortB[nEqual + ib] = i;
            ib++;
        }
        //printf("%li, %li | %li, %li \n", aVdx[i], bVdx[i], aVdxsorted[i],  bVdxsorted[i] );
        //printf("%i, %i \n", argSortA[i], argSortB[i]);
    }
    setupElement(mesh, aVdxsorted, aTsorted);
    setupElement(mesh, bVdxsorted, bTsorted);
    //abort();
    return nEqual;
}



// [End] Helpers Peridynamics ---------------------------------------------------------------------------------------
#endif //NONLOCAL_ASSEMBLY_INTEGRATION_CPP
//
// Created by klar on 14.09.20.
//
#include "MeshTypes.h"
#include "mathhelpers.h"

// ### BASIS FUNCTION ##################################################################################################

void model_basisFunction(const double * p, double *psi_vals){
    psi_vals[0] = 1 - p[0] - p[1];
    psi_vals[1] = p[0];
    psi_vals[2] = p[1];
}
//TODO Is basically a OnePointFunction
// Should take Point and return an output
void model_basisFunction(const double * p, const int dim, double *psi_vals){
    int i=0;

    psi_vals[0] = 1;
    for (i=0; i<dim; i++){
        psi_vals[0] -= p[i];
        psi_vals[i+1] = p[i];
    }
}

void model_basisFunction_substracted(const double * alpha, const int dim, double *psi_vals){
    psi_vals[0] = -alpha[0] - alpha[1] + alpha[2] + alpha[3];
    psi_vals[1] = alpha[0] - alpha[2];
    psi_vals[2] = alpha[1] - alpha[3];
}

double traffoCommonVertex0(double * alpha){
    //xi, eta1, eta2, eta3 = alpha;
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi;
    alpha[1] = eta1 * xi;
    alpha[2] = eta2 * xi;
    alpha[3] = eta2 * eta3 * xi;
    return pow(xi,3)*eta2;
}
double traffoCommonVertex1(double * alpha){
    //xi, eta1, eta2, eta3 = alpha;
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];
    alpha[0] = xi*eta2;
    alpha[1] = xi*eta2*eta3;
    alpha[2] = xi;
    alpha[3] = xi*eta1;

    return pow(xi,3)*eta2;
}

double traffoCommonEdge0( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi;
    alpha[1] = xi*eta1*eta3;
    alpha[2] = xi*(1. - eta1*eta2);
    alpha[3] = xi*eta1*(1.-eta2);
    return pow(xi,3)*pow(eta1,2);
}

double traffoCommonEdge1( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi;
    alpha[1] = xi*eta1;
    alpha[2] = xi*(1. - eta1*eta2*eta3);
    alpha[3] = xi*eta1*eta2*(1.-eta3);
    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoCommonEdge2( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi*(1. - eta1*eta2);
    alpha[1] = xi*eta1*(1. - eta2);
    alpha[2] = xi;
    alpha[3] = xi*eta1*eta2*eta3;

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoCommonEdge3( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi*(1. - eta1*eta2*eta3);
    alpha[1] = xi*eta1*eta2*(1. - eta3);
    alpha[2] = xi;
    alpha[3] = xi*eta1;

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoCommonEdge4( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi*(1. - eta1*eta2*eta3);
    alpha[1] = xi*eta1*(1. - eta2*eta3);
    alpha[2] = xi;
    alpha[3] = xi*eta1*eta2;

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoIdentical0( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi;
    alpha[1] = xi*(1. - eta1 + eta1*eta2);
    alpha[2] = xi*(1. - eta1*eta2*eta3);
    alpha[3] = xi*(1. - eta1);

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoIdentical1( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi*(1. - eta1*eta2*eta3);
    alpha[1] = xi*(1. - eta1);
    alpha[2] = xi;
    alpha[3] = xi*(1. - eta1 + eta1* eta2);

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoIdentical2( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi;
    alpha[1] = xi*eta1*(1. - eta2 + eta2*eta3);
    alpha[2] = xi*(1. - eta1*eta2);
    alpha[3] = xi*eta1*(1. - eta2);

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoIdentical3( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi*(1. - eta1*eta2);
    alpha[1] = xi*eta1*(1. - eta2);
    alpha[2] = xi;
    alpha[3] = xi*eta1*(1. - eta2 + eta2*eta3);

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoIdentical4( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi*(1. - eta1*eta2*eta3);
    alpha[1] = xi*eta1*(1. - eta2*eta3);
    alpha[2] = xi;
    alpha[3] = xi*eta1*(1. - eta2);

    return pow(xi,3)*pow(eta1,2)*eta2;
}

double traffoIdentical5( double * alpha){
    double  xi = alpha[0],
            eta1 = alpha[1],
            eta2 = alpha[2],
            eta3 = alpha[3];

    alpha[0] = xi;
    alpha[1] = xi*eta1*(1. - eta2);
    alpha[2] = xi*(1. - eta1*eta2*eta3);
    alpha[3] = xi*eta1*(1. - eta2*eta3);

    return pow(xi,3)*pow(eta1,2)*eta2;
}

int getElement(ElementClass &element){
    printf("Dimension is %d \n", element.dim);
    return 0;
}

void initializeElement(const int Tdx, const MeshType & mesh, ElementType & T){
    /*
    Copy coordinates of Triange b to bTE.

     Tempalte of Triangle Point data.
     2D Case, a, b, c are the vertices of a triangle
     T.E -> | a1 | a2 | b1 | b2 | c1 | c2 |
     Hence, if one wants to put T.E into col major order matrix it would be of shape\
                                 | a1 | b1 | c1 |
     M(mesh.dim, mesh.dVerts) =  | a2 | b2 | c2 |
    */

    int j, k, Vdx;
    //T.matE = arma::vec(dim*(dim+1));
    for (k=0; k<mesh.dVertex; k++) {
        //Vdx = mesh.ptrTriangles[(mesh.dVertex+1)*Tdx + k+1];
        Vdx = mesh.Triangles(k, Tdx);
        for (j=0; j<mesh.dim; j++){
            T.matE[mesh.dim * k + j] = mesh.Verts(j, Vdx);
            //printf ("%3.2f ", T.matE[mesh.dim * k + j]);
            //T.matE[mesh.dim * k + j] = mesh.ptrVerts[ mesh.dim*Vdx + j];
        }
        //printf("\n");
    }
    // Initialize Struct
    T.E = T.matE.memptr();
    T.absDet = absDet(T.E, mesh.dim);
    T.signDet = static_cast<int>(signDet(T.E, mesh));
    T.label = mesh.LabelTriangles(Tdx);

    //T.label = mesh.ptrTriangles[(mesh.dVertex+1)*Tdx];
    //T.dim = dim;
    T.dim = mesh.dim;
    T.Tdx = Tdx;
}

void initializeQuadrule(QuadratureType & quadRule, const MeshType & mesh){
    std::list<double (*)(double *)> traffoList[3];
    traffoList[0] = traffoCommonVertex;
    traffoList[1] = traffoCommonEdge;
    traffoList[2] = traffoIdentical;
    const double SCALEDET = 0.0625;

    for(int nEqual_case=0; nEqual_case < 3; nEqual_case++){
        arma::Mat<double> * traffodetWeakCanceled;
        arma::Mat<double> * traffodetFractionalCanceled;
        arma::Cube<double> * alpha;
        arma::Cube<double> * alphaCanceled;

        if (nEqual_case == 0){
            traffodetWeakCanceled = &quadRule.traffodetWeakCanceled_CommonVertex;
            traffodetFractionalCanceled = &quadRule.traffodetFractionalCanceled_CommonVertex;
            alpha = &quadRule.alpha_CommonVertex;
            alphaCanceled = &quadRule.alphaCanceled_CommonVertex;

        } else if (nEqual_case == 1) {
            traffodetWeakCanceled = &quadRule.traffodetWeakCanceled_CommonEdge;
            traffodetFractionalCanceled = &quadRule.traffodetFractionalCanceled_CommonEdge;
            alpha = &quadRule.alpha_CommonEdge;
            alphaCanceled = &quadRule.alphaCanceled_CommonEdge;

        } else if (nEqual_case == 2) {
            traffodetWeakCanceled = &quadRule.traffodetWeakCanceled_Identical;
            traffodetFractionalCanceled = &quadRule.traffodetFractionalCanceled_Identical;
            alpha = &quadRule.alpha_Identical;
            alphaCanceled = &quadRule.alphaCanceled_Identical;

        }
        int traffoCounter=0;
        for(auto & traffo : traffoList[nEqual_case]) {
            //cout << "Traffo Number " << traffoCounter << ": ";
            for (int k = 0; k < quadRule.nPg; k++) {
                //cout << "   k: " << k << endl;
                for (int j = 0; j < 4; j++) {
                    // Copy and rescale quadrature points from tensor Gauss quadrature on the
                    // square (-1.,1)^2 -> (0,1)^2
                    (*alpha)(j, k, traffoCounter) = quadRule.Pg[4 * k + j];
                    //alpha[j + k*quadRule.nPg + traffoCounter*2*quadRule.nPg] = quadRule.Pg[4 * k + j]*0.5 + 0.5;

                    //cout << "      j : " << j << ",\t" << (*alpha)(j, k, traffoCounter) << endl;

                    (*alphaCanceled)(j, k, traffoCounter) = quadRule.Pg[4 * k + j];
                    //alphaCanceled[j + k*quadRule.nPg + traffoCounter*2*quadRule.nPg] = quadRule.Pg[4 * k + j]*0.5 + 0.5;
                }
                scale(&(*alpha)(0, k, traffoCounter));
                scale(&(*alphaCanceled)(0, k, traffoCounter));

                //cout << "traffodetCanceled [k, traffoC]" << (*traffodetCanceled)(k, traffoCounter) << endl;

                // Depending on the defree of the singularity we have to cancel out some parts of the
                // determinant. This happens here, based on the correct power of the variable xi = alpha[0] (before mirror(alpha))
                (*traffodetWeakCanceled)(k, traffoCounter)  = pow((*alphaCanceled)(0, k, traffoCounter), 3-2-2*mesh.fractional_s); // Det of rescaling from above.
                (*traffodetFractionalCanceled)(k, traffoCounter) = pow((*alphaCanceled)(0, k, traffoCounter), 3-2*mesh.fractional_s); // Det of rescaling from above.
                (*alphaCanceled)(0, k, traffoCounter) = 1.0;

                traffo(&(*alpha)(0, k, traffoCounter));
                double traffoDet = traffo(&(*alphaCanceled)(0, k, traffoCounter)) ;
                (*traffodetWeakCanceled)(k, traffoCounter) *= traffoDet * SCALEDET;
                (*traffodetFractionalCanceled)(k, traffoCounter) *= traffoDet * SCALEDET;

                //cout << "traffodetCanceled [k, traffoC]" << (*traffodetCanceled)(k, traffoCounter) << endl;

                mirror(&(*alpha)(0, k, traffoCounter));
                mirror(&(*alphaCanceled)(0, k, traffoCounter));
                //cout << "alpha [0, k, traffoC]" << (*alpha)(0, k, traffoCounter) << endl;
                //cout << "alphaCanceled [0, k, traffoC]" << (*alphaCanceled)(0, k, traffoCounter) << endl;
                //model_basisFunction(&(*alpha)(0, k, traffoCounter), 2, psix);
                //model_basisFunction(&(*alpha)(2, k, traffoCounter), 2, psiy);
            }
            traffoCounter++;
        }
    }

}
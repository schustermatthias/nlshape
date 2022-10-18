/**
    Helper functions.

    @file mathhelpers.cpp
    @author Manuel Klar
    @version 0.1 25/08/20
*/
#include <cmath>
#include <iostream>

#include "mathhelpers.h"
using namespace std;


// ___ MATH HELERPS IMPLEMENTATION _____________________________________________________________________________________

// Miscellaneous helpers ###############################################################################################

double evaluateZeta(const long * ptr_indices,
                    const long * ptr_indptr,
                    const long nColumns,
                    const long aTdx, const long bTdx){
    auto aBegin = ptr_indptr[aTdx];
    auto aStep = ptr_indptr[aTdx + 1] - aBegin;
    auto bBegin = ptr_indptr[bTdx];
    auto bStep = ptr_indptr[bTdx + 1] - bBegin;

    auto row_a = std::vector<long> (&ptr_indices[aBegin],
                                    &ptr_indices[aBegin] + aStep);
    auto row_b = std::vector<long> (&ptr_indices[bBegin],
                                    &ptr_indices[bBegin] + bStep);
    std::vector<long> intersection(nColumns);
    auto ls = std::set_intersection(row_a.begin(), row_a.end(),
                                    row_b.begin(), row_b.end(),
                                    intersection.begin());

    return (double) (ls -  intersection.begin());
}


void solve2x2(const double * A, const double * b, double * x){
    int dx0 = 0, dx1 = 1;
    double l=0, u=0;

    // Column Pivot Strategy
    if (absolute(A[0]) < absolute(A[2])){
        dx0 = 1;
        dx1 = 0;
    }

    // Check invertibility
    if (double_eq(A[2*dx0], 0)) {
        cout << "Error in solve2x2. Matrix not invertible." << endl;
        abort();
    }

    // LU Decomposition
    l = A[2*dx1]/A[2*dx0];
    u = A[2*dx1+1] - l*A[2*dx0+1];

    // Check invertibility
    if (u == 0){
        // raise LinAlgError("in solve2x2. Matrix not invertible.")
        cout << "in solve2x2. Matrix not invertible." << endl;
        abort();
    }

    // LU Solve
    x[1] = (b[dx1] - l*b[dx0])/u;
    x[0] = (b[dx0] - A[2*dx0+1]*x[1])/A[2*dx0];
}

// Normal which looks to the right w.r.t the vector from y0 to y1.
void rightNormal(const double * y0, const double * y1, const double orientation, double * normal){
    normal[0] = y1[1] - y0[1];
    normal[1] = y0[0] - y1[0];
    doubleVec_scale(orientation, normal, normal, 2);
}

int faculty(int n){
    int fac=1;
    for(int i=n; i>1; i--){
        fac*=i;
    }
    return fac;
}

void scale(double * alpha){
    for(int k=0; k<4; k++){
        alpha[k] = alpha[k]*0.5 + 0.5;
    }
}

void mirror(double * alpha){
    alpha[0] = alpha[0] - alpha[1];
    alpha[2] = alpha[2] - alpha[3];
}

// ### MATRIX OPERATIONS ###############################################################################################

double absDet(const double * E){
    double M[2][2];
    int i=0;
    for (i=0; i< 2; i++){
        M[i][0] = E[2*1+i] - E[2*0+i];
        M[i][1] = E[2*2+i] - E[2*0+i];
    }
    return absolute(M[0][0]*M[1][1] - M[0][1]*M[1][0]);
}
//TODO Check performance...
//TODO ARMADILLO_DEP
double absDet(const double * E, const int dim){
    // Let a,b,c,d be the Verticies of a Tetrahedon (3D)
    // Then M will be the 3x3 matrix containg [b-a,c-a,d-a]
    int dVertex = dim+1;

    arma::mat M(dim, dim, arma::fill::zeros);

    // Copy Values
    int i=0, j=0;
    for (i = 1; i < dVertex; i++) {
        for (j = 0; j < dim; j++) {
            M(j,i-1) += (E[j + dim*i] - E[j + 0]);
            //printf("%3.2f ", M(j,i-1));
        }
        //printf("\n");
    }
    return absolute(arma::det(M));
}

double signDet(const double * E){
    double M[2][2], det;
    int i=0;
    for (i=0; i< 2; i++){
        M[i][0] = E[2*1+i] - E[2*0+i];
        M[i][1] = E[2*2+i] - E[2*0+i];
    }
    det = (M[0][0]*M[1][1] - M[0][1]*M[1][0]);
    if (det > 0){
        return 1.;
    } else if ( det < 0){
        return -1.;
    } else {
        cout << "Warning in signDet(): Determinant is 0" << endl;
        return 0.0;
    }
}

double signDet(const double * E, const MeshType & mesh){
    // Let a,b,c,d be the Verticies of a Tetrahedon (3D)
    // Then M will be the 3x3 matrix containg [b-a,c-a,d-a]
    arma::mat M(mesh.dim, mesh.dim, arma::fill::zeros);
    double det;
    // Copy Values
    int i=0, j=0;
    for (i = 1; i < mesh.dVertex; i++) {
        for (j = 0; j < mesh.dim; j++) {
            M(j,i-1) += (E[j + mesh.dim*i] - E[j + 0]);
        }
    }
    det = arma::det(M);
    if (det > 0){
        return 1.;
    } else if ( det < 0){
        return -1.;
    } else {
        //cout << "Warning in signDet() 3D: Determinant is 0" << endl;
        return 0.0;
    }
}

void baryCenter(const double * E, double * bary){
    int i=0;
    bary[0] = 0;
    bary[1] = 0;
    for  (i=0; i< 3; i++){
        bary[0] += E[2*i+0];
        bary[1] += E[2*i+1];
    }
    bary[0] = bary[0]/3;
    bary[1] = bary[1]/3;
}

void baryCenter(const int dim, const double * E, double * bary){
    int i=0, j=0;
    int dVert = dim+1;
    doubleVec_tozero(bary, dim);

    for (j=0; j<dim; j++){
        for  (i=0; i<dVert; i++) {
            bary[j] += E[dim * i + j];
            //bary[1] += E[2*i+1];
        }
    bary[j] = bary[j]/ static_cast<double>(dVert);
    }
}

void baryCenter_polygone(const double * P, const int nVerticies, double * bary){
    int k=0;
    bary[0] = 0;
    bary[1] = 0;
    for (k=0; k<nVerticies; k++){
        bary[0] += P[2*k+0];
        bary[1] += P[2*k+1];
    }
    bary[0] = bary[0]/nVerticies;
    bary[1] = bary[1]/nVerticies;
}

void toPhys(const double * E, const double * p, double * out_x){
    int i=0;
    for (i=0; i<2;i++){
        out_x[i] = (E[2*1+i] - E[2*0+i])*p[0] + (E[2*2+i] - E[2*0+i])*p[1] + E[2*0+i];
    }
}

void toPhys(const double * E, const double * p, const int dim, double * out_x) {
    doubleVec_tozero(out_x, dim);
    for (int i=0; i<dim;i++){
        for(int j=0; j<dim;j++){
            out_x[i] += p[j]*(E[dim*(j+1)+i] - E[i]);
        }
        out_x[i] += E[i];
    }
}

void toRef(const double * E, const double * phys_x, double * ref_p){
    double M[2*2];
    double b[2];

    M[0] = E[2] - E[0];
    M[1] = E[4] - E[0];
    M[2] = E[3] - E[1];
    M[3] = E[5] - E[1];

    b[0] = phys_x[0] - E[0];
    b[1] = phys_x[1] - E[1];

    solve2x2(&M[0], &b[0], &ref_p[0]);
}

// get vector that contains the divergence of transformed basis functions(2D)
void get_div(const double * E, double * div){
    double gradient[6] = {-1.0, -1.0, 1.0, 0.0, 0.0, 1.0};
    double M[2*2];
    double inverse[2*2];
    double det;

    M[0] = E[2] - E[0];
    M[1] = E[4] - E[0];
    M[2] = E[3] - E[1];
    M[3] = E[5] - E[1];

    det = 1.0/(M[0]*M[3] - M[1]*M[2]);
    inverse[0] = det * M[3];
    inverse[1] = -det * M[1];
    inverse[2] = -det * M[2];
    inverse[3] = det * M[0];

    for (int i=0; i < 3; i++) {
        div[2*i] = gradient[2*i]*inverse[0] + gradient[2*i + 1]*inverse[2];
        div[2*i + 1] = gradient[2*i]*inverse[1] + gradient[2*i + 1]*inverse[3];
    }
}

// ### VECTOR OPERATIONS ###############################################################################################

// Check whether any, or all elements of a vector are zero -----------
int doubleVec_any(const double * vec, const int len){
    int i=0;
    for (i=0; i < len; i++){
        if (abs(vec[i]) > EPSILON){
            return 1;
        }
    }
    return 0;
}

double vec_dot(const double * x, const double * y, const int len){
    double r=0;
    auto ity = y;
    for (auto itx=x; itx < x+len; itx++){
        r += (*(itx)) * (*(ity));
        ity++;
    }
    return r;
}

double vec_sum(const double *x, const int len){
    double sum=0;
    for(int i=0; i<len; i++){
        sum+=x[i];
    }
    return sum;
}
double vec_sqL2dist(const double * x, const double * y, const int len){
    double r=0;
    int i=0;
    for (i=0; i<len; i++){
        r += pow((x[i] - y[i]), 2);
    }
    return r;
}

double vec_LInfdist(const double * x, const double * y, const int len){
    double r=0;
    int i=0;
    for (i=0; i<len; i++){
        r = max(absolute(x[i] - y[i]), r);
    }
    return r;
}

void doubleVec_tozero(double * vec, const int len){
    for (auto entry = vec; entry < vec+len; entry++){
        *entry  = 0.0;
    }
}

void doubleVec_midpoint(const double * vec1, const double * vec2, double * midpoint, const int len){
    int i = 0;
    for (i=0; i < len; i++){
        midpoint[i]  = (vec1[i] + vec2[i])/2;
    }
}

void doubleVec_subtract(const double * vec1, const double * vec2, double * out, const int len){
    int i=0;
    for (i=0; i < len; i++){
        out[i]  = vec1[i] - vec2[i];
    }
}

void doubleVec_add(const double * vec1, const double * vec2, double * out, const int len){
    int i=0;
    for (i=0; i < len; i++){
        out[i]  = vec1[i] + vec2[i];
    }
}

void doubleVec_scale(const double lambda, const double * vec, double * out, const int len){
    int i=0;
    for (i=0; i < len; i++){
        out[i]  = vec[i]*lambda;
    }
}

void doubleVec_copyTo(const double * input, double * output, const int len){
    int i=0;
    for (i=0; i<len; i++){
        output[i] = input[i];
    }
}
// Long

int longVec_all(const long * vec, const int len){
    int i=0;
    for (i=0; i<len; i++){
        if (vec[i] == 0){
            return 0;
        }
    }
    return 1;
}

int longVec_any(const long * vec, const int len){
    int i=0;
    for (i=0; i<len; i++){
        if (vec[i] != 0){
            return 1;
        }
    }
    return 0;
}

// Int

// Set Vectors to Zero -------------------------------------------------
void intVec_tozero(int * vec, const int len){
    for (auto it= vec; it < vec + len; it++){
        *it = 0;
    }
}
// Scalar --------------------------------------------------------

double absolute(const double value){
    if (value < 0){
        return - value;
    } else {
        return value;
    }
}

bool double_eq(double x, double y, const double eps){
    double diff = absolute(x-y);
    return (diff < eps);
}

double scal_sqL2dist(const double x, const double y){
    return pow((x-y), 2);
}

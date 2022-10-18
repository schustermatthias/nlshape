//
// Created by klar on 19.12.19.
//
#ifndef NONLOCAL_ASSEMBLY_MESHBUILDER_H
#define NONLOCAL_ASSEMBLY_MESHBUILDER_H
#include "armadillo"
#include "MeshTypes.h"
using namespace std;

static int max(int a, int b){
    if (a>=b){
        return a;
    }
    else {
        return b;
    }
}

static long factorial(const long n){
    long i=0, fac=1;
    if (n==0){
        return fac;
    }
    for (i=1; i<=n; i++){
        fac *= i;
    }
    return fac;
}

class Grid{
public:
    const long N_Omega;
    const double delta;
    const long dim;
    const int N_OmegaI = max(floor(delta*(N_Omega-1)) + 1, 2);

    Grid(long dim, long N_Omega, double delta):N_Omega(N_Omega), delta(delta), dim(dim){
        //cout << "Hello Class" << endl;
    }
    long clip(const long a, const long bound) const;
    bool inBounds(const long a, const long bound) const;
    // Guarantee that 0 and 1 are part of baseGrid
    // Create grid fom 0 to 1 with spacing 1/(N_Omega-1)
    const arma::vec omegaGrid = arma::linspace(0,1, N_Omega);
    // Create grid fom -delta to 0 with spacing of at most 1/(N_Omega-1)
    const arma::vec leftOmegaIGrid = arma::linspace(-delta, 0, N_OmegaI);
    // Create grid fom 1 to 1+delta with spacing of at most 1/(N_Omega-1)
    const arma::vec rightOmegaIGrid = arma::linspace(1, 1+delta, N_OmegaI);
    // Concatenate Grids
    const arma::vec baseGrid =  arma::join_cols(
                            arma::join_cols(
                                    leftOmegaIGrid.head(N_OmegaI-1),
                                    omegaGrid),
                                    rightOmegaIGrid.tail(N_OmegaI-1));

    const long N = baseGrid.n_elem; // Number of Nodes per Dimension
    const long M = N-1; // Number of 1D-Elements per Dimension

    void setVertexLabel(arma::vec & gridPoint) const;
    const int nV = pow(N, dim); // Number of Verticies
    const int nE = pow(M, dim) * factorial(dim); // Number of Elements
    long toVdx(arma::Col<long> vIndex) const;
    virtual arma::vec Vertex(const long Vdx) const = 0;
    virtual arma::Col<long> Element(const long i) const = 0;
    arma::vec DataVdx = arma::vec(nV, arma::fill::zeros);

};

class Grid1D : public Grid {
public:
    Grid1D(long N_Omega, double delta) :
    Grid(1, N_Omega, delta){
    }
private:
    arma::vec meshGrid(const long i) const;
    arma::Col<long> baseLine = arma::Col<long> {0, 1};
public:
    arma::vec Vertex(const long Vdx) const;
    arma::Col<long> Element(const long i) const;
    arma::vec Data = arma::vec(DataVdx.memptr(), nV, false, true);
    Grid1D refine() const;
};

class Grid2D : public Grid {
public:
    Grid2D(long N_Omega, double delta) :
            Grid(2, N_Omega, delta){
    }
    Grid2D(const Grid2D * grid); // Copy Constructor performs refinement!

private:
    arma::vec meshGrid(const long i, const long j) const;
    arma::Mat<long> baseIndexTriangle =
            arma::Mat<long> {
                    { 0,        0 },
                    { 1,        1 + N },
                    { 1 + N,    N }
            };
public:
    arma::vec Vertex(const long Vdx) const;
    arma::Col<long> Element(const long i) const;
    long ElementLabel(long Tdx) const;
    arma::Col<long> toijk(const long Vdx) const;
    arma::Col<long> Neighbour(long Tdx) const;
    arma::mat Data = arma::mat(DataVdx.memptr(), N,N, false, true);
    arma::vec refine() const;
    arma::uvec sortIndex; // initilialized in in mesh()
    arma::uvec invSortIndex;
    MeshType mesh(bool setNeighbours);
    QuadratureType quadrule() const;
    int save( string name);
};

class Grid3D : public Grid {
public:
    Grid3D(long N_Omega, double delta) :
            Grid(3, N_Omega, delta){

    }
private:
    arma::vec meshGrid(const long i, const long j, const long k) const;
    // The algorithm to find the tetrahedons in a cube is the same as in the case of the retriangulations.
    arma::Mat<long> baseIndexTetrahedon =
            arma::Mat<long> {
                    { 0,            0,            0,           0,           0,           0              },
                    { 1,            1+ N,         1 + N + N*N, 1,           1 + N*N,     1 + N + N*N    },
                    { 1 + N,        1 + N + N*N,  N,           1 + N*N,     1 + N + N*N, N*N            },
                    { 1 + N + N*N,  N,            N*N + N,     1 + N + N*N, N*N,         N*N + N        }
            };
public:
    arma::vec Vertex(const long Vdx) const;
    arma::Col<long> Element(const long i) const;
    long ElementLabel(long Tdx) const;
    arma::Col<long> toijk(long Vdx) const;
    arma::Col<long> Neighbour(long Tdx) const;
    //long toVdx(const long i, const long j, const long k) const;
    arma::cube Data = arma::cube(DataVdx.memptr(), N,N,N, false, true);
    Grid3D refine() const;
    arma::uvec sortIndex; // initilialized in in mesh()
    MeshType mesh();
    QuadratureType quadrule() const;
    int save();
};

#endif //NONLOCAL_ASSEMBLY_MESHBUILDER_H

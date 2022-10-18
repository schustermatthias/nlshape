//
// Created by klar on 16.03.20.
//
#ifndef NONLOCAL_ASSEMBLY_MESHTYPES_H
#define NONLOCAL_ASSEMBLY_MESHTYPES_H
#include "armadillo"
#include "cstring"
#include "metis.h"
#include <list>
using namespace std;

struct ElementStruct
/*!
 * Struct containing all necessary data of one finite element.
 *
 *  Template of Triangle Point data
 *
 *   2D Case, a, b, c are the vertices of a triangle
 *
*/
{
    /*!
     *
     *   T.E is of ordered the following way:
     *   | 0 | 1 | 2 | 3 | 4 | 5 |
     *   | -- | -- | -- | -- | -- | -- |
     *   | a1 | a2 | b1 | b2 | c1 | c2 |
     *
     *   Hence, if one wants to put T.E into cloumn major order matrix it would be of shape
     *   *M(mesh.dim, mesh.dVerts)* =
     *
     *    | 0   | 1   | 2   |
     *    | --- | --- | --- |
     *    | a1  | b1  | c1  |
     *    | a2  | b2  | c2  |
     */
    arma::vec matE;
    double * E;
    int dim;
    long label;
    double absDet;
    int signDet;
    long Tdx=0;
};
typedef ElementStruct ElementType;

class ElementClass {
    public:
        double * E = nullptr;
        int dim = 0;
        long label = -1;
        double absDet = 0.;
        int signDet = 0.;
        long Tdx=0;
        arma::vec matE;

        ElementClass(){
        };
        ElementClass(int dim_):
                dim(dim_){
            matE = arma::vec(this->dim*(dim+1), arma::fill::zeros);
            E = matE.memptr();
        };
        ~ElementClass () {};
};
int getElement(ElementClass &element);

struct ConfigurationStruct {
    const string path_spAd;
    const string path_fd;
    const string model_kernel;
    const string model_f;
    const string integration_method_remote;
    const string integration_method_close;
    const bool is_placePointOnCap;
    bool is_singularKernel;
    const bool is_fullConnectedComponentSearch;
    const int verbose;
    const int is_ShapeDerivative;
};
typedef ConfigurationStruct ConfigurationType;
//TODO This class is too complex and initialization is extremly error prone
struct MeshStruct{
    const int K_Omega;
    const int K;
    const long * ptrTriangles{};
    const long * ptrLabelTriangles{};
    const double * ptrVerts{};
    const long * ptrLabelVerts{};
    // Number of Triangles and number of Triangles in Omega
    const int nE;
    const int nE_Omega;
    // Number of vertices (in case of CG = K and K_Omega)
    const int nV;
    const int nV_Omega;
    const double delta;
    const double sqdelta;
    //const long * ptrNeighbours;
    //const int nNeighbours;

    const int is_DiscontinuousGalerkin;
    const int is_NeumannBoundary;

    const int dim;
    const int outdim;
    const int dVertex;

    // Weights for Domain decomposition (optional)
    const long * ptrZetaIndicator_indices{};
    const long * ptrZetaIndicator_indptr{};
    const long nZeta; // Should be set to 0
    // Optional Argument Mesh Diameter
    const double maxDiameter; // Should be set to 0 if unused.

    const double fractional_s;
    // This corresponds to a kernel with no singularity.
    // Therefore the integration for fractional kernels can be used for non-singular kernels as well
    // Makes sense for tests e.g.
    // TODO Write a wrapper for METIS
    idx_t *xadj{};
    idx_t *adjncy{};
    idx_t *eptr{};
    idx_t *eind{};

    const arma::Mat<double> Verts{arma::Mat<double>(this->ptrVerts, this->dim, this->nV)};
    //const arma::Mat<long> Neighbours{arma::Mat<long>(this->ptrNeighbours, this->nNeighbours, this->nE)};
    //TODO Refactor to Elements
    //TODO static vectors ...
    const arma::Mat<long> Triangles{arma::Mat<long>(this->ptrTriangles, this->dVertex, this->nE)};
    // Label of Triangles inside Omega = 1
    // Label of Triangles in OmegaI = 2
    const arma::Col<long> LabelTriangles{arma::Col<long>(this->ptrLabelTriangles, this->nE)};
    // TODO Struct of Points mit static vectors Verts, LabelVerts
    const arma::Col<long> LabelVerts{arma::Col<long>(this->ptrLabelVerts, this->nV)};

};
typedef MeshStruct MeshType;


struct QuadratureStruct{
    // Quadrature rule for the non-singular case
    //TODO Should contain Points, Weights
    const double * Px;
    const double * Py;
    const double * dx;
    const double * dy;

    const int nPx;
    const int nPy;
    const int dim;

    // Tensor Gauss-quadrature rule for regularizing integral transforms
    const double * Pg;
    const double * dg;
    const int tensorGaussDegree;
    const int nPg = pow(tensorGaussDegree, dim * 2);

    //Precomputed values of ansatz functions (non-singular case)
    arma::Mat<double> psix{arma::Mat<double>(this->dim +1, this->nPx)};
    arma::Mat<double> psiy{arma::Mat<double>(this->dim +1, this->nPy)};

    // Precomputed values for regularizing integral transforms

    arma::Mat<double> traffodetFractionalCanceled_CommonVertex{arma::Mat<double>( this->nPg, 2)};
    arma::Mat<double> traffodetFractionalCanceled_CommonEdge{arma::Mat<double>( this->nPg, 5)};
    arma::Mat<double> traffodetFractionalCanceled_Identical{arma::Mat<double>( this->nPg, 6)};

    arma::Mat<double> traffodetWeakCanceled_CommonVertex{arma::Mat<double>( this->nPg, 2)};
    arma::Mat<double> traffodetWeakCanceled_CommonEdge{arma::Mat<double>( this->nPg, 5)};
    arma::Mat<double> traffodetWeakCanceled_Identical{arma::Mat<double>( this->nPg, 6)};

    arma::Cube<double> alpha_CommonVertex{arma::Cube<double>(4, this->nPg, 2)};
    arma::Cube<double> alpha_CommonEdge{arma::Cube<double>(4, this->nPg, 5)};
    arma::Cube<double> alpha_Identical{arma::Cube<double>(4, this->nPg, 6)};

    arma::Cube<double> alphaCanceled_CommonVertex{arma::Cube<double>(4, this->nPg, 2)};
    arma::Cube<double> alphaCanceled_CommonEdge{arma::Cube<double>(4, this->nPg, 5)};
    arma::Cube<double> alphaCanceled_Identical{arma::Cube<double>(4, this->nPg, 6)};

    //arma::Cube<double> psialpha{arma::Cube<double>(this->dim , this->nPg, 3)};
    //arma::Cube<double> psialphaCenceled{arma::Cube<double>(this->dim , this->nPg, 3)};
};
typedef QuadratureStruct QuadratureType;

struct entryStruct{
    unsigned long dx;
    double value;

    bool operator<(const entryStruct &other) const{
        return this->dx < other.dx;
    }
    bool operator>(const entryStruct &other) const{
        return this->dx > other.dx;
    }
    bool operator==(const entryStruct &other) const{
        return this->dx == other.dx;
    }
};
typedef entryStruct entryType;

/**
 * @brief Definition of basis function (deprecated).
 * @param p Quadrature point
 * @param psi_vals Value of 3 basis functions.
 */
void model_basisFunction(const double * p, double *psi_vals);

/**
 * @brief  Definition of basis function.
 * @param p Quadrature point
 * @param dim Dimension of domain (2,3).
 * @param psi_vals  Value of 3 or 4 basis functions, depending on the dimension.
 */
void model_basisFunction(const double * p, int dim, double *psi_vals);
/**
 * @brief  Definition of basis function terms as they appear in the representation necessary for
 * the fractional laplacian.
 *
 * @param alpha Quadrature point of dim 4.
 * @param dim Dimension of domain (2).
 * @param psi_vals  Value of 3 basis functions.
 */
void model_basisFunction_substracted(const double * alpha, const int dim, double *psi_vals);

double traffoCommonVertex0(double * alpha);
double traffoCommonVertex1(double * alpha);

double traffoCommonEdge0( double * alpha);
double traffoCommonEdge1(double * alpha);
double traffoCommonEdge2( double * alpha);
double traffoCommonEdge3( double * alpha);
double traffoCommonEdge4( double * alpha);

double traffoIdentical0( double * alpha);
double traffoIdentical1( double * alpha);
double traffoIdentical2( double * alpha);
double traffoIdentical3( double * alpha);
double traffoIdentical4( double * alpha);
double traffoIdentical5( double * alpha);

const std::list<double(*)(double *)> traffoCommonVertex = {traffoCommonVertex0,
                                                           traffoCommonVertex1};
const std::list<double(*)(double *)> traffoCommonEdge = {traffoCommonEdge0,
                                                         traffoCommonEdge1,
                                                         traffoCommonEdge2,
                                                         traffoCommonEdge3,
                                                         traffoCommonEdge4};
const std::list<double(*)(double *)> traffoIdentical = {traffoIdentical0,
                                                        traffoIdentical1,
                                                        traffoIdentical2,
                                                        traffoIdentical3,
                                                        traffoIdentical4,
                                                        traffoIdentical5};
void initializeElement(int Tdx, const MeshType & mesh, ElementType & T);
//TODO It is good, that this is a free function, but the initialization should be inside of Quadrule somewhere
void initializeQuadrule(QuadratureType & quadRule, const MeshType & mesh);
#endif //NONLOCAL_ASSEMBLY_MESHTYPES_H

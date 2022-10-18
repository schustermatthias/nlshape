//
// Created by klar on 19.12.19.
//

#include "MeshBuilder.h"
#include "iostream"
#include "armadillo"
using namespace std;

arma::Col<long> Grid2D::toijk(const long Vdx) const {
    arma::Col<long> vIndex(2);
    // Think in Column Major order.
    // 2D Case - extract Cols
    vIndex(1) = Vdx / N;
    // 1D Case - extract Row
    vIndex(0) =  Vdx % N;;
    // Prune to correct dimension
    return vIndex;
}
arma::Col<long> Grid3D::toijk(long Vdx) const {
    arma::Col<long> vIndex(3);
    // Think in Column Major order.
    // i.e. The indices jump first in Row then Col then Tube.
    // 3D Case - extract Tubes
    vIndex(2) = Vdx / (long) pow(N,2);
    Vdx = Vdx % (long) pow(N,2);
    // 2D Case - extract Cols
    vIndex(1) = Vdx / N;
    Vdx = Vdx % N;
    // 1D Case - extract Row
    vIndex(0) = Vdx;
    // Prune to correct dimension
    return vIndex;
}

long Grid::toVdx(arma::Col<long> vIndex) const {
    long d = 0, Vdx = 0;
    for (d = 0; d < dim; d++) {
        Vdx += vIndex(d) * ((long) pow(N, d));
    }
    return Vdx;
}
/*
long Grid2D::toVdx(const long i, const long j) const {
    return Grid::toVdx(arma::Col<long> {i,j});
}

long Grid3D::toVdx(const long i, const long j, const long k) const {
    return Grid::toVdx(arma::Col<long> {i,j,k});
}
*/
void  Grid::setVertexLabel(arma::vec & gridPoint) const {
    if (arma::all(gridPoint.head(dim) > 1e-12) && arma::all(gridPoint.head(dim) < 1- 1e-12)){
        // Label of Omega is 1.
        gridPoint(dim) = 1;
    } else {
        // Label of OmegaI is 2.
        gridPoint(dim) = 2;
    }
}
arma::vec Grid1D::meshGrid(const long i) const {
    arma::vec vertexi(1);
    vertexi(0) = baseGrid(i);
    return vertexi;
}
arma::vec Grid2D::meshGrid(const long i, const long j)  const {
    arma::vec vertexij(2);
    vertexij(0) = baseGrid(i);
    vertexij(1) = baseGrid(j);
    return vertexij;
}
arma::vec Grid3D::meshGrid(const long i, const long j, const long k)  const {
    arma::vec vertexijk(3);
    vertexijk(0) = baseGrid(i);
    vertexijk(1) = baseGrid(j);
    vertexijk(2) = baseGrid(k);
    return vertexijk;
}
/**
* @brief Returns grid point Vdx and its label. Label consquently is of type double!
* @param Vdx Index of grid point.
* @return arma::vec, length dim + 1 [grid point, label]
*/
arma::vec Grid1D::Vertex(long Vdx) const {
    arma::vec gridPoint(2, arma::fill::zeros);
    gridPoint.head(1) = meshGrid(Vdx);
    setVertexLabel(gridPoint);
    return gridPoint;
}

arma::vec Grid2D::Vertex(long Vdx) const {
    arma::Col<long> vIndex(2);
    arma::vec gridPoint(3, arma::fill::zeros);
    if(Vdx <0 || Vdx > nV){
        cout << "Error in Grid::Vertex, Vertex does not exist." << endl;
        abort();
    }
    vIndex = toijk(Vdx);
    gridPoint.head(2) = meshGrid(vIndex(0), vIndex(1));
    setVertexLabel(gridPoint);
    return gridPoint;
}
arma::vec Grid3D::Vertex(long Vdx) const {
    arma::Col<long> vIndex(3);
    arma::vec gridPoint(4, arma::fill::zeros);
    if(Vdx <0 || Vdx > nV){
        cout << "Error in Grid::Vertex, Vertex does not exist." << endl;
        abort();
    }
    vIndex = toijk(Vdx);
    gridPoint.head(3) = meshGrid(vIndex(0), vIndex(1), vIndex(2));
    setVertexLabel(gridPoint);
    return gridPoint;
}

arma::Col<long> Grid1D::Element(const long Ldx) const {
    if (Ldx >= nE){
        cout << "Error in Grid::Line, Element does not exist." << endl;
        abort();
    }
    return baseLine + Ldx;
}

arma::Col<long> Grid2D::Element(long Tdx) const {
    /*arma::Col <long> triangleIndex(3);
    long cornerVdx=0;
    long trNumber; // 0, 1

    if (Tdx >= nE){
        cout << "Error in Grid::Triangle, Element does not exist." << endl;
        abort();
    }
    triangleIndex(2) = Tdx / (2*M);
    Tdx = Tdx % (2*M);
    triangleIndex(1) = Tdx / 2;
    Tdx = Tdx % 2;
    triangleIndex(0) = Tdx;

    // Corner Vertex = triangleIndex(1), triangleIndex(2)
    cornerVdx = triangleIndex(2)*N + triangleIndex(1);
    trNumber = triangleIndex(0);
    arma::Col<long> Vdx = baseIndexTriangle.col(trNumber);
    Vdx += cornerVdx;
    return  Vdx;
*/
    long trdx2, trdx1, cornerVdx;

    trdx2 = Tdx / (2*M);
    Tdx = Tdx % (2*M);
    trdx1 =  Tdx / 2;
    cornerVdx = trdx2*N + trdx1;

    if (Tdx % 2) {
        return arma::Col<long>({cornerVdx, cornerVdx+1 + N, cornerVdx + N});
    }
    else {
        return arma::Col<long> ({cornerVdx, cornerVdx+1, cornerVdx+1 + N});
    }
}

bool Grid::inBounds(const long a, const long bound) const{
    return ((a<bound) && (a >= 0));
}

arma::Col<long> Grid2D::Neighbour(long Tdx) const{
    arma::Col <long> neighbourIndex(3, arma::fill::ones);
    neighbourIndex.fill(nE);
    arma::Col <long> aT=Element(Tdx);
    arma::Col <long> bT(3);
    int bTdx, i, sum, k=0;

    for (bTdx=0; bTdx<nE; bTdx++){
        bT = Element(bTdx);
        sum=0;
        for (i=0; i<9; i++){
            sum += (aT(i/3) == bT(i%3));
        }
        if (sum == 2){
            neighbourIndex(k) = bTdx;
            k++;
        }
        if (k == 3){
            bTdx = nE;
        }
    }
    return neighbourIndex;
}

arma::Col<long> Grid3D::Neighbour(long Tdx) const{
    cout << "Not Thread Safe!" << endl;
    abort();
    arma::Col <long> neighbourIndex(4, arma::fill::ones);
    neighbourIndex.fill(nE);
    arma::Col <long> aT=Element(Tdx);
    arma::Col <long> bT(4);
    int bTdx, i, sum, k=0;

    for (bTdx=0; bTdx<nE; bTdx++){
        bT = Element(bTdx);
        sum=0;
        for (i=0; i<16; i++){
            sum += aT(i/4) == bT(i%4);
        }
        if (sum == 3){
            neighbourIndex(k) = bTdx;
            k++;
        }
        if (k == 4){
            bTdx = nE;
        }
    }
    return neighbourIndex;
}

arma::Col<long> Grid3D::Element(long Tdx) const {
    arma::Col <long> tetrahedonIndex(4);
    arma::Col<long> Vdx(4);
    long cornerVdx=0;
    long tetNumber; // 0 to 5

    if (Tdx >= nE){
        cout << "Error in Grid::Triangle, Element does not exist." << endl;
        abort();
    }
    tetrahedonIndex(3) = Tdx / (6*M*M);
    Tdx = Tdx % (6*M*M);
    tetrahedonIndex(2) = Tdx / (6*M);
    Tdx = Tdx % (6*M);
    tetrahedonIndex(1) = Tdx / 6;
    Tdx = Tdx % 6;
    tetrahedonIndex(0) = Tdx;

    // Corner Vertex = tetrahedonIndex(3), tetrahedonIndex(2), tetrahedonIndex(1)
    cornerVdx = tetrahedonIndex(3)*N*N + tetrahedonIndex(2)*N + tetrahedonIndex(1);
    tetNumber = tetrahedonIndex(0);
    // Shift baseIndexTetrahedon with number tetNumber to corner Vertex.
    Vdx = baseIndexTetrahedon.col(tetNumber);
    Vdx += cornerVdx;
    return  Vdx;
}

long Grid2D::ElementLabel(long Tdx) const {
    arma::Col<long> Vdx = Element(Tdx);
    int k;
    for (k=0; k<3; k++){
        if (1 == Vertex(Vdx(k))(2)){
            return 1;
        }
    }
    // Label for Triangles in Omega is 1, for Triangles in OmegaI is 2.
    // If there is any vertex in Omega then the trianlge is in Omega.
    return  2;
}

long Grid3D::ElementLabel(long Tdx) const {
    arma::Col<long> Vdx = Element(Tdx);
    arma::Col <long> vertexLabels(4);
    int k;
    for (k=0; k<4; k++){
        vertexLabels(k) = (long) Vertex(Vdx(k))(3);
    }
    // Label for Triangles in Omega is 1, for Triangles in OmegaI is 2.
    // If there is any vertex in Omega then the trianlge is in Omega.
    return  vertexLabels.min();
}

int Grid2D::save( string name ){
    MeshType mesh = Grid2D::mesh(false);
    //mesh.Neighbours.save("data/"+to_string(N)+name+"nbr", arma::arma_ascii);
    mesh.Verts.save("data/"+to_string(N)+name+"vrt", arma::arma_ascii);
    mesh.Triangles.save("data/"+to_string(N)+name+"tri", arma::arma_ascii);
    mesh.LabelTriangles.save("data/"+to_string(N)+name+"lab", arma::arma_ascii);
    // As Mesh is a sorted version of the points, DataVdx has to be sorted accordingly!
    arma::vec orderedData = DataVdx(sortIndex);
    orderedData.save("data/"+to_string(N)+name+"dat", arma::arma_ascii);
    return 0;
}

Grid1D Grid1D::refine() const{
    long i;
    Grid1D fineGrid(N_Omega*2-1, delta);
    if(fineGrid.N_OmegaI !=  2*(this->N_OmegaI)-1 ){
        cout << endl << "ERROR in Grid1D::refine: Number of points did not double in OmegaI. Please choose smaller grid size." << endl;
        abort();
    }
    // Scale
    for (i=0; i<N; i++){
        fineGrid.DataVdx(2*i) = DataVdx(i);
    }
    // Interpolate 1D
    for (i=0; i<N-1; i++){
        fineGrid.DataVdx(2*i+1) = (fineGrid.DataVdx(2*i) + fineGrid.DataVdx(2*i+2))/2;
    }
    return fineGrid;
}

Grid2D:: Grid2D(const Grid2D * coarseGrid):
Grid(2, 2*coarseGrid->N_Omega-1, coarseGrid->delta)
{
    //cout << "Note: Copy Constructor performs refinement." << endl;
    long i, j;
    if(N_OmegaI !=  2*(coarseGrid->N_OmegaI)-1 ){
        cout << endl << "ERROR in Grid2D::refine: Number of points did not double in OmegaI. Please choose smaller grid size!" << endl;
        abort();
    }

    // Scale
    for (i=0; i<coarseGrid->N; i++) {
        for (j = 0; j < coarseGrid->N; j++) {
            Data(2*i, 2*j) = coarseGrid->Data(i, j);
        }
    }

    // Interpolate 1D
    for (i=0; i<coarseGrid->N-1; i++){
        for (j = 0; j<coarseGrid->N; j++) {
            Data(2*i+1, 2*j) = (
                                                Data(2*i, 2*j) +
                                                Data(2*i + 2, 2*j)
                                        ) / 2;
        }
    }

    // Interpolate 2D
    for (i=0; i<2*coarseGrid->N-1; i++){
        for (j = 0; j<coarseGrid->N-1; j++) {
            Data(i, 2*j+1) = (
                                              Data(i, 2*j) +
                                              Data(i, 2*j+2)
                                      ) / 2;
        }
    }
}

arma::vec Grid2D::refine() const{
    long i, j;
    Grid2D fineGrid(N_Omega*2-1, delta);
    if(fineGrid.N_OmegaI !=  2*(this->N_OmegaI)-1 ){
        cout << endl << "ERROR in Grid1D::refine: Number of points did not double in OmegaI. Please choose smaller grid size." << endl;
        abort();
    }
    // Scale
    for (i=0; i<N; i++) {
        for (j = 0; j < N; j++) {
            fineGrid.Data(2*i, 2*j) = Data(i, j);
        }
    }

    // Interpolate 1D
    for (i=0; i<N-1; i++){
        for (j = 0; j<N; j++) {
            fineGrid.Data(2*i+1, 2*j) = (
                    fineGrid.Data(2*i, 2*j) +
                    fineGrid.Data(2*i + 2, 2*j)
                    ) / 2;
        }
    }

    // Interpolate 2D
    for (i=0; i<2*N-1; i++){
        for (j = 0; j<N-1; j++) {
            fineGrid.Data(i, 2*j+1) = (
                      fineGrid.Data(i, 2*j) +
                      fineGrid.Data(i, 2*j+2)
                                   ) / 2;
        }
    }
    return fineGrid.DataVdx;
}

Grid3D Grid3D::refine() const{
    long i,j,k, Nfine = 2*N_Omega-1;
    long counter, endLoop;
    int coarsenV = nV;
    Grid3D fineGrid(Nfine, delta);
    if(fineGrid.N_OmegaI !=  2*(this->N_OmegaI)-1 ){
        cout << endl << "ERROR in Grid1D::refine: Number of points did not double in OmegaI. Please choose compatible grid size." << endl;
        abort();
    }
    //arma::Cube<double> fineData(fineGrid.DataVdx.memptr(), Nfine, Nfine, Nfine, false, true);

    // Scale
    i=0; j=-1; k=-1;
    for (counter=0; counter<coarsenV; counter++) {
        // Think in Column-Major order.
        // see https://stackoverflow.com/questions/12710138/c-fastest-loop-scan-over-3d-vector
        // i is a row index, so i is between 0,...,N-1
        i = counter % N;
        // j is increased whenever i is 0, but we have to start with j=-1.
        j += !i;
        // j is a Column index, so between 0,...,N-1
        j %= N;
        // k is increased whenever i and j are 0, but we have to start with k=-1.
        k += !j*!i;
        fineGrid.Data(2 * i, 2 * j, 2 * k) = Data(i, j, k);

    }

    // Interpolate 1D
    i=0; j=-1; k=-1;
    endLoop = (N-1)*N*N;
    for (counter=0; counter<endLoop; counter++) {
        // Get index ijk;
        // We send our loop over every every spacing between columns i (N-1) and every second column j (N) and slice k (N).
        i = counter % (N-1); j += !i; j %= N; k += !j*!i;
        // Interpolate
        fineGrid.Data(2 * i + 1, 2 * j, 2 * k) = (
                 fineGrid.Data(2 * i, 2 * j, 2 * k) +
                 fineGrid.Data(2 * i + 2, 2 * j, 2 * k)
                                            ) / 2;
    }

    // Interpolate 2D
    i=0; j=-1; k=-1;
    // See loop description.
    endLoop = Nfine*(N-1)*N;
    for (counter=0; counter<endLoop; counter++) {
        // Get index ijk;
        // We send our loop over every Row i (Nfine), every spacing between columns j (N-1) and every second slice k (N).
        i = counter % Nfine; j += !i; j %= (N-1); k += !j*!i;
        // Interpolate
        fineGrid.Data(i, 2 * j + 1, 2*k) = (
                               fineGrid.Data(i, 2 * j, 2*k) +
                               fineGrid.Data(i, 2 * j + 2, 2*k)
                                 ) / 2;
    }

    // Interpolate 3D
    i=0; j=-1; k=-1;
    // See loop description.
    endLoop = Nfine*Nfine*(N-1);
    for (counter=0; counter<endLoop; counter++) {
        // Get index ijk;
        // We send our loop over every Row i (Nfine), every column j (Nfine) and the spacing between slices k (N-1).
        i = counter % Nfine; j += !i; j %= Nfine; k += !j*!i;
        // Interpolate
        fineGrid.Data(i, j, 2*k+1) = (fineGrid.Data(i, j, 2*k) + fineGrid.Data(i, j, 2*k+2)) / 2;
    }
    return fineGrid;
}

static arma::uvec invert_permutation(arma::uvec p) {
    int i;
    int n = p.n_elem;
    arma::uvec s(n);
    for (i=0; i<n; i++){
        s(p(i)) = i;
    }
    return s;
}

MeshType Grid2D::mesh(bool setNeighbours=true) {
    arma::Col<long> indices(3);
    arma::Mat<double> Verts(dim, nV);
    arma::Col<long> VertLabels(nV);
    arma::Mat<long> Elements(dim+1, nE);
    arma::Col<long> LabelElements(nE);
    arma::Mat<long> Neighbours(dim+1, nE);
    Neighbours.fill(nE);

    int i,k, nE_Omega, nV_Omega;
    for (i=0; i<nV; i++){
        Verts.col(i) = (Vertex(i)).head(2);
        VertLabels(i) = (Vertex(i))(2);
    }
    for (i=0; i<nE; i++){
        Elements.col(i) = Element(i);
        LabelElements(i) = ElementLabel(i);
    }

    if(setNeighbours){
        for (i=0; i<nE; i++) {
            Neighbours.col(i) = Neighbour(i);
        }
    }

    //arma::uvec sortIndex;
    sortIndex = arma::sort_index(VertLabels);
    invSortIndex =  invert_permutation(sortIndex);
    nE_Omega = 0;
    for (i=0; i< nE; i++){
        nE_Omega += (LabelElements(i) ==1);
    }
    nV_Omega = 0;
    for (i=0; i< nV; i++){
        nV_Omega += (VertLabels(i) == 1);
    }

    // Permute Vertices
    Verts = Verts.cols(sortIndex);
    VertLabels = VertLabels.elem(sortIndex);
    // Change entries in Elements accordingly (inverse permutation)
    for (i=0; i<nE; i++){
        indices = Elements.col(i);
        for (k=0; k<3; k++){
            Elements(k, i) = invSortIndex(Elements(k, i));
        }
    }

    //Permute Triangles as well
    arma::uvec sortIndexElements = arma::sort_index(LabelElements);
    arma::uvec invIndexElements =  invert_permutation(sortIndexElements);
    Elements = Elements.cols(sortIndexElements);
    LabelElements = LabelElements.elem(sortIndexElements);
    Neighbours = Neighbours.cols(sortIndexElements);
    // Both, the order and the entries of Neighbours, have to be changed accordingly
    for (i=0; i<nE; i++){
        for (k=0; k<3; k++){
            if (Neighbours(k, i) < nE) {
                Neighbours(k, i) = invIndexElements(Neighbours(k, i));
            }
        }
    }

    MeshType mesh = {nV_Omega, nV,
                     Elements.memptr(),
                     LabelElements.memptr(),
                     Verts.memptr(),
                     nE, nE_Omega,
                     nV, nV_Omega,
                     pow(delta, 2),
                     Neighbours.memptr(),
                     false,
                     false,
                     2, 3};
    return mesh;
}


QuadratureType Grid2D::quadrule() const{
    arma::vec dx, dy;
    dx.load("conf/dx16.txt", arma::arma_ascii);
    arma::mat Px, Py;
    Px.load("conf/dx16.txt", arma::arma_ascii);
    QuadratureType quadRule = {Px.memptr(),
                               Px.memptr(),
                               dx.memptr(),
                               dx.memptr(),
                               (int) Px.n_rows,
                               (int) Px.n_rows,
                               2};
    return quadRule;
}
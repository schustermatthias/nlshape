/**
    Contains assembly algorithm for nonlocal stiffnes matrix and forcing function.
    @file Cassemble.cpp
    @author Manuel Klar
    @version 0.1 25/08/20
*/

#include <iostream>
#include <cmath>
#include <queue>
#include <armadillo>
#include <map>
#include "metis.h"
#include <Cassemble.h>
#include <omp.h>
#include "integration.h"
#include "mathhelpers.h"
#include "model.h"
#include "checks.cpp"
using namespace std;
/**
 * This function looks up the configuration. It has to be updated whenever a new kernel,
 * forcing function or integration routine is added in order to make the option available.
 *
 * @param conf ConfigurationType
 * @param verbose
 */
void lookup_configuration(ConfigurationType & conf, int verbose=0){
    //TODO Should be in model kernel and model forcing

    // Lookup right hand side ------------------------------------------------------------------------------------------
    if (verbose) cout << "Right hand side: " << conf.model_f << endl;
    //void (*model_f)(const double * x, double * forcing_out);
    map<string, void (*)(const double * x, double * forcing_out)> lookup_f = {
            {"linear", f_linear},
            {"gaussian", f_gaussian},
            {"jump", f_jump},
            {"linear1D", f_linear1D},
            {"linear3D", f_linear3D},
            {"constant", f_constant},
            {"tensorsin", f_tensorsin},
            {"linearField", fField_linear},
            {"linearField3D", fField_linear3D},
            {"constantRightField", fField_constantRight},
            {"constantDownField", fField_constantDown},
            {"constantBothField", fField_constantBoth}
    };
    map<string, void (*)(const double * x, double * forcing_out)>::iterator it_f;
    it_f = lookup_f.find(conf.model_f);
    if (it_f != lookup_f.end()){
        if (verbose) cout << "Forcing: " << conf.model_f << endl;
        model_f = lookup_f[conf.model_f];
    } else {
        if (verbose) cout << "No forcing function chosen." << endl;
        model_f = ERROR_wrongAccess;
    }

    map<string, void (*)(const double * x, long labelx, const double * y, long labely,
                         const MeshType &mesh, double * kernel_val)> lookup_kernel = {
            {"constantTruncated", kernel_constantTruncated},
            {"notch", kernel_notch},
            {"labeledNotch", kernel_labeledNotch},
            {"constant", kernel_constant},
            {"constantLinf2D", kernel_constantLinf2D},
            {"labeled", kernel_labeled},
            {"integrable_sym", kernel_integrable_sym},
			{"integrable_unsym", kernel_integrable_unsym},
            {"constant3D", kernel_constant3D},
            {"constant1D", kernel_constant1D},
            {"antisymmetric1D", kernel_antisymmetric1D},
            {"parabola", kernel_parabola},
            {"linearPrototypeMicroelastic", kernel_linearPrototypeMicroelastic},
            {"linearPrototypeMicroelasticField", kernelField_linearPrototypeMicroelastic},
            {"linearPrototypeMicroelasticField3D", kernelField_linearPrototypeMicroelastic3D},
            {"constantField", kernelField_constant},
            {"fractional", kernel_fractional},
            {"fractional_shape", kernel_fractional_shape},
            {"labeledValve", kernel_labeledValve}
    };

    map<string, void (*)(const double * x, long labelx, const double * y, long labely,
                         const MeshType &mesh, double * kernel_val)>::iterator it_kernel;

    it_kernel = lookup_kernel.find(conf.model_kernel);
    if (it_kernel != lookup_kernel.end()){
        if (verbose) cout << "Kernel: " << conf.model_kernel << endl;
        model_kernel = lookup_kernel[conf.model_kernel];
    } else {
        if (verbose) cout << "No kernel chosen." << endl;
        model_kernel = ERROR_wrongAccess;
    }

    map<string, bool> lookup_singularKernels = {
            {"linearPrototypeMicroelastic", true},
            {"linearPrototypeMicroelasticField", true},
            {"fractional", true},
            {"fractional_shape", true}
    };
    conf.is_singularKernel = lookup_singularKernels[conf.model_kernel];


    // Lookup integration method  --------------------------------------------------------------------------------------
    if (conf.is_ShapeDerivative) {
        map<string, int (*)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                            const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                            double *termLocalPrime, double *termNonlocPrime, double *state_nodal_values, double *adjoint_nodal_values)> lookup_integrate {
                {"retriangulate_shape", integrate_retriangulate_shape},
                {"fractional_shape", integrate_fractional_shape},
                {"ERROR_wrongAccess_shape", ERROR_wrongAccess_shape}
        };
        map<string, int (*)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                            const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                            double *termLocalPrime, double *termNonlocPrime, double *state_nodal_values, double *adjoint_nodal_values)>::iterator it_integrator;
        it_integrator = lookup_integrate.find(conf.integration_method_remote);
        if (it_integrator != lookup_integrate.end()){
            if (verbose) cout << "Integration Method Remote: " << conf.integration_method_remote << endl;
            integrate_remote_shape = lookup_integrate[conf.integration_method_remote];
        } else {
            if (verbose) cout << "No integration routine for remote elements chosen." << endl;
            integrate_remote_shape = ERROR_wrongAccess_shape;
        }

        it_integrator = lookup_integrate.find(conf.integration_method_close);
        if (it_integrator != lookup_integrate.end()){
            if (verbose) cout << "Integration Method Close: " << conf.integration_method_close << endl;
            integrate_close_shape = lookup_integrate[conf.integration_method_close];
        } else {
            if (verbose) cout << "No integration routine for close elements chosen." << endl;
            integrate_close_shape = ERROR_wrongAccess_shape;
        }
    }
    else {
        map<string, int (*)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                            const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                            double *termLocalPrime, double *termNonlocPrime)> lookup_integrate {
                {"baryCenter", integrate_baryCenter},
                {"subSetBall", integrate_subSuperSetBalls},
                {"averageBall", integrate_subSuperSetBalls},
                {"baryCenterRT", integrate_baryCenterRT},
                {"retriangulate", integrate_retriangulate},
                {"retriangulate_unsymm", integrate_retriangulate_unysmm},
                {"retriangulate_unsymmLinfty", integrate_retriangulate_unysmm},
                {"retriangulateLinfty", integrate_retriangulate},
                {"exactBall", integrate_exact},
                {"noTruncation", integrate_fullyContained},
                {"fractional", integrate_fractional},
                {"weakSingular", integrate_weakSingular},
                {"ERROR_wrongAccess", ERROR_wrongAccess}
        };
        map<string, int (*)(const ElementType &aT, const ElementType &bT, const QuadratureType &quadRule, const MeshType &mesh,
                            const ConfigurationType &conf, bool is_firstbfslayer, double *termLocal, double *termNonloc,
                            double *termLocalPrime, double *termNonlocPrime)>::iterator it_integrator;

        it_integrator = lookup_integrate.find(conf.integration_method_remote);
        if (it_integrator != lookup_integrate.end()){
            if (verbose) cout << "Integration Method Remote: " << conf.integration_method_remote << endl;
            integrate_remote = lookup_integrate[conf.integration_method_remote];
        } else {
            if (verbose) cout << "No integration routine for remote elements chosen." << endl;
            integrate_remote = ERROR_wrongAccess;
        }

        it_integrator = lookup_integrate.find(conf.integration_method_close);
        if (it_integrator != lookup_integrate.end()){
            if (verbose) cout << "Integration Method Close: " << conf.integration_method_close << endl;
            integrate_close = lookup_integrate[conf.integration_method_close];
        } else {
            if (verbose) cout << "No integration routine for close elements chosen." << endl;
            integrate_close = ERROR_wrongAccess;
        }
    }

    if (conf.integration_method_remote == "retriangulateLinfty"
        || conf.integration_method_close == "retriangulateLinfty"
        || conf.integration_method_remote == "retriangulate_unsymmLinfty"
        || conf.integration_method_close == "retriangulate_unsymmLinfty") {
        method = method_retriangulateInfty;
    } else {
        // Becomes effective only if integrate_retriangulate is actually used. Has no meaning otherwise.
        method = method_retriangulate;
    }


}
// Compute A and f -----------------------------------------------------------------------------------------------------
void compute_f(     const ElementType & aT,
                    const QuadratureType &quadRule,
                    const MeshType & mesh,
                    double * termf){
    int i,a;
    double x[mesh.dim];
    double forcing_value[mesh.outdim];

    for (a=0; a<mesh.dVertex*mesh.outdim; a++){
        for (i=0; i<quadRule.nPx; i++){
            toPhys(aT.E, &(quadRule.Px[mesh.dim * i]),  mesh.dim,&x[0]);
            model_f(&x[0], forcing_value);
            termf[a] += quadRule.psix(a/mesh.outdim, i) * forcing_value[a%mesh.outdim] * aT.absDet * quadRule.dx[i];
        }
    }
}

void par_evaluateMass(double *vd, const double *ud, long *Elements,
                      const long *ElementLabels, const double *Verts, const long * VertexLabels, int K_Omega, int J, int nP,
                      double *P, const double *dx, const int dim, const int outdim, const int is_DiscontinuousGalerkin) {
    const int dVerts = dim+1;
    double tmp_psi[dVerts];
    auto *psi = (double *) malloc((dVerts)*nP*sizeof(double));

    for(int k=0; k<nP; k++){
        //model_basisFunction(const double * p, const MeshType & mesh, double *psi_vals){
       model_basisFunction(&P[dim*k], dim, &tmp_psi[0]);
       for (int kk=0; kk<dVerts; kk++) {
           psi[nP * kk + k] = tmp_psi[kk];
           //psi[nP * 1 + j] = tmp_psi[1];
           //psi[nP * 2 + j] = tmp_psi[2];
       }
    }
    // For some reason gcc does not allow to put const variables as shared. Const pointers however appear.
    // The constant function parameters can still be accsess inside the environemnt - might be different
    // in other compilers...
    #pragma omp parallel default(none) shared(ud, ElementLabels, Verts, VertexLabels, dx, nP, J, vd, psi, Elements, dVerts, dim, is_DiscontinuousGalerkin, outdim)
    {
    //private(aAdx, a, b, aTE, aTdet, j)
        double aTdet;
        long * aAdx;
        long aDGdx[dVerts]; // Index for discontinuous Galerkin

        double aTE[dim*(dVerts)];
        #pragma omp for
        for (int aTdx=0; aTdx < J; aTdx++){
            if (ElementLabels[aTdx] > 0) {
                // Get index of ansatz functions in matrix compute_A.-------------------

                // Discontinuous Galerkin
                if (is_DiscontinuousGalerkin) {
                    // Discontinuous Galerkin
                    //aDGdx[0] = (dVertex+1)*aTdx+1; aDGdx[1] = (dVertex+1)*aTdx+2; aDGdx[2] = (dVertex+1)*aTdx+3;
                    for (int j = 0; j < dVerts; j++) {
                        aDGdx[j] = dVerts * aTdx + j;
                    }
                    aAdx = aDGdx;
                } else {
                    // Continuous Galerkin
                    aAdx = &Elements[dVerts* aTdx];
                }

                // Prepare Triangle information aTE and aTdet ------------------
                // Copy coordinates of Triangel a to aTE.
                for (int jj=0; jj<dVerts; jj++){
                    for (int j = 0; j < dim; j++) {
                        aTE[dim * jj + j] = Verts[dim * Elements[dVerts * aTdx + jj] + j];
                    }
                    //aTE[2 * 0 + j] = Verts[2 * Elements[4 * aTdx + 1] + j];
                    //aTE[2 * 1 + j] = Verts[2 * Elements[4 * aTdx + 2] + j];
                    //aTE[2 * 2 + j] = Verts[2 * Elements[4 * aTdx + 3] + j];
                }
                // compute Determinant
                aTdet = absDet(&aTE[0], dim);

                for (int a = 0; a < dVerts; a++) {
                    if(is_DiscontinuousGalerkin || VertexLabels[ aAdx[a] ] > 0) {
                        for (int aOut = 0; aOut < outdim; aOut++) {
                            for (int b = 0; b < dVerts; b++) {
                                if(is_DiscontinuousGalerkin || VertexLabels[ aAdx[b] ] > 0) {
                                    for (int bOut = 0; bOut < outdim; bOut++) {
                                        for (int j = 0; j < nP; j++) {
                                            // Evaluation
#pragma omp atomic update
                                            vd[outdim*aAdx[a] + aOut] +=
                                                    psi[nP * a + j] * psi[nP * b + j] * aTdet * dx[j] * ud[outdim*aAdx[b] + bOut];
                                        }
                                    }
                                }
                            }
                        }
                   }
                }
            }
        }
    } // Pragma Omp Parallel
    free(psi);
}
//TODO This algorithm is very, very similar to the main algorithm...!
// Can BFS Algorithms should somehow be separated..?
void estimateNNZperRow(const MeshType & mesh, const ConfigurationType & conf){
    const int sampleSize = 3;
    int indexList[sampleSize];
    int nE = mesh.nE;
    bool isDG = mesh.is_DiscontinuousGalerkin;

    ElementType aT, bT;
    aT.matE = arma::vec(mesh.dim*(mesh.dim+1));
    bT.matE = arma::vec(mesh.dim*(mesh.dim+1));
    
    arma::Col<int> rowCount(mesh.K, arma::fill::zeros);
    arma::Col<int> bvertexVisited(mesh.nV, arma::fill::zeros);
    arma::Col<int> avertexVisited(mesh.nV, arma::fill::zeros);

    // Traverse Triangles by Neighbours
    for (int & k : indexList){
        k = rand() % nE;
    }
    // Breadth First Search --------------------------------------
    // Every Thread owns its copy!
    arma::Col<int> visited(nE, arma::fill::zeros);

    // Queue for Breadth first search
    queue<int> queue;
    // List of visited neighbours
    //const long *NTdx;
    //for (int aTdx=0; aTdx<mesh.nE; aTdx++)
    for (int & aTdx : indexList){
        initializeElement(aTdx, mesh, aT);
        // Intialize search queue with current outer triangle
        queue.push(aTdx);
        // Initialize vector of visited triangles with 0
        visited.zeros();

        while (!queue.empty()) {
            // Get and delete the next Triangle indexList of the queue. The first one will be the triangle aTdx itself.
            int sTdx = queue.front();
            queue.pop();
            // Get all the neighbours of sTdx.
            //NTdx = &mesh.Neighbours(0, sTdx);
            idx_t startNdx = mesh.xadj[sTdx];
            //cout << "startNdx " << startNdx << endl;
            idx_t endNdx = mesh.xadj[sTdx+1];
            //cout << "endNdx "<< endNdx << endl;
            idx_t nTdx = -1;

            // Run through the list of neighbours.
            while (startNdx + nTdx < endNdx) {
            //for (int j = 0; j < mesh.nNeighbours; j++) {

                // The next valid neighbour is our candidate for the inner Triangle b.
                // int bTdx = NTdx[j];
                int bTdx;
                if (nTdx == -1) {bTdx = aTdx;} // Make sure, that the first 'neighbour' is the element itself
                //TODO Narrowing conversion!
                else {bTdx = mesh.adjncy[startNdx + nTdx];}
                // Check how many neighbours sTdx has.
                // In order to be able to store the list as contiguous array we fill
                // up the empty spots with the number nE
                // i.e. the total number of Triangles (which cannot be an indexList).
                //if (bTdx < mesh.nE) {
                    // Check whether bTdx is already visited.
                    if (!visited(bTdx)){
                        initializeElement(bTdx, mesh, bT);

                        double abary[mesh.dim], bbary[mesh.dim];
                        baryCenter(mesh.dim, aT.E, abary);
                        baryCenter(mesh.dim, bT.E, bbary);
                        double dist = sqrt(vec_sqL2dist(abary, bbary, mesh.dim));
                        if (dist < mesh.delta + mesh.maxDiameter){
                            queue.push(bTdx);

                            if (isDG) {
                                rowCount(3*aTdx) += 3;
                                rowCount(3*aTdx+1) += 3;
                                rowCount(3*aTdx+2) += 3;
                            } else {
                                for (int i = 0; i < mesh.dVertex; i++) {
                                    //TODO Narrowing conversion!
                                    int bVdx = mesh.Triangles(i, bTdx);
                                    if (!bvertexVisited(bVdx)) {
                                        for (int k = 0; k < mesh.dVertex; k++) {
                                            //TODO Narrowing conversion!
                                            int aVdx = mesh.Triangles(k, aTdx);
                                            if (!avertexVisited(aVdx)) {
                                                //Vdx = mesh.ptrTriangles[(mesh.dVertex+1)*Tdx + k+1];
                                                rowCount(aVdx) += 1;
                                            }
                                        } // End for (int k = 0; k < mesh.dVertex; k++)
                                        bvertexVisited(bVdx) = 1;
                                    } // End if (!bvertexVisited(bVdx))
                                } // End for (int i = 0; i < mesh.dVertex; i++)
                            } // End if (isDG)

                        } // End  if (dist < mesh.delta + mesh.maxDiameter)
                    } // End visited[bTdx] = 1;
                    visited[bTdx] = 1;
                //} // if (bTdx < mesh.nE)
                nTdx++;
            } // End for (int j = 0; j < mesh.nNeighbours; j++)
        } // End while (!queue.empty())
        for (int k = 0; k < mesh.dVertex; k++) {
            int aVdx = mesh.Triangles(k, aTdx);
            avertexVisited(aVdx) = 1;
        }
    } // End for (int aTdx=0; aTdx<mesh.nE; aTdx++)

    printf("C++ Estimated Row NNZ %i\n", arma::max(rowCount));
} // End estimateNNZperRow

// Assembly algorithm with BFS -----------------------------------------------------------------------------------------
void par_assemble(const string compute, const string path_spAd, const string path_fd, const int K_Omega, const int K,
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
                  int is_ShapeDerivative, const double *state, const double *adjoint
                  ) {

    MeshType mesh = {K_Omega, K, ptrTriangles, ptrLabelTriangles, ptrVerts, ptrLabelVerts, nE, nE_Omega,
                     nV, nV_Omega, sqrt(sqdelta), sqdelta,
                     //ptrNeighbours, nNeighbours,
                     is_DiscontinuousGalerkin,
                     is_NeumannBoundary, dim, outdim, dim+1,
                     ptrZetaIndicator_indices, ptrZetaIndicator_indptr, nZeta, maxDiameter, fractional_s};

    ConfigurationType conf = {path_spAd, path_fd, str_model_kernel, str_model_f,
                              str_integration_method_remote, str_integration_method_close,
                              static_cast<bool>(is_PlacePointOnCap),
                              false,
                              static_cast<bool>(is_fullConnectedComponentSearch),
                              verbose,
                              is_ShapeDerivative};

    if (verbose) printf("Constructing adjacency graph...");
    int dVertex = mesh.dim+1;
    idx_t nE_metis=mesh.nE;
    idx_t nV_metis=mesh.nV;
    idx_t ncommon=1;
    idx_t eind[mesh.nE * dVertex];
    idx_t eptr[mesh.nE + 1];
    idx_t numflag=0;
    idx_t *xadj;
    idx_t *adjncy;

    for (int k=0; k<mesh.nE; k++){
        //TODO Be careful with the casts etc..
        eptr[k] =  (idx_t)(k*mesh.dVertex);
        for (int l=0; l<mesh.dVertex; l++){
            eind[mesh.dVertex*k + l] = (idx_t)(mesh.Triangles[mesh.dVertex*k + l]);
        }
    }
    eptr[mesh.nE] =  (idx_t)(mesh.nE*mesh.dVertex);

    // Compute Adjacency Graph with METIS
    //TODO Write a wrapper for metis somehow...
    METIS_MeshToDual(&nE_metis, &nV_metis, eptr, eind, &ncommon, &numflag, &xadj, &adjncy);
    if (verbose) printf("Done. \n");

    mesh.xadj = xadj;
    mesh.adjncy = adjncy;
    mesh.eptr = eptr;
    mesh.eind = eind;

    QuadratureType quadRule = {Px, Py, dx, dy, nPx, nPy, dim, Pg, dg, degree};
    initializeQuadrule(quadRule, mesh);

    if (compute=="system") {
        //map<unsigned long, double> Ad;
        chk_Mesh(mesh, verbose);
        chk_Conf(mesh, conf, quadRule);

        //estimateNNZperRow(mesh, conf);
        par_system(mesh, quadRule, conf, state, adjoint);
        /*
        if (verbose) cout << "K_Omega " << mesh.K_Omega << endl;
        if (verbose) cout << "K " << mesh.K << endl;

        int nnz_total = static_cast<int>(Ad.size());
        arma::vec values_all(nnz_total);
        arma::umat indices_all(2, nnz_total);
        if (verbose) cout << "Total NNZ " << nnz_total << endl;

        int k = 0;
        for (auto &it : Ad) {
            unsigned long adx = it.first;
            double value = it.second;
            values_all(k) = value;
            // column major format of transposed matrix Ad
            indices_all(0, k) = adx % mesh.K;
            indices_all(1, k) = adx / mesh.K;
            //printf("Index a %llu, b %llu, k %i\n", indices_all(0, k), indices_all(1, k), k);
            k++;
        }
        arma::sp_mat sp_Ad(true, indices_all, values_all, mesh.K, mesh.K);

        sp_Ad.save(conf.path_spAd);
        */
    }

    if (compute=="forcing") {
        par_forcing(mesh, quadRule, conf);
    }

    METIS_Free(xadj);
    METIS_Free(adjncy);

}

void par_system(MeshType &mesh, QuadratureType &quadRule, ConfigurationType &conf,
                const double *state, const double *adjoint) {

    const int verbose = conf.verbose;
    if (verbose) printf("Function: par_system\n");
    if (verbose) printf("Ansatz Space: %s\n", mesh.is_DiscontinuousGalerkin ? "DG" : "CG");
    if (verbose) printf("Mesh dimension: %i\n", mesh.dim);
    if (verbose) printf("Output dimension: %i\n", mesh.outdim);
    if (verbose) printf("Recieved Zeta for DD: %s\n", (mesh.nZeta) ? "true" : "false");
    lookup_configuration(conf, verbose);
    if (verbose) printf("Quadrule outer: %i\n", quadRule.nPx);
    if (verbose) printf("Quadrule inner: %i\n", quadRule.nPy);
    if (verbose) printf("Full Graph Search: %i\n", conf.is_fullConnectedComponentSearch);
    arma::vec values_all;
    arma::umat indices_all(2,0);
    int nnz_total = 0;
    idx_t nn_idxt = mesh.nV;
    idx_t ne_idxt = mesh.nE;
    //cout << "nV " << mesh.nV << " nE " << mesh.nE << endl;
    //cout << "nV " << nn_idxt << " nE " << ne_idxt << endl;
    idx_t nparts = 1;
    idx_t ncommon = 2;
    idx_t objval=0;
    idx_t epart[ne_idxt];
    idx_t npart[nn_idxt];
    //printf("npart %ld\n", nparts);
    #pragma omp parallel
    {

        #pragma omp single
        {
            nparts = omp_get_num_threads();
        }
    }
    //printf("npart %ld -> METIS \n", nparts);
    if (nparts > 1) {
        METIS_PartMeshDual(&ne_idxt, &nn_idxt, mesh.eptr, mesh.eind,
                           nullptr, nullptr, &ncommon,
                           &nparts, nullptr, nullptr,
                           &objval, epart, npart);
    }
    //cout << epart << endl;
    //cout << npart << endl;
    //TODO psix should be a OnePointData (like OnePointFunction but for precomputed data)
    // Basically OnePointData is a performance tweak for ReferencePoints (in opposition ti PhysicalPoints)
    //can they also depend on points only ?
    for(int h=0; h<quadRule.nPx; h++){
        // This works due to Column Major ordering of Armadillo Matricies!
        model_basisFunction(& quadRule.Px[mesh.dim*h], mesh.dim, & quadRule.psix[mesh.dVertex * h]);
    }
    for(int h=0; h<quadRule.nPy; h++){
        // This works due to Column Major ordering of Armadillo Matricies!
        model_basisFunction(& quadRule.Py[mesh.dim*h], mesh.dim,& quadRule.psiy[mesh.dVertex * h]);
    }
    chk_BasisFunction(quadRule);
    auto zeros = arma::fill::zeros;
    #pragma omp parallel default(none) shared(mesh, quadRule, conf, values_all, indices_all, nnz_total, epart, nparts, verbose, zeros, state, adjoint)
    {
    map<unsigned long, double> Ad;
    unsigned long Adx;
    idx_t tid = 0;
    if (nparts > 1){
        tid = omp_get_thread_num();
    }
    if (verbose && !tid) printf("Thread %ld of %ld  (num partitions)\n", tid, nparts);
    // Breadth First Search --------------------------------------
    arma::Col<int> visited(mesh.nE, zeros);

    // Queue for Breadth first search
    //queue<idx_t> queue;
    queue<int> queue;
    // List of visited triangles
    //const long *NTdx;

    // Vector containing the coordinates of the vertices of a Triangle
    ElementType aT, bT;
    aT.matE = arma::vec(mesh.dim*(mesh.dim+1));
    bT.matE = arma::vec(mesh.dim*(mesh.dim+1));

    // Integration information ------------------------------------
    // (Pointer to) Vector of indices of Basisfuntions (Adx) for triangle a and b
    const long * aAdx;
    const long * bAdx;
    long aDGdx[mesh.dVertex]; // Index for discontinuous Galerkin
    long bDGdx[mesh.dVertex];

    // Buffers for integration solutions
    double termLocal[mesh.dVertex*mesh.dVertex*mesh.outdim*mesh.outdim];
    double termNonloc[mesh.dVertex*mesh.dVertex*mesh.outdim*mesh.outdim];
    double termLocalPrime[mesh.dVertex*mesh.dVertex*mesh.outdim*mesh.outdim];
    double termNonlocPrime[mesh.dVertex*mesh.dVertex*mesh.outdim*mesh.outdim];
    //#pragma omp for
    for (int aTdx=0; aTdx<mesh.nE; aTdx++) {
        //printf("epart[aTdx] %i", epart[aTdx]);

        const bool is_mypart = (nparts == 1) || (tid == epart[aTdx]);
        if (is_mypart && (mesh.LabelTriangles[aTdx])) {

            //[DEBUG]
            /*
            if (aTdx==7){
                cout << endl << "Total Local Term" << endl ;
                printf ("[%17.16e, %17.16e, %17.16e] \n", DEBUG_termTotalLocal[0], DEBUG_termTotalLocal[1], DEBUG_termTotalLocal[2]);
                printf ("[%17.16e, %17.16e, %17.16e] \n", DEBUG_termTotalLocal[3], DEBUG_termTotalLocal[4], DEBUG_termTotalLocal[5]);
                printf ("[%17.16e, %17.16e, %17.16e] \n", DEBUG_termTotalLocal[6], DEBUG_termTotalLocal[7], DEBUG_termTotalLocal[8]);

                cout << endl << "Total Nonlocal Term" << endl ;
                printf ("[%17.16e, %17.16e, %17.16e] \n", DEBUG_termTotalNonloc[0], DEBUG_termTotalNonloc[1], DEBUG_termTotalNonloc[2]);
                printf ("[%17.16e, %17.16e, %17.16e] \n", DEBUG_termTotalNonloc[3], DEBUG_termTotalNonloc[4], DEBUG_termTotalNonloc[5]);
                printf ("[%17.16e, %17.16e, %17.16e] \n", DEBUG_termTotalNonloc[6], DEBUG_termTotalNonloc[7], DEBUG_termTotalNonloc[8]);

                abort();
            }
            */
            //[DEBUG]

            // Get index of Ansatz functions in matrix compute_A.-------------------
            if (mesh.is_DiscontinuousGalerkin) {
                // Discontinuous Galerkin
                for (int j = 0; j < mesh.dVertex; j++) {
                    aDGdx[j] = mesh.dVertex * aTdx + j;
                }
                aAdx = aDGdx;
            } else {
                // Continuous Galerkin
                aAdx = &mesh.Triangles(0, aTdx);
            }
            // Prepare Triangle information aTE and aTdet ------------------
            initializeElement(aTdx, mesh, aT);

            double state_nodal_values[6];
            double adjoint_nodal_values[6];
            if (conf.is_ShapeDerivative) {
                for (int i=0; i < mesh.dVertex; i++) {
                    state_nodal_values[i] = state[aAdx[i]];
                    adjoint_nodal_values[i] = adjoint[aAdx[i]];
                }
            }

            // BFS -------------------------------------------------------------
            // Intialize search queue with current outer triangle
            queue.push(aTdx);
            // Initialize vector of visited triangles with 0
            visited.zeros();

            // Tells that we are in the first layer of the BFS
            bool is_firstbfslayer = conf.is_singularKernel;
            // Check whether BFS is completed.
            while (!queue.empty()) {
                // Get and delete the next Triangle index of the queue. The first one will be the triangle aTdx itself.
                int sTdx = queue.front();
                //cout << "sTdx " << sTdx << endl;
                queue.pop();
                // Get all the neighbours of sTdx.
                //NTdx = &mesh.Neighbours(0, sTdx);
                idx_t startNdx = mesh.xadj[sTdx];
                //cout << "startNdx " << startNdx << endl;
                idx_t endNdx = mesh.xadj[sTdx+1];
                //cout << "endNdx "<< endNdx << endl;
                idx_t nTdx = -1;

                // Run through the list of neighbours.
                while (startNdx + nTdx != endNdx) {
                //for (int j = 0; j < mesh.nNeighbours; j++) {
                    // The next valid neighbour is our candidate for the inner Triangle b.
                    //int bTdx = NTdx[j];
                    int bTdx;
                    if (nTdx == -1) {bTdx = aTdx;} // Make sure, that the first 'neighbour' is the element itself
                    else {bTdx = mesh.adjncy[startNdx + nTdx];}

                    // Check how many neighbours sTdx has.
                    // In order to be able to store the list as contiguous array we fill
                    // up the empty spots with the number nE
                    // i.e. the total number of Triangles (which cannot be an index).
                    //if (bTdx < mesh.nE) {
                        // Check whether bTdx is already visited.
                        if (!visited[bTdx]) {
                            // Check whether bTdx is part of the discretization
                            // otherwise it is just appended to the queue
                            if (mesh.LabelTriangles[bTdx]) {
                                // Prepare Triangle information bTE and bTdet ------------------
                                initializeElement(bTdx, mesh, bT);
                                // Retriangulation and integration -----------------------------
                                if (mesh.is_DiscontinuousGalerkin) {
                                    // Discontinuous Galerkin
                                    for (int jj = 0; jj < mesh.dVertex; jj++) {
                                        bDGdx[jj] = mesh.dVertex * bTdx + jj;
                                    }
                                    bAdx = bDGdx;
                                } else {
                                    // Get (pointer to) index of basis function (in Continuous Galerkin)
                                    bAdx = &mesh.Triangles(0, bTdx);
                                }

                                // Domain decomposition. If Zeta is empty, the weight is set to 1.
                                double weight = 1.;
                                //printf("aTdx %i, bTdx %i \n", aTdx, bTdx);
                                if (mesh.nZeta) {
                                    double zeta = evaluateZeta(mesh.ptrZetaIndicator_indices,
                                                 mesh.ptrZetaIndicator_indptr,
                                                 mesh.nZeta,
                                                 aTdx, bTdx);
                                    weight = 1./zeta;
                                }

                                // Assembly of matrix ---------------------------------------
                                doubleVec_tozero(termLocal, mesh.dVertex * mesh.dVertex * mesh.outdim *
                                                            mesh.outdim); // Initialize Buffer
                                doubleVec_tozero(termNonloc, mesh.dVertex * mesh.dVertex * mesh.outdim *
                                                             mesh.outdim); // Initialize Buffer
                                doubleVec_tozero(termLocalPrime, mesh.dVertex * mesh.dVertex * mesh.outdim *
                                                                 mesh.outdim); // Initialize Buffer
                                doubleVec_tozero(termNonlocPrime, mesh.dVertex * mesh.dVertex * mesh.outdim *
                                                                  mesh.outdim); // Initialize Buffer

                                // Compute integrals and write to buffer
                                int doesInteract = 0;
                                if (conf.is_ShapeDerivative) {
                                    for (int i=0; i < mesh.dVertex; i++) {
                                        state_nodal_values[i + mesh.dVertex] = state[bAdx[i]];
                                        adjoint_nodal_values[i + mesh.dVertex] = adjoint[bAdx[i]];
                                    }
                                    doesInteract = integrate_shape(aT, bT, quadRule, mesh, conf, is_firstbfslayer,
                                                             termLocal, termNonloc, termLocalPrime, termNonlocPrime,
                                                             &state_nodal_values[0], &adjoint_nodal_values[0]);
                                }
                                else {
                                    doesInteract = integrate(aT, bT, quadRule, mesh, conf, is_firstbfslayer,
                                                             termLocal, termNonloc, termLocalPrime, termNonlocPrime);
                                }
                                if (doesInteract) {

                                    // If bT interacts it will be a candidate for our BFS, so it is added to the queue

                                    //[DEBUG]
                                    /*
                                    if (is_firstbfslayer){//aTdx == 0 && bTdx == 0){

                                    printf("aTdx %i\ndet %17.16e, label %li \n", aTdx, aT.absDet, aT.label);
                                    printf ("aTE\n[%17.16e, %17.16e]\n[%17.16e, %17.16e]\n[%17.16e, %17.16e]\n", aT.E[0],aT.E[1],aT.E[2],aT.E[3],aT.E[4],aT.E[5]);
                                    printf("bTdx %i\ndet %17.16e, label %li \n", bTdx, bT.absDet, bT.label);
                                    printf ("bTE\n[%17.16e, %17.16e]\n[%17.16e, %17.16e]\n[%17.16e, %17.16e]\n", bT.E[0],bT.E[1],bT.E[2],bT.E[3],bT.E[4],bT.E[5]);

                                    cout << endl << "Local Term" << endl ;
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termLocal[0], termLocal[1], termLocal[2]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termLocal[3], termLocal[4], termLocal[5]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termLocal[6], termLocal[7], termLocal[8]);

                                    cout << endl << "Nonlocal Term" << endl ;
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termNonloc[0], termNonloc[1], termNonloc[2]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termNonloc[3], termNonloc[4], termNonloc[5]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termNonloc[6], termNonloc[7], termNonloc[8]);

                                    cout << endl << "Local Term Prime" << endl ;
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termLocalPrime[0], termLocalPrime[1], termLocalPrime[2]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termLocalPrime[3], termLocalPrime[4], termLocalPrime[5]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termLocalPrime[6], termLocalPrime[7], termLocalPrime[8]);

                                    cout << endl << "Nonlocal Term Prime" << endl ;
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termNonlocPrime[0], termNonlocPrime[1], termNonlocPrime[2]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termNonlocPrime[3], termNonlocPrime[4], termNonlocPrime[5]);
                                    printf ("[%17.16e, %17.16e, %17.16e] \n", termNonlocPrime[6], termNonlocPrime[7], termNonlocPrime[8]);
                                    //if (aTdx == 10){
                                    //    abort();
                                    //}

                                    }
                                    if ((aTdx == 1) && !is_firstbfslayer) abort();
                                    */
                                    //[End DEBUG]

                                    queue.push(bTdx);
                                    // We only check whether the integral
                                    // (termLocal, termNonloc) are 0, in which case we dont add bTdx to the queue.
                                    // However, this works only if we can guarantee that interacting triangles do actually
                                    // also contribute a non-zero entry, i.e. the Kernel as to be > 0 everywhere on its support for example.
                                    // The effect (in speedup) of this more precise criterea depends on delta and meshsize.

                                    if (conf.is_ShapeDerivative) {
                                        double factor = 1.0;
                                        if (aTdx == bTdx) {
                                            factor = 0.5;
                                        }
                                        double  div[6];
                                        get_div(aT.E, &div[0]);
                                        for (int a = 0; a < mesh.dVertex; a++) {
                                            if (mesh.LabelVerts[aAdx[a]] > 0) {
                                                Ad[aAdx[a]] += 2 * termLocal[0] * div[2*a] + factor * termNonloc[a];
                                                Ad[aAdx[a] + mesh.K] += 2 * termLocal[0] * div[2*a + 1]
                                                        + factor * termNonlocPrime[a];
                                            }
                                        }

                                    }
                                    else{
                                        for (int a = 0; a < mesh.dVertex * mesh.outdim; a++) {
                                            for (int b = 0; b < mesh.dVertex * mesh.outdim; b++) {
                                                //printf("Local: a %i, b %i, ker (%i, %i) \nAdx %lu \n", a/mesh.outdim, b/mesh.outdim, a%mesh.outdim, b%mesh.outdim, Adx);
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

                                                if (mesh.is_DiscontinuousGalerkin ||
                                                    (mesh.LabelVerts[aAdx[a / mesh.outdim]] > 0)) {
                                                    Adx = (mesh.outdim * aAdx[a / mesh.outdim] + a % mesh.outdim) * mesh.K
                                                            + mesh.outdim * aAdx[b / mesh.outdim] + b % mesh.outdim;
//#pragma omp critical
                                                    {
                                                        Ad[Adx] += termLocal[mesh.dVertex * mesh.outdim * a + b] * weight;
                                                    }

                                                    Adx = (mesh.outdim * aAdx[a / mesh.outdim] + a % mesh.outdim) * mesh.K
                                                          + mesh.outdim * bAdx[b / mesh.outdim] + b % mesh.outdim;
//#pragma omp critical
                                                    {
                                                        Ad[Adx] += -termNonloc[mesh.dVertex * mesh.outdim * a + b] * weight;
                                                    }
                                                }
                                                if (mesh.is_DiscontinuousGalerkin ||
                                                    (mesh.LabelVerts[bAdx[b / mesh.outdim]] > 0)) {
                                                    Adx = (mesh.outdim * bAdx[b / mesh.outdim] + b % mesh.outdim) * mesh.K
                                                          + mesh.outdim * bAdx[a / mesh.outdim] + a % mesh.outdim;
//#pragma omp critical
                                                    {
                                                        Ad[Adx] +=
                                                                termLocalPrime[mesh.dVertex * mesh.outdim * a + b] * weight;
                                                    }

                                                    Adx = (mesh.outdim * bAdx[b / mesh.outdim] + b % mesh.outdim) * mesh.K
                                                          + mesh.outdim * aAdx[a / mesh.outdim] + a % mesh.outdim;
//#pragma omp critical
                                                    {
                                                        Ad[Adx] += -termNonlocPrime[mesh.dVertex * mesh.outdim * a + b] *
                                                                   weight;
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // Copy buffer into matrix. Again solutions which lie on the boundary are ignored (in Continuous Galerkin)

                                }
                            }// End if(mesh.LabelTriangles[bTdx])
                            else {
                                queue.push(bTdx);
                            }
                        }// End if BFS (visited[bTdx] == 0)
                        // Mark bTdx as visited
                        visited[bTdx] = 1;
                    //}// End if BFS (bTdx < mesh.nE)
                    nTdx++;
                }//End for loop BFS (j = 0; j < mesh.nNeighbours; j++)
                is_firstbfslayer = false;
            }//End while loop BFS (!queue.empty())
       }// End if LabelTriangles != 0
    }// End parallel for


    /*
    arma::mat dualGraph(15, mesh.nE);
    dualGraph.fill(0);

    for (idx_t i = 0; i < mesh.nE; i++){
        idx_t s = xadj[i];
        idx_t j = 0;
        printf("i %i \n", i);
        printf("bound: %i \n", xadj[i+1]);
        while(s + j < xadj[i+1]){
            printf("j %i, ", j);
            dualGraph(j, i) = adjncy[s + j];
            j++;
        }
    }
    dualGraph.save("data/result.dual_par_system", arma::arma_binary);
    */
    int nnz_start = 0;
    #pragma omp critical
        {
            int nnz_current = static_cast<int>(Ad.size());
            //printf("NNZ of Thread %i is %i\n", omp_get_thread_num(), nnz_current);
            //int estimatedNNZ = chunkSize * pow(2*ceil(mesh.delta / mesh.maxDiameter + 1), mesh.dim);
            //if (conf.verbose) printf("Estimated NNZ is %i\n", estimatedNNZ);
            nnz_start = nnz_total;
            nnz_total += nnz_current;
            //cout << "Thread "<< omp_get_thread_num() << ", start  " << nnz_start << endl;
        }
    #pragma omp barrier
    #pragma omp single
        {
            //if (conf.verbose) cout << "NNZ (total of all threads) "<< nnz_total << endl;
            //if (conf.verbose) cout << "Try to allocate ";
            values_all.set_size(nnz_total);
            indices_all.reshape(2, nnz_total);
            //if (conf.verbose) cout << " - Done - " << endl;
        }
    #pragma omp barrier
    #pragma omp critical
        {
            //cout << "Thread "<< omp_get_thread_num() << ", start  " << nnz_start << endl;
            int k = 0;
            for (auto & it : Ad) {
                unsigned long adx = it.first;
                double value = it.second;
                values_all(nnz_start + k) = value;
                // column major format of transposed matrix Ad
                if (conf.is_ShapeDerivative) {
                    indices_all(0, nnz_start + k) = 0;
                    indices_all(1, nnz_start + k) = adx;
                } else {
                    indices_all(0, nnz_start + k) = adx % mesh.K;
                    indices_all(1, nnz_start + k) = adx / mesh.K;
                }
                //printf("Index a %llu, b %llu\n", indices_all(0, nnz_start + k), indices_all(1, nnz_start + k));
                k++;
            }
        }

    }// End pragma omp parallel
    //indices_all.save("indices_all");
    //values_all.save("values_all");
    if (conf.verbose) cout << "K_Omega " << mesh.K_Omega << endl;
    if (conf.verbose) cout << "K " << mesh.K << endl;
    if (conf.verbose) cout << "NNZ (total of all threads) " << nnz_total << endl;
    //cout << arma::max(indices_all.row(1)) << endl;
    //TODO ARMADILLO_DEP
    int num_rows = mesh.K;
    int num_cols = mesh.K;
    if (conf.is_ShapeDerivative) {
        num_rows = 1;
        num_cols = 2*mesh.K;
    }
    arma::sp_mat sp_Ad(true, indices_all, values_all, num_rows, num_cols);
    sp_Ad.save(conf.path_spAd);
    //cout << "Data saved." << endl;


}// End function par_system

void par_forcing(MeshType &mesh, QuadratureType &quadRule, ConfigurationType &conf) {
    arma::vec fd(mesh.K, arma::fill::zeros);
    const int verbose = conf.verbose;
    if (verbose) printf("Function: par_forcing\n");
    //printf("Mesh dimension: %i\n", mesh.dim);
    //printf("Recieved Zeta for DD: %s\n", (mesh.ptrZeta) ? "true" : "false");
    lookup_configuration(conf);
    //printf("Quadrule outer: %i\n", quadRule.nPx);
    //printf("Quadrule inner: %i\n", quadRule.nPy);

    for (int h = 0; h < quadRule.nPx; h++) {
        // This works due to Column Major ordering of Armadillo Matricies!
        model_basisFunction(&quadRule.Px[mesh.dim * h], mesh.dim, &quadRule.psix[mesh.dVertex * h]);
    }

    #pragma omp parallel
    {
        // General Loop Indices ---------------------------------------
        // Vector containing the coordinates of the vertices of a Triangle
        ElementType aT;
        aT.matE = arma::vec(mesh.dim * (mesh.dim + 1));
        // (Pointer to) Vector of indices of Basisfuntions (Adx) for triangle a and b
        const long *aAdx;
        long aDGdx[mesh.dVertex]; // Index for discontinuous Galerkin
        // Buffers for integration solutions
        double termf[mesh.dVertex*mesh.outdim];
        #pragma omp for
        for (int aTdx = 0; aTdx < mesh.nE; aTdx++) {
            if (mesh.LabelTriangles[aTdx] > 0) {
                // Get index of ansatz functions in matrix compute_A.-------------------
                if (mesh.is_DiscontinuousGalerkin) {
                    // Discontinuous Galerkin
                    for (int j = 0; j < mesh.dVertex; j++) {
                        aDGdx[j] = mesh.dVertex * aTdx + j;
                    }
                    aAdx = aDGdx;
                } else {
                    // Continuous Galerkin
                    aAdx = &mesh.Triangles(0, aTdx);
                }
                // Prepare Triangle information aTE and aTdet ------------------
                initializeElement(aTdx, mesh, aT);
                // Assembly of right side ---------------------------------------
                // We unnecessarily integrate over vertices which might lie on the boundary of Omega for convenience here.
                doubleVec_tozero(termf, mesh.dVertex*mesh.outdim); // Initialize Buffer

                compute_f(aT, quadRule, mesh, termf); // Integrate and fill buffer
                // Add content of buffer to the right side.

                // Domain decomposition. If Zeta is empty, the weight is set to 1.
                double weight = 1.;
                //printf("aTdx %i, bTdx %i \n", aTdx, bTdx);
                if (mesh.nZeta) {
                    double zeta = evaluateZeta(mesh.ptrZetaIndicator_indices,
                                               mesh.ptrZetaIndicator_indptr,
                                               mesh.nZeta,
                                               aTdx, aTdx);
                    weight = 1./zeta;
                }

                for (int a = 0; a < mesh.dVertex*mesh.outdim; a++) {
                    if (mesh.is_DiscontinuousGalerkin || (mesh.LabelVerts[aAdx[a/mesh.outdim]] > 0)) {
                        #pragma omp atomic update
                        fd[mesh.outdim*aAdx[a/mesh.outdim] + a%mesh.outdim] += termf[a]*weight;
                    }
                }// end for rhs

            }// end outer if (mesh.LabelTriangles[aTdx] > 0)
        }// end outer for loop (aTdx=0; aTdx<mesh.nE; aTdx++)
    }// end pragma omp parallel
    fd.save(conf.path_fd);
}// end par_righthandside

// [DEBUG] _____________________________________________________________________________________________________________
/*
double compute_area(double * aTE, double aTdet, long labela, double * bTE, double bTdet, long labelb, double * P, int nP, double * dx, double sqdelta){
    double areaTerm=0.0;
    int rTdx, Rdx;
    double * x;
    double physical_quad[2];
    double reTriangle_list[9*3*2];
    const MeshType mesh = {K_Omega, K, ptrTriangles, ptrVerts, nE, nE_Omega,
                                        L, L_Omega, sqdelta, ptrNeighbours, is_DiscontinuousGalerkin,
                                        is_NeumannBoundary, dim, dim+1};;
    x = &P[0];
    toPhys(aTE, x, physical_quad);
    Rdx = retriangulate(physical_quad, bTE, mesh, reTriangle_list, true);
    for (rTdx=0; rTdx < Rdx; rTdx++){
        areaTerm += absDet(&reTriangle_list[2*3*rTdx])/2;
    }
    return areaTerm;
}
*/

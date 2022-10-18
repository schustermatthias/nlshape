#-*- coding:utf-8 -*-
#distutils: include_dirs = include

#cython: language_level=3
#cython: boundscheck=False, wraparound=False, cdivision=True

# Assembly routine
from libcpp.string cimport string
cimport Cassemble

import numpy as np
import time
from libc.math cimport pow
import scipy.sparse as sparse
cimport numpy as c_np
import datetime

def timestamp():
    """
    Returns current timestamp as string.

    :return: string, format %m%d_%H-%M-%S
    """
    # Link to strftime Doc
    # http://strftime.org/
    return datetime.datetime.now().strftime("%m%d_%H-%M-%S")

def read_arma_mat(path, is_verbose=False):
    """
    Read armadillo vector format from file.

    :param path: string, Path to file.
    :param is_verbose: bool, Verbose mode.
    :return: np.ndarray, Vector of double.
    """
    import numpy as np

    sizeof_double = 8

    f = open(path, "rb")
    # Read Armadillo header
    arma_header = f.readline()
    if arma_header != b'ARMA_MAT_BIN_FN008\n':
        raise ValueError("in read_arma_mat(), input file is of wrong format.")
    # Get shape of sparse matrix
    arma_shape = f.readline()
    n_rows, n_cols = tuple([int(x) for x in arma_shape.decode("utf-8").split()])
    if is_verbose: print("Shape (", n_rows, ", ", n_cols, ")", sep="")
    # Raw binary of sparse Matrix in csc-format
    b_data = f.read()
    f.close()

    b_values = b_data[:sizeof_double * n_rows * n_cols]
    values = np.array(np.frombuffer(b_values)).reshape((n_rows, n_cols), order="F")
    if is_verbose: print("Values ", values)
    if is_verbose: print(values)

    return values

def read_arma_spMat(path, is_verbose=False):
    """
    Read sparse Matrix in armadillo spMat format from file.

    :param path: string, Path to file.
    :param is_verbose: bool, Verbose mode.
    :return: scipy.csr_matrix, Matrix of double.
    """

    import scipy.sparse as sp
    import numpy as np

    sizeof_double = 8

    f = open(path, "rb")
    # Read Armadillo header
    arma_header = f.readline()
    if arma_header != b'ARMA_SPM_BIN_FN008\n':
        raise ValueError("in read_arma_spMat(), input file is of wrong format.")
    # Get shape of sparse matrix
    arma_shape = f.readline()
    n_rows, n_cols, n_nonzero = tuple([int(x) for x in arma_shape.decode("utf-8").split()])
    if is_verbose: print("Shape (", n_rows, ", ", n_cols, ")", sep="")
    # Raw binary of sparse Matrix in csc-format
    b_data = f.read()
    b_values = b_data[:sizeof_double * n_nonzero]
    b_pointers = b_data[sizeof_double * n_nonzero:]
    f.close()

    values = np.frombuffer(b_values)
    if is_verbose: print("Values ", values)

    pointers = np.frombuffer(b_pointers, dtype=np.uint)
    row_index = pointers[:n_nonzero]
    if is_verbose: print("Row index", row_index)
    col_pointer = pointers[n_nonzero:]
    if is_verbose: print("Column pointer", col_pointer)

    A = sp.csc_matrix((values, row_index, col_pointer), shape=(n_rows, n_cols)).transpose()
    A = A.tocsr() # This is efficient, linearly in n_nonzeros.
    if is_verbose: print(A.todense())
    return A

def remove_arma_tmp(path):
    import os
    os.remove(path)

# Tensor Gauss Quadrature
class tensorgauss:
    def __init__(self, deg, dim=4):
        self.deg = deg
        self.dim = dim
        p, w = np.polynomial.legendre.leggauss(deg)
        #print("Length", len(p), "\np", p, "w", w)
        self.N = deg**dim

        self.weights = np.ones(deg**dim)
        self.points = np.zeros((deg**dim, dim))
        l = 0
        k = np.zeros(self.dim, dtype=int)

        point_coord = np.array([row.flatten() for row in np.meshgrid(*([np.arange(self.deg)]*self.dim),
        indexing = "ij")]).T
        for k, dxRow in enumerate(point_coord):
            self.points[k] = p[dxRow]
            for dx in dxRow:
                self.weights[k] *= w[dx]

def stiffnessMatrix(
        # Mesh information ------------------------------------
        mesh,
        kernel,
        configuration,
        shape_derivative_data
    ):
    """ Computes the stiffness matrix corresponding to the nonlocal operator

    .. math::

      \mathcal{L}(\mathbf{u})(\mathbf{x}) = p.v. \int_{B_{\delta}(\mathbf{x}) \cap \widetilde{\Omega}}(\mathbf{C}_\delta(\mathbf{x}, \mathbf{y})  \mathbf{u}(\mathbf{x}) - \mathbf{C}_\delta(\mathbf{y}, \mathbf{x})\mathbf{u}(\mathbf{y}))  d\mathbf{y},

    on a given mesh. The input parameters are expected to be dictionaries. Find more information on the expected content in the `examples/Test2D/conf/testConfFull.py`. This function is a wrapper for the function par_assemble.

    :param mesh: Dictionary containing the mesh information (i.e. the keys ``"elements"``, ``"elementLabels"``, ``"vertices"``, and ``"vertexLabels"``). The values in the arrays ``"elements"``, ``"elementLabels"`` and ``"vertexLabels"`` are expected to be of datatype ``numpy.int``. The values in ``"vertices"`` are supposed to be of type ``numpy.float64``. Elements in the domain have positive labels. Elements in the nonlocal Dirichlet boundary have negative labels. For this purpose it does not matter which positive or negative number is used, and the kernels can depend on the element labels. The labels of the elements should be consistent with the vertex labels. In the case of Discontinuous Galerkin this means that the signs of the element labels and corresponding vertex labels coincide. In case of Continuous Galerkin this means that elements with negative label have only vertices with negative label.

    :param kernel: The kernel is assumed to be a dictionary containing the keys ``"function"``, and ``"outputdim"``. The value of ``"function"`` is a string which contains the name of the kernel (e.g. "constant"). You find all available options in the function *src/Cassemble.cpp/lookup_configuration()*. The key ``"outputdim"`` describes the output dimension of the kernel function. In case of a scalar valued kernel this would be 1. In case of the kernel "linearPrototypeMicroelasticField" this would be 2 (find an example in *examples/Test2D/testConfFull.py*). Note, that the value of ``"outputdim"`` does not define the output dimension of the kernel, but *describes* it. The value has to be consistent with the definition of the corresponding kernel in *src/model.cpp*.

    :param configuration: Dictionary containing the configuration (find an example in *examples/Test2D/conf/testConfFull.py*).

    :return: Matrix A in scipy csr-sparse format of shape K x K where K = nVerts * outdim (Continuous Galerkin) or K = nElems * (dim+1) * outdim (Discontinuous Galerkin).
    """
    tmstmp = timestamp()
    path_spAd = configuration.get("savePath", f"tmp_spAd_{tmstmp}")
    cdef string path_spAd_ = path_spAd.encode('UTF-8')

    #cdef long[:] neighbours = mesh["neighbours"].flatten()#nE*np.ones((nE*dVertex), dtype=int)
    #cdef int nNeighbours = mesh["neighbours"].shape[1]

    cdef long dim = mesh["elements"].shape[1] - 1

    cdef long nE = mesh["elements"].shape[0]
    cdef long nE_Omega = np.sum(mesh["elementLabels"] > 0)
    cdef long[:] elements = mesh["elements"].flatten()
    cdef long[:] elementLabels = mesh["elementLabels"].flatten()

    cdef long nV = mesh["vertices"].shape[0]
    cdef long nV_Omega = np.sum(mesh["vertexLabels"] > 0)
    cdef double[:] vertices = mesh["vertices"].flatten()
    cdef long[:] vertexLabels = mesh["vertexLabels"].flatten()

    cdef string model_kernel_ = kernel["function"].encode('UTF-8')
    cdef string integration_method_remote = configuration["approxBalls"]["method"].encode('UTF-8')
    cdef string integration_method_close = configuration.get("closeElements", configuration["approxBalls"]["method"]).encode('UTF-8')
    cdef int is_PlacePointOnCap_ = configuration["approxBalls"].get("is_PlacePointOnCap", 1)
    cdef int verbose = configuration.get("verbose", 0)
    if (verbose): print("stiffnessMatrix(): verbose mode.")

    cdef double [:] ptrPx = configuration["quadrature"]["outer"]["points"].flatten()
    cdef double [:] ptrPy = configuration["quadrature"]["inner"]["points"].flatten()
    cdef double [:] ptrdx = configuration["quadrature"]["outer"]["weights"].flatten()
    cdef double [:] ptrdy = configuration["quadrature"]["inner"]["weights"].flatten()
    cdef int nPx = configuration["quadrature"]["outer"]["weights"].shape[0]
    cdef int nPy = configuration["quadrature"]["inner"]["weights"].shape[0]

    cdef double [:] Pg
    cdef const double * ptrPg = NULL
    cdef double [:] dg
    cdef const double * ptrdg = NULL

    tensorGaussDegree = configuration["quadrature"].get("tensorGaussDegree", 0)
    if tensorGaussDegree != 0:
        quadgauss = tensorgauss(tensorGaussDegree)
        Pg = quadgauss.points.flatten()
        ptrPg = &Pg[0]
        dg = quadgauss.weights.flatten()
        ptrdg = &dg[0]
    else:
        if (verbose): print("stiffnessMatrix(): tensorGaussDegree not found or 0 (default).")

    cdef long [:] ZetaIndicator_indices
    cdef long [:] ZetaIndicator_indptr

    cdef long * ptrZetaIndicator_indices = NULL
    cdef long * ptrZetaIndicator_indptr = NULL
    cdef long nZeta

    try:
        nZeta = mesh["nZeta"]
        if nZeta > 0:
            ZetaIndicator_indices = np.array(mesh["ZetaIndicator"].indices, dtype=np.int)
            ZetaIndicator_indptr = np.array(mesh["ZetaIndicator"].indptr, dtype=np.int)
            ptrZetaIndicator_indices = &ZetaIndicator_indices[0]
            ptrZetaIndicator_indptr = &ZetaIndicator_indptr[0]
    except KeyError:
        nZeta = 0
        if (verbose): print("stiffnessMatrix(): Zeta not found. nZeta set to 0.")

    cdef int is_fullConnectedComponentSearch = configuration.get("is_fullConnectedComponentSearch", 0)
    if (verbose and not is_fullConnectedComponentSearch): print("stiffnessMatrix(): is_fullConnectedComponentSearch not found or 0 (default).")

    cdef int is_ShapeDerivative = configuration.get("is_ShapeDerivative", 0)
    if (verbose and not is_ShapeDerivative): print("stiffnessMatrix(): is_ShapeDerivative not found or 0 (default).")

    cdef int outdim_ = kernel["outputdim"]
    if mesh["outdim"] != kernel["outputdim"]:
        raise ValueError("The output dimension of the mesh has to be equal to the output dimension of the kernel.")

    cdef double fractional_s  = kernel.get("fractional_s", -1.0)
    if (verbose and (fractional_s==-1.0)): print("stiffnessMatrix(): fractional_s not found or -1.0 (default).")

    maxDiameter = mesh.get("diam", 0.0)
    if (verbose and (not maxDiameter)): print("stiffnessMatrix(): diam not found or 0.0 (default).")

    is_DG = configuration.get("ansatz", "CG") == "DG"
    if (verbose and (not maxDiameter)): print("stiffnessMatrix(): ansatz not found or CG (default).")

    cdef long K_Omega, K
    if is_DG:
        K = nE * (dim+1) * outdim_
        K_Omega = nE_Omega * (dim+1) * outdim_
    else:
        K = nV * outdim_
        K_Omega = nV_Omega * outdim_

    cdef double [:] state = np.zeros(nV)
    cdef double [:] adjoint = np.zeros(nV)
    if is_ShapeDerivative:
        state = shape_derivative_data['state']
        adjoint = shape_derivative_data['adjoint']

    # Compute Assembly --------------------------------
    start = time.time()
    Cassemble.par_assemble( "system".encode('UTF-8'), path_spAd_,
                            "".encode('UTF-8'),
                            K_Omega, K,
                            &elements[0], &elementLabels[0], &vertices[0], &vertexLabels[0],
                            nE , nE_Omega,
                            nV, nV_Omega,
                            &ptrPx[0], nPx, &ptrdx[0],
                            &ptrPy[0], nPy, &ptrdy[0],
                            kernel["horizon"]**2,
                            is_DG,
                            0,
                            &model_kernel_[0],
                            "".encode('UTF-8'),
                            &integration_method_remote[0],
                            &integration_method_close[0],
                            is_PlacePointOnCap_,
                            dim, outdim_,
                            ptrZetaIndicator_indices,
                            ptrZetaIndicator_indptr,
                            nZeta,
                            ptrPg, tensorGaussDegree, ptrdg,
                            maxDiameter,
                            fractional_s,
                             is_fullConnectedComponentSearch,
                             verbose,
                             is_ShapeDerivative, &state[0], &adjoint[0])

    total_time = time.time() - start
    if (verbose): print("Assembly Time\t", "{:1.2e}".format(total_time), " Sec")

    Ad = read_arma_spMat(path_spAd)
    if path_spAd == f"tmp_spAd_{tmstmp}":
        remove_arma_tmp(path_spAd)
    return Ad

def loadVector(
        # Mesh information ------------------------------------
        mesh,
        load,
        configuration
    ):
    """
    Computes load vector.
    The input parameters are expected to be dictionaries. Find more information
    on the expected content in the example/Test2D/testConfFull.py

    :param mesh: Dictionary containing the mesh information (i.e. the keys ``"elements"``, ``"elementLabels"``, ``"vertices"``, ``"vertexLabels"`` and ``"outdim"``), where ``"outdim"`` is the output-dimension of the kernel. In this function the output dimension is read from ``"outdim"``, as no kernel is expected. The load cannot yet depend on label information.

    :param load: Dictionary containing the load information (find an example in testConfFull.py).
    :param configuration: Dictionary containing the configuration (find an example in *examples/Test2D/conf/testConfFull.py*).

    :return: Vector f of shape K = nVerts * outdim (Continuous Galerkin) or K = nElems * (dim+1) * outdim (Discontinuous Galerkin)
    """
    tmstmp = timestamp()
    path_fd = configuration.get("savePath", f"tmp_fd_{tmstmp}")
    cdef string path_fd_ = path_fd.encode('UTF-8')

    #cdef long[:] neighbours = mesh["neighbours"].flatten()#nE*np.ones((nE*dVertex), dtype=int)
    #cdef int nNeighbours = mesh["neighbours"].shape[1]

    cdef long dim = mesh["elements"].shape[1] - 1
    cdef long nE = mesh["elements"].shape[0]
    cdef long nE_Omega = np.sum(mesh["elementLabels"] > 0)
    cdef long[:] elements = mesh["elements"].flatten()
    cdef long[:] elementLabels = mesh["elementLabels"].flatten()

    cdef long nV = mesh["vertices"].shape[0]
    cdef long nV_Omega = np.sum(mesh["vertexLabels"] > 0)
    cdef double[:] vertices = mesh["vertices"].flatten()
    cdef long[:] vertexLabels = mesh["vertexLabels"].flatten()

    cdef string model_load_ = load["function"].encode('UTF-8')
    cdef string integration_method_remote = configuration["approxBalls"]["method"].encode('UTF-8')
    cdef string integration_method_close = configuration.get("closeElements", configuration["approxBalls"]["method"]).encode('UTF-8')
    cdef int is_PlacePointOnCap_ = configuration["approxBalls"].get("is_PlacePointOnCap", 1)
    cdef int verbose = configuration.get("verbose", 0)

    cdef double [:] ptrPx = configuration["quadrature"]["outer"]["points"].flatten()
    #cdef double [:] ptrPy = configuration["quadrature"]["inner"]["points"].flatten()
    cdef double [:] ptrdx = configuration["quadrature"]["outer"]["weights"].flatten()
    #cdef double [:] ptrdy = configuration["quadrature"]["inner"]["weights"].flatten()
    cdef int nPx = configuration["quadrature"]["outer"]["weights"].shape[0]
    #cdef int nPy = configuration["quadrature"]["inner"]["weights"].shape[0]

    cdef double [:] Pg
    cdef const double * ptrPg = NULL
    cdef double [:] dg
    cdef const double * ptrdg = NULL
    tensorGaussDegree = 0

    #cdef long [:] ZetaIndicator_data
    cdef long [:] ZetaIndicator_indices
    cdef long [:] ZetaIndicator_indptr

    #cdef long * ptrZetaIndicator_data = NULL
    cdef long * ptrZetaIndicator_indices = NULL
    cdef long * ptrZetaIndicator_indptr = NULL
    cdef long nZeta=0

    try:
        nZeta = mesh["nZeta"]
        if nZeta > 0:
            ZetaIndicator_indices = np.array(mesh["ZetaIndicator"].indices, dtype=np.int)
            ZetaIndicator_indptr = np.array(mesh["ZetaIndicator"].indptr, dtype=np.int)
            ptrZetaIndicator_indices = &ZetaIndicator_indices[0]
            ptrZetaIndicator_indptr = &ZetaIndicator_indptr[0]
    except KeyError:
       nZeta = 0
       if (verbose): print("stiffnessMatrix(): Zeta not found. nZeta set to 0.")

    cdef int outdim_ = mesh["outdim"]

    is_DG = configuration.get("ansatz", "CG") == "DG"

    cdef long K_Omega, K
    if is_DG:
        K = nE * (dim+1) * outdim_
        K_Omega = nE_Omega * (dim+1) * outdim_
    else:
        K = nV * outdim_
        K_Omega = nV_Omega * outdim_

    # Compute Assembly --------------------------------
    start = time.time()
    Cassemble.par_assemble( "forcing".encode('UTF-8'), "".encode('UTF-8'),
                            path_fd_,
                            K_Omega, K,
                            &elements[0], &elementLabels[0],
                            &vertices[0], &vertexLabels[0],
                            nE, nE_Omega,
                            nV, nV_Omega,
                            &ptrPx[0], nPx, &ptrdx[0],
                            &ptrPx[0], nPx, &ptrdx[0],
                            0, is_DG, 0, "".encode('UTF-8'),
                            model_load_,
                            "".encode('UTF-8'), "".encode('UTF-8'), False,
                            dim, outdim_,
                            ptrZetaIndicator_indices,
                            ptrZetaIndicator_indptr,
                            nZeta,
                            NULL, 0, NULL, 0.0, -1.0, 0, verbose, 0, NULL, NULL)

    total_time = time.time() - start
    if (verbose): print("Assembly Time\t", "{:1.2e}".format(total_time), " Sec")

    fd = read_arma_mat(path_fd)[:,0]
    if path_fd == f"tmp_fd_{tmstmp}":
        remove_arma_tmp(path_fd)
    return fd


def get_vertexLabel(elements, elementLabels, vertices):
    nV = vertices.shape[0]
    vertexLabels = np.zeros(nV, dtype=np.int)
    label_list = np.sort(np.unique(elementLabels))
    #print(label_list)
    for label in label_list[::-1]:
        if label:
            label_idx = np.where(elementLabels == label)
            vertex_idx = np.unique(elements[label_idx].ravel())
            vertexLabels[vertex_idx] = label
    return vertexLabels

def get_diam(elements, vertices):
    dim = vertices.shape[1]
    nVerts = dim + 1
    diam = 0.0
    T = vertices[elements]
    for t in T:
        for k in range(nVerts):
            diff = t[k] - t[(k+1)%nVerts]
            dist = np.sqrt(diff**2)[0]
            if dist > diam:
                diam = dist
    return diam

def meshFromArrays(elements, elementLabels, vertices, outputdim=1):
    vertexLabels = get_vertexLabel(elements, elementLabels, vertices)
    diam = get_diam(elements, vertices)

    mesh = {
        "dim": elements.shape[1] - 1,
        "nE": elements.shape[0],
        "nV": vertices.shape[0],
        "nE_Omega": np.sum(elementLabels > 0),
        "nV_Omega": np.sum(vertexLabels > 0),
        "elements": np.array(elements, dtype=np.int),
        "elementLabels": elementLabels,
        "vertices": np.array(vertices),
        "vertexLabels": vertexLabels,
        "diam": diam,
        "outdim": outputdim
    }
    return mesh

def setSolutionLabels(mesh, ansatz="CG", outputdim=1):
    mesh["ansatz"] = ansatz
    mesh["outputdim"] = outputdim
    if ansatz == "CG":
        mesh["solutionLabels"] = np.repeat(mesh["vertexLabels"], outputdim)
        mesh["K"] = mesh["nV"]*outputdim
        mesh["K_Omega"] = mesh["nV_Omega"]*outputdim
    if ansatz == "DG":
        mesh["solutionLabels"] = np.repeat(mesh["elementLabels"], 3*outputdim)
        mesh["K"] = 3*mesh["nE"]*outputdim
        mesh["K_Omega"] = 3*mesh["nE_Omega"]*outputdim


def stiffnessMatrix_fromArray(
        elements, elementLabels, vertices,
        kernel,
        configuration,
        shape_derivative_data
    ):
    """ Computes the stiffness matrix corresponding to the nonlocal operator

    .. math::

      \mathcal{L}(\mathbf{u})(\mathbf{x}) = p.v. \int_{B_{\delta}(\mathbf{x}) \cap \widetilde{\Omega}}(\mathbf{C}_\delta(\mathbf{x}, \mathbf{y})  \mathbf{u}(\mathbf{x}) - \mathbf{C}_\delta(\mathbf{y}, \mathbf{x})\mathbf{u}(\mathbf{y}))  d\mathbf{y},

    on a given mesh. This function is a wrapper around ``stiffnessMatrix()``.

    :param elements: Numpy array containing the information about the elements. The values in the array ``"elements"`` are expected to be of datatype ``numpy.int``.


    :param elementLabels: The values in the array ``"elementLabels"`` are expected to be of datatype ``numpy.int``. Elements in the domain have positive labels. Elements in the nonlocal Dirichlet boundary have negative labels. For this purpose it does not matter which positive or negative number is used, and the kernels can depend on the element labels. The routine does not compute contributions of elements with label zero, i.e. they are ignored. The label zero can therefore be added to connect disconnected domains. The labels of the vertices are automatically derived from the element labels. In the case of Discontinuous Galerkin this means that the signs of the element labels and corresponding vertex labels coincide. In case of Continuous Galerkin this means that a vertex gets the smallest label of the non-zero labels of the adjacent elements.

    :param vertices: The values in the array ``"vertices"`` are supposed to be of type ``numpy.float64``.

    :param kernel: The kernel is assumed to be a dictionary containing the keys ``"function"``, and ``"outputdim"``. The value of ``"function"`` is a string which contains the name of the kernel (e.g. "constant"). You find all available options in the function *src/Cassemble.cpp/lookup_configuration()*. The key ``"outputdim"`` describes the output dimension of the kernel function. In case of a scalar valued kernel this would be 1. In case of the kernel "linearPrototypeMicroelasticField" this would be 2 (find an example in *examples/Test2D/testConfFull.py*). Note, that the value of ``"outputdim"`` does not define the output dimension of the kernel, but *describes* it. The value has to be consistent with the definition of the corresponding kernel in *src/model.cpp*.

    :param configuration: Dictionary containing the configuration (find an example in *examples/Test2D/conf/testConfFull.py*).

    :return mesh: Dictionary containing the mesh information. Can also be used for ``stiffnessMatrix()`` and ``loadVector()``.

    :return A: Matrix in scipy csr-sparse format of shape K x K where K = nVerts * outdim (Continuous Galerkin) or K = nElems * (dim+1) * outdim (Discontinuous Galerkin).
    """

    mesh = meshFromArrays(elements, elementLabels, vertices, kernel["outputdim"])
    setSolutionLabels(mesh, configuration["ansatz"], kernel["outputdim"])
    A = stiffnessMatrix(mesh, kernel, configuration, shape_derivative_data)
    return mesh, A

def assemble(
        # Mesh information ------------------------------------
        mesh,
        Px,
        Py,
        # Weights for quadrature rule
        dx,
        dy,
        double delta,
        path_spAd=None,
        path_fd=None,
        compute="systemforcing", # "forcing", "system"
        model_kernel="constant",
        model_f = "constant",
        integration_method = "retriangulate",
        is_PlacePointOnCap = 1,
        tensorGaussDegree=0
    ):
    """
    Computes stiffness matrix and load vector. This function is deprecated. It expects a class mesh
    and some other input parameters. It returns a load vector of shape K_Omega, instead of K, which assumes
    that the nodes are ordered such that Dirichlet nodes appear as last values. In case of Continuous Galerkin
    this means that vertices with positve label appear last. In case of Discontinuous Galerkin this means
    that elements with negative label appear last.

    :return: Vector f of shape K = nVerts * outdim (Continuous Galerkin) or K = nElems * (dim+1) * outdim
    (Discontinuous Galerkin)
    """
    print("This function is deprecated. Please use stiffnessMatrix()")
    #raise KeyboardInterrupt
    cdef int verbose  = 1
    try:
        verbose = mesh.verbose
    except AttributeError:
        print("Verbose Mode")

    is_tmpAd = False
    if path_spAd is None:
        is_tmpAd = True
        path_spAd = "tmp_spAd"
    cdef string path_spAd_ = path_spAd.encode('UTF-8')
    Ad = None

    is_tmpfd = False
    if path_fd is None:
        is_tmpfd = True
        path_fd = "tmp_fd"
    cdef string path_fd_ = path_fd.encode('UTF-8')
    fd = None

    #cdef long[:] neighbours = mesh.neighbours.flatten()#nE*np.ones((nE*dVertex), dtype=int)
    #cdef int nNeighbours = mesh.neighbours.shape[1]
    cdef long[:] elements = mesh.elements.flatten()
    cdef long[:] elementLabels = mesh.elementLabels.flatten()
    cdef double[:] vertices = mesh.vertices.flatten()
    cdef long[:] vertexLabels = mesh.vertexLabels.flatten()
    #cdef double[:] ptrAd = Ad
    cdef double[:] ptrfd = fd
    cdef string model_kernel_ = model_kernel.encode('UTF-8')
    cdef string model_f_ = model_f.encode('UTF-8')
    cdef string integration_method_remote = integration_method.encode('UTF-8')
    cdef string integration_method_close = integration_method.encode('UTF-8')

    cdef string compute_system_ = "system".encode('UTF-8')
    cdef string compute_forcing_ = "forcing".encode('UTF-8')
    cdef int is_PlacePointOnCap_ = is_PlacePointOnCap


    cdef double [:] ptrPx = Px.flatten()
    cdef double [:] ptrPy = Py.flatten()
    cdef double [:] ptrdx = dx.flatten()
    cdef double [:] ptrdy = dy.flatten()

    cdef long nZeta

    try:
        outdim = mesh.outdim
    except AttributeError:
        outdim = 1
        if (verbose): print("Mesh out dim not found. outdim set to 1.")

    cdef double [:] Pg
    cdef const double * ptrPg = NULL
    cdef double [:] dg
    cdef const double * ptrdg = NULL

    cdef double maxDiameter = 0.0
    try:
        maxDiameter = mesh.diam
    except AttributeError:
        if (verbose): print("Element diameter not found. maxDiameter set to 0.0.")

    if tensorGaussDegree != 0:
        quadgauss = tensorgauss(tensorGaussDegree)
        Pg = quadgauss.points.flatten()
        ptrPg = &Pg[0]
        dg = quadgauss.weights.flatten()
        ptrdg = &dg[0]

    cdef int is_fullConnectedComponentSearch = 0
    cdef int is_ShapeDerivative = 0
    try:
        is_fullConnectedComponentSearch = mesh.is_fullConnectedComponentSearch
    except AttributeError:
        if (verbose): print("is_fullConnectedComponentSearch not found. Set to 0.")

    # Compute Assembly --------------------------------
    if (compute=="system" or compute=="systemforcing"):
        start = time.time()
        Cassemble.par_assemble( compute_system_, path_spAd_, path_fd_, mesh.K_Omega, mesh.K,
                            &elements[0], &elementLabels[0], &vertices[0], &vertexLabels[0],
                            mesh.nE , mesh.nE_Omega, mesh.nV, mesh.nV_Omega,
                            &ptrPx[0], Px.shape[0], &ptrdx[0],
                            &ptrPy[0], Py.shape[0], &ptrdy[0],
                            delta**2,
                            mesh.is_DiscontinuousGalerkin,
                            mesh.is_NeumannBoundary,
                            &model_kernel_[0],
                            &model_f_[0],
                            &integration_method_remote[0],
                            &integration_method_close[0],
                            is_PlacePointOnCap_,
                            mesh.dim, outdim, NULL, NULL, 0,
                            ptrPg, tensorGaussDegree, ptrdg, maxDiameter, -1.0, is_fullConnectedComponentSearch,
                            verbose, is_ShapeDerivative, NULL, NULL)

        total_time = time.time() - start
        if (verbose): print("Assembly Time\t", "{:1.2e}".format(total_time), " Sec")

        Ad = read_arma_spMat(path_spAd)[:mesh.K_Omega]
        if is_tmpAd:
            remove_arma_tmp(path_spAd)

    if (compute=="forcing" or compute =="systemforcing"):
        if (verbose): print("")
        Cassemble.par_assemble( compute_forcing_, path_spAd_, path_fd_, mesh.K_Omega, mesh.K,
                            &elements[0], &elementLabels[0], &vertices[0], &vertexLabels[0],
                            mesh.nE , mesh.nE_Omega, mesh.nV, mesh.nV_Omega,
                            &ptrPx[0], Px.shape[0], &ptrdx[0],
                            &ptrPy[0], Py.shape[0], &ptrdy[0],
                            delta**2,
                            mesh.is_DiscontinuousGalerkin,
                            mesh.is_NeumannBoundary,
                            &model_kernel_[0],
                            &model_f_[0],
                            &integration_method_remote[0],
                            &integration_method_close[0],
                            is_PlacePointOnCap_,
                            mesh.dim, outdim,  NULL, NULL, 0,
                            ptrPg, tensorGaussDegree, ptrdg, maxDiameter, -1.0, is_fullConnectedComponentSearch,
                            verbose, is_ShapeDerivative, NULL, NULL)

        fd = read_arma_mat(path_fd)[:mesh.K_Omega,0]
        if is_tmpfd:
            remove_arma_tmp(path_fd)

    return Ad, fd

def evaluateMass(
      # Mesh information ------------------------------------
            mesh,
            ud,
            Px,
            # Weights for quadrature rule
            dx
        ):
    vd = np.zeros(mesh.K)
    cdef double[:] ptrvd = vd
    cdef double[:] ptrud = ud.flatten()

    cdef long[:] elements = mesh.elements.flatten()
    cdef long [:] elementLabels = mesh.elementLabels.flatten()
    cdef double[:] vertices = mesh.vertices.flatten()
    cdef long [:] vertexLabels = mesh.vertexLabels.flatten()

    cdef double[:] ptrPx = Px.flatten()
    cdef double[:] ptrdx = dx.flatten()

    cdef int isDG = mesh.is_DiscontinuousGalerkin;

    cdef long outdim = 1
    try:
        outdim = mesh.outdim
    except AttributeError:
        pass

    Cassemble.par_evaluateMass(
            &ptrvd[0],
            &ptrud[0],
            &elements[0],
            &elementLabels[0],
            &vertices[0],
            &vertexLabels[0],
            mesh.K_Omega,
            mesh.nE_Omega,
            Px.shape[0], &ptrPx[0], &ptrdx[0], mesh.dim, outdim, isDG)
    return vd

def constructAdjaciencyGraph(long[:,:] elements, ncommon = 1):
    print("Constructing adjaciency graph...")
    print("This function is deprecated, and its use is not required anymore.")
    nE = elements.shape[0]
    nV = np.max(elements)+1
    cdef int dVerts = elements.shape[1]
    cdef int dim = dVerts-1

    #neigs = np.ones((nE, dim+1), dtype=np.int)*nE
    grph_elements = sparse.lil_matrix((nV, nE), dtype=np.int)

    for Tdx, Vdx in enumerate(elements):
        for d in range(dVerts):
            grph_elements[Vdx[d], Tdx] = 1
    #grph_neigs = ((grph_elements.transpose() @ grph_elements) == dim)
    grph_neigs = ((grph_elements.transpose() @ grph_elements) > (ncommon-1))

    rowsum = np.sum(grph_neigs, axis= 0)
    nCols = np.max(rowsum)
    neigs = np.ones((nE, nCols), dtype=np.int)*nE
    elemenIndices, neighbourIndices = grph_neigs.nonzero()

    neigs[elemenIndices[0],0] = neighbourIndices[0]
    cdef int colj = 0
    cdef int k

    for k in range(1, len(elemenIndices)):
        colj *= ((elemenIndices[k-1]-elemenIndices[k])==0)
        colj += ((elemenIndices[k-1]-elemenIndices[k])==0)
        neigs[elemenIndices[k], colj] =  neighbourIndices[k]
    return neigs

# DEBUG Helpers - -----------------------------------------------------------------------------------------------------
from Cassemble cimport get_div
def py_get_div(
        double [:] TE):
        cdef double div[6]
        get_div(&TE[0], &div[0]);
        return div

from Cassemble cimport method_retriangulate, method_retriangulateInfty
def py_retriangulate(
        double [:] x_center,
        double [:] TE,
        double delta,
        int is_placePointOnCaps,
        pp
    ):

    TriangleList =  np.zeros(9*3*2)
    cdef:
        double sqdelta = pow(delta,2)
        double [:] cTriangleList = TriangleList
    Rdx = method_retriangulate(&x_center[0], &TE[0], sqdelta, &cTriangleList[0], is_placePointOnCaps)

    return Rdx, TriangleList
def py_retriangulateInfty(
        double [:] x_center,
        double [:] TE,
        double delta,
        int is_placePointOnCaps,
        pp
    ):

    TriangleList =  np.zeros(9*3*2)
    cdef:
        double sqdelta = pow(delta,2)
        double [:] cTriangleList = TriangleList
    Rdx = method_retriangulateInfty(&x_center[0], &TE[0], sqdelta, &cTriangleList[0], is_placePointOnCaps)
    return Rdx, TriangleList

from Cassemble cimport toRef
def py_toRef(
    double [:] TE,
    double [:] phys_x):
    ref_p = np.zeros(2)
    cdef double [:] cref_p = ref_p
     # void toRef(const double * E, const double * phys_x, double * ref_p);
    toRef(&TE[0], &phys_x[0], &cref_p[0]);
    return ref_p

from Cassemble cimport toPhys
def py_toPhys(
    double [:] TE,
    double [:] p):
    out_x = np.zeros(2)
    cdef double [:] cout_x = out_x
     # void toPhys(const double * E, const double * p, int dim, double * out_x)
    toPhys(&TE[0], &p[0], 2, &cout_x[0]);
    return out_x

from Cassemble cimport solve2x2
def py_solve2x2(
    double [:] A,
    double [:] b
    ):
    x = np.zeros(2)
    cdef double [:] cx = x
    # void solve2x2(const double * A, const double * b, double * x)
    solve2x2(&A[0], &b[0], &cx[0])
    return x
"""
from Cassemble cimport retriangulate
from Cassemble cimport toRef, model_basisFunction
import matplotlib
from libcpp.queue cimport queue
import matplotlib.pyplot as plt

def py_retriangulate(
        double [:] x_center,
        double [:] TE,
        double delta,
        int is_placePointOnCaps,
        pp
    ):

    TriangleList =  np.zeros(9*3*2)
    cdef:
        double sqdelta = pow(delta,2)
        double [:] cTriangleList = TriangleList
    Rdx = retriangulate(&x_center[0], &TE[0], sqdelta, &cTriangleList[0], is_placePointOnCaps)

    return Rdx, TriangleList

def py_check_par_assemble(
        # Mesh information ------------------------------------
        mesh,
        double [:,:] py_P,
        # Weights for quadrature rule
        double [:] dx,
        double delta,
        **kwargs
    ):

    pp = kwargs.get("pp", None)

    K = mesh.nE

    py_fd = np.zeros(K).flatten("C")

    cdef:
        long nE = mesh.nE
        long nV = mesh.nV
        long nE_Omega = mesh.nE_Omega
        long nV_Omega = mesh.nV_Omega

        int nP = py_P.shape[0] # Does not differ in inner and outer integral!
        long [:] Triangles = (np.array(mesh.triangles, int)).flatten("C")
        double [:] Verts = (np.array(mesh.vertices, float)).flatten("C")
        double[:] P = (np.array(py_P, float)).flatten("C")

        # Cython interface of C-aligned arrays of solution and right side
        double[:] fd = py_fd
        # List of neighbours of each triangle, Neighbours[Tdx] returns row with Neighbour indices
        long[:] Neighbours = mesh.nE*np.ones((mesh.nE*dVertex), dtype=int)
        # Squared interaction horizon
        double sqdelta = pow(delta,2)

    # Setup adjaciency graph of the mesh --------------------------
    neigs = []
    for aTdx in range(mesh.nE):
        neigs = get_neighbour(mesh.nE, dVertex, &Triangles[0], &Triangles[(dVertex+1)*aTdx])
        n = len(neigs)
        for i in range(n):
            Neighbours[4*aTdx + i] = neigs[i]

    start = time.time()

    cdef int  h=0
    ## General Loop Indices ---------------------------------------
    cdef int bTdx=0

    ## Breadth First Search --------------------------------------
    cdef visited = np.zeros(mesh.nE)#(int *) malloc(nE*sizeof(int));

    ## Loop index of current outer triangle in BFS
    cdef int sTdx=0
    ## Queue for Breadth first search
    cdef queue[int] queue
    ## List of visited triangles
    cdef long *NTdx
    ## Determinant of Triangle a and b.
    cdef  double aTdet, bTdet, termarea
    ## Vector containing the coordinates of the vertices of a Triangle
    cdef double aTE[2*3]
    cdef double bTE[2*3]
    ## Integration information ------------------------------------
    ##// Loop index of basis functions
    cdef int a=0, b=0, aAdxj =0
    cdef long labela=0, labelb=0
    #cdef double [:] Verts = mesh.vertices.flatten("C")
    #cdef long [:] Triangles = mesh.triangles.flatten("C")

    for aTdx in [7,16, 41]:#, 95, 98, 100]:#(mesh.nE-10, mesh.nE):
        fig = plt.figure()  # create a figure object
        ax = fig.add_subplot(1, 1, 1)
        ## Prepare Triangle information aTE and aTdet ------------------
        ## Copy coordinates of Triange a to aTE.
        ## this is done fore convenience only, actually those are unnecessary copies!
        print(aTdx)
        for j in range(2):
            aTE[2*0+j] = Verts[2*Triangles[(dVertex+1)*aTdx+1] + j]
            aTE[2*1+j] = Verts[2*Triangles[(dVertex+1)*aTdx+2] + j]
            aTE[2*2+j] = Verts[2*Triangles[(dVertex+1)*aTdx+3] + j]

        ## compute Determinant
        aTdet = absDet(aTE)
        labela = Triangles[(dVertex+1)*aTdx]

        ## Of course some uneccessary computation happens but only for some verticies of thos triangles which lie
        ## on the boundary. This saves us from the pain to carry the information (a) into the integrator compute_f.

        ## BFS -------------------------------------------------------------
        ## Intialize search queue with current outer triangle
        queue.push(aTdx)
        ## Initialize vector of visited triangles with 0
        visited.fill(0)
        h=0
        ## Check whether BFS is over.
        while ( not queue.empty()):
            h+=1
            ## Get and delete the next Triangle index of the queue. The first one will be the triangle aTdx itself.
            sTdx = queue.front()
            queue.pop()
            ## Get all the neighbours of sTdx.
            NTdx =  &Neighbours[4*sTdx]
            ## Run through the list of neighbours. (4 at max)
            for j in range(4):
                ## The next valid neighbour is our candidate for the inner Triangle b.
                bTdx = NTdx[j]

                ## Check how many neighbours sTdx has. It can be 4 at max. (Itself, and the three others)
                ## In order to be able to store the list as contiguous array we fill up the empty spots with the number nE
                ## i.e. the total number of Triangles (which cannot be an index).
                if (bTdx < mesh.nE):

                    ## Prepare Triangle information bTE and bTdet ------------------
                    ## Copy coordinates of Triange b to bTE.
                    ## again this is done fore convenience only, actually those are unnecessary copies!
                    for i in range(2):
                        bTE[2*0+i] = Verts[2*Triangles[(dVertex+1)*bTdx+1] + i]
                        bTE[2*1+i] = Verts[2*Triangles[(dVertex+1)*bTdx+2] + i]
                        bTE[2*2+i] = Verts[2*Triangles[(dVertex+1)*bTdx+3] + i]

                    bTdet = absDet(bTE)
                    labelb = Triangles[(dVertex+1)*bTdx]
                    ## Check wheter bTdx is already visited.
                    if (visited[bTdx]==0):

                        ## Assembly of matrix ---------------------------------------
                        ## Compute integrals and write to buffer
                        #termarea = 1
                        termarea = compute_area(aTE, aTdet, labela, bTE, bTdet, labelb, P, nP, dx, sqdelta, pp, ax, aTdx, bTdx, h)
                        ## If bT interacts it will be a candidate for our BFS, so it is added to the queue
                        if (termarea > 0):
                            queue.push(bTdx)
                            fd[aTdx] += termarea
                    visited[bTdx] = 1
        if pp is not None:
            plt.savefig(pp, format='pdf')
            plt.close()

    total_time = time.time() - start
    if (verbose): print("Assembly Time\t", "{:1.2e}".format(total_time), " Sec")


    return py_fd

cdef absDet(double [:] E):
    cdef double [:,:] M = np.zeros((2,2))
    cdef int i=0
    for i in range(2):
        M[i][0] = E[2*1+i] - E[2*0+i]
        M[i][1] = E[2*2+i] - E[2*0+i]

    return np.abs(M[0][0]*M[1][1] - M[0][1]*M[1][0])

cdef toPhys(double [:] E, double * p, double * out_x):
    cdef int i=0
    for i in range(2):
        out_x[i] = (E[2*1+i] - E[2*0+i])*p[0] + (E[2*2+i] - E[2*0+i])*p[1] + E[2*0+i]

cdef double compute_area(double [:] aTE, double aTdet, long labela, double [:] bTE, double bTdet, long labelb, double [:] P, int nP, double [:] dx, double sqdelta, pp, ax, aTdx, bTdx, h):
    cdef:
        double areaTerm=0.0
        int rTdx, Rdx, s_max=15
        double * x
        double physical_quad[2], reference_quad[2], psi_value[3]
        double [:] reTriangle_list = np.zeros([9*3*2])

    x = &P[2*15]
    toPhys(&aTE[0], x, &physical_quad[0])
    py_IntegrationPoint = np.array(physical_quad)
    TE = np.array(bTE)
    TE = TE.reshape((3,2))
    if pp is not None:
        ax.fill(TE[:, 0], TE[:, 1], color="b", alpha=1, fill=False,  lw=.1)
        ax.scatter(py_IntegrationPoint[0], py_IntegrationPoint[1], c = "blue", s=2)
    is_placePointOnCaps = 1
    Rdx = retriangulate(physical_quad, &bTE[0], sqdelta, &reTriangle_list[0], is_placePointOnCaps)
    for rTdx in range(Rdx):
        areaTerm += absDet(reTriangle_list[2*3*rTdx:])/2
        TE = np.array(reTriangle_list[2*3*rTdx:(2*3*rTdx + 6)])
        TE = TE.reshape((3,2))
        if pp is not None:
            plt.title(aTdx)
            plt.gca().set_aspect('equal')
            ax.fill(TE[:, 0], TE[:, 1], color="blue", alpha=.3)
            ax.fill(TE[:, 0], TE[:, 1], color="blue", fill=False,  lw=.7, alpha=.3)
            if bTdx == aTdx:
                ax.fill(TE[:, 0], TE[:, 1], color="blue", lw=.7, alpha=.3)
                physical_quad_List = np.zeros((nP,2))
                psi_value_List = np.zeros(nP)
                for i in range(0):
                    toPhys(&reTriangle_list[2*3*rTdx], &P[2*i], &physical_quad[0])
                    py_physical_quad = np.array(physical_quad)
                    ax.scatter(py_physical_quad[0], py_physical_quad[1], c="black", s=.1, alpha=.5)
                    ax.annotate(i, py_physical_quad,alpha=.7, fontsize=2)
                    toRef(&bTE[0], &physical_quad[0], &reference_quad[0])
                    model_basisFunction(&reference_quad[0], &psi_value[0])
                    #print(np.array(psi_value))
                    #ax.scatter(py_physical_quad[0], py_physical_quad[1], c="black", s=psi_value[1]*s_max, alpha=.9)
                    #ax.tricontour(py_physical_quad[0], py_physical_quad[1], Z=psi_value[1])
                    py_bTE = np.array(bTE)
                    py_bTE = py_bTE.reshape((3,2))
                    ax.scatter(py_bTE[:,0], py_bTE[:,1], s=3, alpha=1)
                    for k in range(3):
                        ax.annotate("Psi"+str(k) ,py_bTE[k], fontsize=3, alpha=.1)
                    physical_quad_List[i] = py_physical_quad
                    psi_value_List[i] =  psi_value[0]
                #ax.tricontourf(physical_quad_List[:,0], physical_quad_List[:,1], psi_value_List,  5, cmap=plt.cm.get_cmap('rainbow'), vmin=0, vmax=1)
                #ax, _ = matplotlib.colorbar.make_axes(plt.gca(), shrink=.7)
                #matplotlib.colorbar.ColorbarBase(ax, cmap=plt.cm.get_cmap('rainbow'),
                #                             norm=matplotlib.colors.Normalize(vmin=0, vmax=1))
    if pp is not None:
        ax.annotate(bTdx, baryCenter(bTE), fontsize=5)

    return areaTerm

cdef double [:] baryCenter(double [:] E):
    cdef:
        int i=0
        double [:] bary = np.zeros(2)
    for i in range(3):
        bary[0] += E[2*i]
        bary[1] += E[2*i+1]
    bary[0] = bary[0]/3
    bary[1] = bary[1]/3
    return bary

"""
# [END] DEBUG Helpers - ------------------------------------------------------------------------------------------------
import numpy as np
import scipy
import scipy.sparse as ss

import fenics


def l2_norm_2d(x):
    return np.sqrt(x[0]**2 + x[1]**2)


def basis(x):
    return np.array([1 - x[0] - x[1], x[0], x[1]])


def get_function_value(point, values):
    return np.dot(basis(point), values)


def get_transformation_matrix(vertices):
    transformation = np.array([vertices[1] - vertices[0], vertices[2] - vertices[0]]).transpose()
    det_transformation = transformation[0, 0] * transformation[1, 1] - transformation[1, 0] * transformation[0, 1]
    inverse_transformation = 1. / det_transformation * np.array([[transformation[1, 1], -transformation[0, 1]],
                                                                 [-transformation[1, 0], transformation[0, 0]]])
    return transformation, det_transformation, inverse_transformation


def convert_nonlocal_function_into_fenics_function_cg(mesh, u, theta, vertexLabels):
    V = fenics.FunctionSpace(mesh, 'CG', 1)
    num_vertices = mesh.num_vertices()
    u_sd_temp = np.zeros(num_vertices, dtype=float)
    indices_theta = 0
    indices_u = 0
    for i in range(0, num_vertices):
        if vertexLabels[i] > 0:
            u_sd_temp[i] = u[indices_u]
            indices_u += 1
        else:
            u_sd_temp[i] = theta[indices_theta]
            indices_theta += 1
    u_sd = fenics.Function(V)
    u_sd.vector().set_local(u_sd_temp)
    return u_sd


def convert_nonlocal_function_into_fenics_function_dg(mesh, u, g, elementLabels):
    V = fenics.FunctionSpace(mesh, 'DG', 1)
    num_elements = mesh.num_cells()
    u_sd_temp = np.zeros(3*num_elements, dtype=float)
    indices_g = 0
    indices_u = 0
    for i in range(0, num_elements):
        if elementLabels[i] > 0:
            u_sd_temp[3*i] = u[indices_u]
            u_sd_temp[3*i + 1] = u[indices_u + 1]
            u_sd_temp[3*i + 2] = u[indices_u + 2]
            indices_u += 3
        else:
            u_sd_temp[3*i] = g[indices_g]
            u_sd_temp[3*i + 1] = g[indices_g + 1]
            u_sd_temp[3*i + 2] = g[indices_g + 2]
            indices_g += 3
    u_sd = fenics.Function(V)
    u_sd.vector().set_local(u_sd_temp)
    return u_sd


def petsc_to_csr_scipy_matrix(petsc_matrix):
    indptr, indices, data = petsc_matrix.mat().getValuesCSR()
    n = petsc_matrix.size(0)
    m = petsc_matrix.size(1)
    csr_matrix = scipy.sparse.csr_matrix((data, indices, indptr), shape=(n, m))
    csr_matrix.eliminate_zeros()
    return csr_matrix


def mesh_not_admissible(coordinates, indices_boundary_vertices, tol):
    num_vertices = coordinates.shape[0]
    vertex_indices = range(num_vertices)
    omega_vertex_indices = list(set(vertex_indices) - set(indices_boundary_vertices))

    for vertex_index in omega_vertex_indices:
        vertex = coordinates[vertex_index]
        if min(vertex[0], vertex[1]) < 0 - tol or max(vertex[0], vertex[1]) > 1 + tol:
            return 1
    else:
        return 0

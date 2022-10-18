import numpy as np
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import cg
import scipy.sparse as scs
import time

import fenics
import nlfem

import mesh_data
import helper
import bfgs
import results_data


class NonlocalShapeProblem:
    def __init__(self, conf):
        fenics.parameters["reorder_dofs_serial"] = False  # Order of vertices and dofs coincide

        self.nlfem_shape_conf = conf['nlfem_shape_conf']
        self.kernel_dependent = conf['kernel_dependent']

        self.ansatz = conf['nlfem_conf']['ansatz']
        self.horizon = conf['kernel']['horizon']
        self.kernel = conf['kernel']
        self.kernel_function = conf['kernel']['function']

        self.nlfem_conf = conf['nlfem_conf']
        self.source = conf['source']
        self.file_name = conf['init_shape']
        self.file_name_current = conf['init_shape']
        self.boundary_label = int(conf['boundary_label'])
        self.interface_label = int(conf['interface_label'])
        self.alpha = conf['alpha']
        self.memory_length = conf['memory_length']

        self.c_per = conf['c_per']
        self.remesh = conf['remesh']
        self.element_size = conf['element_size']

        self.lbd = conf['lambda']
        self.mu_min = conf['mu_min']
        self.mu_max = conf['mu_max']
        self.nu = conf['nu']
        self.adapt_mu_max = conf['adapt_mu_max']
        self.lame_threshold_up = conf['lame_threshold_up']
        self.lame_threshold_down = conf['lame_threshold_down']
        self.adapt_up = conf['adapt_up']
        self.adapt_down = conf['adapt_down']

        self.outer_points = conf['nlfem_conf']['quadrature']['outer']['points']
        self.outer_weights = conf['nlfem_conf']['quadrature']['outer']['weights']
        self.outer_basis_values = helper.basis(self.outer_points.transpose())

        self.inner_points = conf['nlfem_conf']['quadrature']['inner']['points']
        self.inner_weights = conf['nlfem_conf']['quadrature']['inner']['weights']
        self.inner_basis_values = helper.basis(self.inner_points.transpose())

        self.target_mesh_data = mesh_data.MeshData(conf["target_shape"], self.interface_label, self.boundary_label)
        self.current_mesh_data = mesh_data.MeshData(conf["init_shape"], self.interface_label, self.boundary_label)

        num_vertices = self.current_mesh_data.mesh.num_vertices()
        self.remesh_counter = 1

        self.armijo_const = conf["armijo_const"]
        self.l_bfgs = conf["l_bfgs"]
        if self.l_bfgs:
            self.l_bfgs_algorithm = bfgs.BFGS(num_vertices, conf["memory_length"])

        self.u_bar = self.solve_state(self.target_mesh_data.mesh, self.target_mesh_data.subdomains, self.c_per)
        self.save_results = conf["save_results"]
        if self.save_results:
            self.results = results_data.Results(conf, self.current_mesh_data.interface,
                                                self.current_mesh_data.subdomains, self.target_mesh_data.interface)

            self.results.save_target_data(self.u_bar, self.target_mesh_data.subdomains)

    def solve_state(self, mesh, subdomains, c_per):
        elements, vertices, elementLabels = mesh_data.convert_mesh(mesh, subdomains)
        empty_dict = {}
        self.nlfem_conf['is_ShapeDerivative'] = 0
        mesh_nl, A = nlfem.stiffnessMatrix_fromArray(elements, elementLabels, vertices, self.kernel, self.nlfem_conf,
                                                     empty_dict)

        vertexLabels = mesh_nl["vertexLabels"]

        if self.ansatz == "DG":
            V = fenics.FunctionSpace(mesh, 'DG', 1)
            v = fenics.TestFunction(V)
            dx = fenics.Measure('dx', domain=mesh, subdomain_data=subdomains)
            l = self.source[0] * v * dx(1) + self.source[1] * v * dx(2) + 0 * v * dx(3)
            b = fenics.PETScVector()
            fenics.assemble(l, tensor=b)
            b_vec = b.get_local()

            nodeLabels = np.repeat(elementLabels, 3)

            if c_per:
                u = fenics.TrialFunction(V)
                a_per = fenics.dot(fenics.grad(u), fenics.grad(v)) * dx
                A_per = fenics.PETScMatrix()
                fenics.assemble(a_per, tensor=A_per)
                A_per = helper.petsc_to_csr_scipy_matrix(A_per)
                A = A + A_per

            A_O = A[nodeLabels > 0][:, nodeLabels > 0]
            A_I = A[nodeLabels > 0][:, nodeLabels < 0]
            f_O = b_vec[nodeLabels > 0]

            def g2(x):
                return 0.0

            num_boundary_elements = np.sum(elementLabels < 0)
            g = np.zeros(3 * num_boundary_elements)

            for index, element in enumerate(elements[elementLabels < 0]):
                for vertex_index in range(3):
                    vertex_coordinates = vertices[element[vertex_index]]
                    g[3 * index + vertex_index] = g2(vertex_coordinates)

            f_O -= A_I @ g.ravel()

            u, info = gmres(A_O, f_O, f_O, tol=1E-10)
            u_fenics = helper.convert_nonlocal_function_into_fenics_function_dg(mesh, u, g, elementLabels)
        else:
            V = fenics.FunctionSpace(mesh, 'CG', 1)
            v = fenics.TestFunction(V)
            dx = fenics.Measure('dx', domain=mesh, subdomain_data=subdomains)
            l = self.source[0] * v * dx(1) + self.source[1] * v * dx(2) + 0 * v * dx(3)
            b = fenics.PETScVector()
            fenics.assemble(l, tensor=b)
            b_vec = b.get_local()

            if c_per:
                u = fenics.TrialFunction(V)
                a_per = fenics.dot(fenics.grad(u), fenics.grad(v)) * dx
                A_per = fenics.PETScMatrix()
                fenics.assemble(a_per, tensor=A_per)
                A_per = helper.petsc_to_csr_scipy_matrix(A_per)
                A = A + A_per

            A_O = A[vertexLabels > 0][:, vertexLabels > 0]
            A_I = A[vertexLabels > 0][:, vertexLabels < 0]
            f_O = b_vec[vertexLabels > 0]

            def g2(x):
                return 0.0

            g = np.apply_along_axis(g2, 1, vertices[vertexLabels < 0])
            f_O -= A_I @ g.ravel()
            # u = spsolve(A_O, f_O)
            u, info = gmres(A_O, f_O, f_O, tol=1E-10)
            if info != 0:
                print('Computing of state equation failed with code ' + str(info))
                raise RuntimeError
            u_fenics = helper.convert_nonlocal_function_into_fenics_function_cg(mesh, u, g, vertexLabels)
        return u_fenics

    def solve_adjoint(self, mesh, subdomains, u, u_bar, c_per):
        elements, vertices, elementLabels = mesh_data.convert_mesh(mesh, subdomains)
        empty_dict = {}
        self.nlfem_conf['is_ShapeDerivative'] = 0
        mesh_nl, A = nlfem.stiffnessMatrix_fromArray(elements, elementLabels, vertices, self.kernel, self.nlfem_conf,
                                                     empty_dict)
        vertexLabels = mesh_nl["vertexLabels"]

        if self.ansatz == "DG":
            V = fenics.FunctionSpace(mesh, 'DG', 1)
            p = fenics.TestFunction(V)
            dx = fenics.Measure('dx', domain=mesh)
            l = -(u - u_bar) * p * dx
            b = fenics.PETScVector()
            b = fenics.assemble(l, tensor=b)
            b_vec = b.get_local()

            nodeLabels = np.repeat(elementLabels, 3)

            if c_per:
                v = fenics.TrialFunction(V)
                a_per = fenics.dot(fenics.grad(v), fenics.grad(p)) * dx
                A_per = fenics.PETScMatrix()
                fenics.assemble(a_per, tensor=A_per)
                A_per = helper.petsc_to_csr_scipy_matrix(A_per)
                A = A + A_per

            A_O = A[nodeLabels > 0][:, nodeLabels > 0]
            A_I = A[nodeLabels > 0][:, nodeLabels < 0]
            f_O = b_vec[nodeLabels > 0]

            def g2(x):
                return 0.0

            num_boundary_elements = np.sum(elementLabels < 0)
            g = np.zeros(3 * num_boundary_elements)

            for index, element in enumerate(elements[elementLabels < 0]):
                for vertex_index in range(3):
                    vertex_coordinates = vertices[element[vertex_index]]
                    g[3 * index + vertex_index] = g2(vertex_coordinates)

            f_O -= A_I @ g.ravel()
            v, info = gmres(A_O.transpose(), f_O, f_O, tol=1E-10)
            v_fenics = helper.convert_nonlocal_function_into_fenics_function_dg(mesh, v, g, elementLabels)
        else:
            V = fenics.FunctionSpace(mesh, 'CG', 1)
            p = fenics.TestFunction(V)
            dx = fenics.Measure('dx', domain=mesh)
            l = -(u - u_bar) * p * dx
            b = fenics.PETScVector()
            b = fenics.assemble(l, tensor=b)
            b_vec = b.get_local()

            if c_per:
                v = fenics.TrialFunction(V)
                a_per = fenics.dot(fenics.grad(v), fenics.grad(p)) * dx
                A_per = fenics.PETScMatrix()
                fenics.assemble(a_per, tensor=A_per)
                A_per = helper.petsc_to_csr_scipy_matrix(A_per)
                A = A + A_per

            A_O = A[vertexLabels > 0][:, vertexLabels > 0]
            A_I = A[vertexLabels > 0][:, vertexLabels < 0]
            f_O = b_vec[vertexLabels > 0]

            def g2(x):
                return 0.0

            g = np.apply_along_axis(g2, 1, vertices[vertexLabels < 0])
            f_O -= A_I @ g.ravel()
            v, info = gmres(A_O.transpose(), f_O, f_O, tol=1E-10)
            if info != 0:
                print('Computing of adjoint equation failed with code ' + str(info))
                raise RuntimeError
            v_fenics = helper.convert_nonlocal_function_into_fenics_function_cg(mesh, v, g, vertexLabels)

        return v_fenics

    def get_shape_gradient(self, mesh_info, u, u_bar, v, nu, mu, l, kernel_dependent):
        number_nodes = mesh_info.mesh.num_vertices()
        W = fenics.VectorFunctionSpace(mesh_info.mesh, 'CG', 1)
        V = fenics.TestFunction(W)
        shape_gradient = fenics.TrialFunction(W)

        dx = fenics.Measure('dx', domain=mesh_info.mesh, subdomain_data=mesh_info.subdomains)

        normal = fenics.FacetNormal(mesh_info.mesh)
        ds = fenics.Measure('dS', domain=mesh_info.mesh, subdomain_data=mesh_info.interface)

        A = self.get_linear_elasticity_matrix(l, mu, shape_gradient, V, dx)
        if kernel_dependent:
            b_vec = self.get_shape_derivative_vector_labeled(mesh_info, u, u_bar, v, V, normal, dx, ds, nu)
        else:
            b_vec = self.get_shape_derivative_vector(u, u_bar, v, V, normal, dx, ds, nu)

        for index in mesh_info.indices_nodes_not_on_interface:
            b_vec[index] = 0.0
            b_vec[index + number_nodes] = 0.0

        shape_derivative = fenics.Function(W)
        shape_derivative.vector()[:] = b_vec

        for index in mesh_info.boundary_vertices:
            A[index, :] = scs.eye(1, 2 * number_nodes, index).tocsr()
            A[:, index] = scs.eye(1, 2 * number_nodes, index).tocsr().transpose()

            A[index + number_nodes, :] = scs.eye(1, 2 * number_nodes, index).tocsr()
            A[:, index + number_nodes] = scs.eye(1, 2 * number_nodes, index).tocsr().transpose()

        shape_gradient, info = gmres(A, b_vec, tol=1E-10)
        if info != 0:
            print('Computing of shape gradient failed with code ' + str(info))
        shape_gradient_fenics = fenics.Function(W)
        shape_gradient_fenics.vector().set_local(shape_gradient)
        # print(fenics.norm(shape_gradient_fenics, "L2"))
        return shape_gradient_fenics, shape_derivative

    def get_linear_elasticity_matrix(self, l, mu, shape_gradient, V, dx):
        linear_elasticity = fenics.inner(l * fenics.tr(fenics.sym(fenics.grad(shape_gradient))) * fenics.Identity(2)
                                         + 2.0 * mu * fenics.sym(fenics.grad(shape_gradient)),
                                         fenics.sym(fenics.grad(V))) * dx
        A = fenics.PETScMatrix()
        fenics.assemble(linear_elasticity, tensor=A)
        return helper.petsc_to_csr_scipy_matrix(A)

    def get_shape_derivative_vector(self, u, u_bar, v, V, normal, dx, ds, nu):

        objective_functional_derivative = 0.5 * fenics.div(V) * (u - u_bar) ** 2 * dx \
                                          - (u - u_bar) * fenics.dot(fenics.grad(u_bar), V) * dx
        shape_derivative = objective_functional_derivative \
                           - self.source[0] * v * fenics.div(V) * dx(1) - self.source[1] * v * fenics.div(V) * dx(2) \
                           + nu * (fenics.div(V("+")) - fenics.dot(fenics.dot(fenics.grad(V("+")), normal("+")),
                                                                   normal("+"))) * ds(12)
        b = fenics.PETScVector()
        fenics.assemble(shape_derivative, tensor=b)
        return b.get_local()

    def get_shape_derivative_vector_labeled(self, mesh_info, u, u_bar, v, V, normal, dx, ds, nu):

        objective_functional_derivative = 0.5 * fenics.div(V) * (u - u_bar) ** 2 * dx \
                                          - (u - u_bar) * fenics.dot(fenics.grad(u_bar), V) * dx
        shape_derivative = objective_functional_derivative \
                           - self.source[0] * v * fenics.div(V) * dx(1) - self.source[1] * v * fenics.div(V) * dx(2) \
                           + nu * (fenics.div(V("+")) - fenics.dot(fenics.dot(fenics.grad(V("+")), normal("+")),
                                                                   normal("+"))) * ds(12)

        if self.c_per:
            shape_derivative += self.c_per * (fenics.inner(fenics.grad(u), fenics.grad(v)) * fenics.div(V)
                                              - 2 * fenics.inner(fenics.grad(u),
                                                                 fenics.dot(fenics.sym(fenics.grad(V)),
                                                                            fenics.grad(v)))) * dx

        start_time = time.time()

        elements, vertices, elementLabels = mesh_data.convert_mesh(mesh_info.mesh, mesh_info.subdomains)
        # self.nlfem_conf['is_ShapeDerivative'] = 1
        shape_derivative_data = {'state': u.vector().get_local(),
                                 'adjoint': v.vector().get_local()}

        mesh_nl, b_nl = nlfem.stiffnessMatrix_fromArray(elements, elementLabels, vertices, self.kernel,
                                                        self.nlfem_shape_conf, shape_derivative_data)

        # Convert nlfem's shape derivative matrix into numpy array
        b_nl = b_nl.toarray().flatten()

        np.set_printoptions(threshold=np.inf)
        print('Assembling of nonlocal shape derivative part: ' + str(time.time() - start_time) + 's')
        b = fenics.PETScVector()
        fenics.assemble(shape_derivative, tensor=b)
        return b_nl + b.get_local()

    def get_mu(self, mesh, interface_vertices, indices_boundary_vertices, mu_min, mu_max):
        number_nodes = mesh.num_vertices()
        V = fenics.FunctionSpace(mesh, "CG", 1)
        v = fenics.TestFunction(V)
        mu = fenics.TrialFunction(V)
        dx = fenics.Measure("dx", domain=mesh)
        f = fenics.Expression("0.0", degree=1)
        lhs = fenics.inner(fenics.grad(mu), fenics.grad(v)) * dx
        rhs = f * v * dx
        A = fenics.PETScMatrix()
        b = fenics.PETScVector()
        fenics.assemble(lhs, tensor=A)
        fenics.assemble(rhs, tensor=b)
        A = helper.petsc_to_csr_scipy_matrix(A)
        b_vec = b.get_local()
        for index in indices_boundary_vertices:
            b_vec[index] = mu_min
            A[index, :] = scs.eye(1, number_nodes, index).tocsr()

        for index in interface_vertices:
            b_vec[index] = mu_max
            A[index, :] = scs.eye(1, number_nodes, index).tocsr()

        mu, info = gmres(A, b_vec, tol=1E-10)
        mu_fenics = fenics.Function(V)
        mu_fenics.vector()[:] = mu[:]
        return mu_fenics

    def line_search_update(self, mesh, subdomains, interface, indices_boundary_vertices, u_bar, direction,
                           shape_derivative, armijo_const, objective_function_value, alpha, nu, c_per):
        # compute sufficient decrease parameter
        sufficient_decrease = armijo_const * np.dot(shape_derivative.vector(), direction.vector())

        W = fenics.VectorFunctionSpace(mesh, "CG", 1)
        deformation = fenics.Function(W)
        c = -1.0
        counter = 1
        mesh_not_admissible = 0
        try:
            print("line_search... ")
            new_value, mesh_not_admissible = self.get_new_target_function_value(mesh, subdomains, interface,
                                                                                indices_boundary_vertices,
                                                                                direction, u_bar, nu, c_per, c)
        except RuntimeError as error:
            # print(error)
            new_value = objective_function_value + 1

        while new_value >= (objective_function_value + c * sufficient_decrease) or mesh_not_admissible:
            c = alpha * c
            try:
                new_value, mesh_not_admissible = self.get_new_target_function_value(mesh, subdomains, interface,
                                                                                    indices_boundary_vertices,
                                                                                    direction, u_bar, nu, c_per, c)
            except RuntimeError as error:
                # print(error)
                new_value = objective_function_value + 1

            counter += 1

        print("Number of line-search iterations: " + str(counter))
        deformation.vector()[:] = c * direction.vector()[:]
        return deformation, counter

    def get_new_target_function_value(self, mesh, subdomains, interface, indices_boundary_vertices, shape_gradient,
                                      u_bar, nu, c_per, c):
        local_mesh = fenics.Mesh(mesh)
        local_subdomains = fenics.cpp.mesh.MeshFunctionSizet(local_mesh, subdomains.dim())
        local_subdomains.set_values(subdomains.array())
        local_interface = fenics.cpp.mesh.MeshFunctionSizet(local_mesh, interface.dim())
        local_interface.set_values(interface.array())

        W = fenics.VectorFunctionSpace(local_mesh, "CG", 1)
        deformation = fenics.Function(W)
        deformation.vector()[:] = c * shape_gradient.vector()[:]

        fenics.ALE.move(local_mesh, deformation)

        info = helper.mesh_not_admissible(local_mesh.coordinates(), indices_boundary_vertices, 1E-8)

        V_new = fenics.FunctionSpace(local_mesh, self.ansatz, 1)
        u_bar_new = fenics.project(u_bar, V_new)
        u_new = self.solve_state(local_mesh, local_subdomains, c_per)
        value = 1. / 2. * fenics.norm(fenics.project(u_new - u_bar_new, V_new), "L2", local_mesh) ** 2
        if nu:
            ones = fenics.Function(V_new)
            ones.vector()[:] = 1.0
            ds = fenics.Measure("dS", domain=local_mesh, subdomain_data=local_interface)
            regularization_value = nu * fenics.assemble(ones('+') * ds(12))
            value += regularization_value
        return value, info

    def get_target_function_value(self, mesh, subdomains, interface, u_bar, nu, c_per):
        V_new = fenics.FunctionSpace(mesh, self.ansatz, 1)
        u_new = self.solve_state(mesh, subdomains, c_per)
        value = 1. / 2. * fenics.norm(fenics.project(u_new - u_bar, V_new), "L2", mesh) ** 2
        if nu:
            ones = fenics.Function(V_new)
            ones.vector()[:] = 1.0
            ds = fenics.Measure("dS", domain=mesh, subdomain_data=interface)
            regularization_value = nu * fenics.assemble(ones('+') * ds(12))
            value += regularization_value
        return value

    def solve_shape_problem(self, num_iterations):
        start_time = time.time()
        if self.l_bfgs:
            last_deformation = np.zeros(2 * self.current_mesh_data.mesh.num_vertices())

        k = 1
        while k <= num_iterations:
            iteration_time = time.time()
            print("Iteration: " + str(k))
            if self.remesh and k > 1 and (k - 1) % 5 == 0 and (k - 1) <= 10:
                self.file_name_current = mesh_data.remesh(self.current_mesh_data.mesh,
                                                          self.current_mesh_data.interface_vertices,
                                                          self.file_name, self.remesh_counter, self.element_size)
                self.remesh_counter += 1
                self.current_mesh_data = mesh_data.MeshData(self.file_name_current, self.interface_label,
                                                            self.boundary_label)

                if self.save_results:
                    self.results.save_interface_and_subdomains(self.current_mesh_data.interface,
                                                               self.current_mesh_data.subdomains)
                if self.l_bfgs:
                    num_vertices = self.current_mesh_data.mesh.num_vertices()
                    last_deformation = np.zeros(2 * num_vertices)
                    self.l_bfgs_algorithm = bfgs.BFGS(num_vertices, self.memory_length)
                print("Domain was remeshed")

            V_current = fenics.FunctionSpace(self.current_mesh_data.mesh, self.ansatz, 1)
            u_bar = fenics.project(self.u_bar, V_current)

            u = self.solve_state(self.current_mesh_data.mesh, self.current_mesh_data.subdomains, self.c_per)

            objective_function_value = self.get_target_function_value(self.current_mesh_data.mesh,
                                                                      self.current_mesh_data.subdomains,
                                                                      self.current_mesh_data.interface,
                                                                      u_bar, self.nu, self.c_per)

            v = self.solve_adjoint(self.current_mesh_data.mesh, self.current_mesh_data.subdomains, u, u_bar, self.c_per)

            mu = self.get_mu(self.current_mesh_data.mesh, self.current_mesh_data.interface_vertices,
                             self.current_mesh_data.boundary_vertices, self.mu_min, self.mu_max)

            shape_gradient, shape_derivative = self.get_shape_gradient(self.current_mesh_data, u, u_bar, v,
                                                                       self.nu, mu, self.lbd, self.kernel_dependent)

            if self.l_bfgs:
                bfgs_update = self.l_bfgs_algorithm.get_bfgs_update(self.current_mesh_data.mesh,
                                                                    last_deformation, shape_gradient, mu)
                deformation, counter = self.line_search_update(
                    self.current_mesh_data.mesh,
                    self.current_mesh_data.subdomains,
                    self.current_mesh_data.interface,
                    self.current_mesh_data.boundary_vertices,
                    self.u_bar, bfgs_update, shape_derivative, self.armijo_const,
                    objective_function_value,
                    self.alpha, self.nu, self.c_per)
                last_deformation[:] = deformation.vector()
            else:
                deformation, counter = self.line_search_update(
                    self.current_mesh_data.mesh,
                    self.current_mesh_data.subdomains,
                    self.current_mesh_data.interface,
                    self.current_mesh_data.boundary_vertices,
                    self.u_bar, shape_gradient, shape_derivative, self.armijo_const,
                    objective_function_value,
                    self.alpha, self.nu, self.c_per)

            if self.save_results:
                self.results.save_u(u)
                self.results.save_v(v)
                self.results.store_additional_results(self.current_mesh_data.mesh, shape_gradient, deformation,
                                                      objective_function_value, shape_derivative)

            fenics.ALE.move(self.current_mesh_data.mesh, deformation)

            if self.adapt_mu_max:
                if counter >= self.lame_threshold_up:
                    self.mu_max *= self.adapt_up
                elif counter <= self.lame_threshold_down:
                    self.mu_max *= self.adapt_down

            if self.save_results:
                self.results.save_interface_and_subdomains(self.current_mesh_data.interface,
                                                           self.current_mesh_data.subdomains)

            print('Iteration ' + str(k) + ' needed ' + str(time.time() - iteration_time) + 's')
            k = k + 1
        end_time = time.time() - start_time
        print('Algorithm ended after ' + str(end_time) + 's')
        if self.save_results:
            V_current = fenics.FunctionSpace(self.current_mesh_data.mesh, self.ansatz, 1)
            u_bar = fenics.project(self.u_bar, V_current)
            u = self.solve_state(self.current_mesh_data.mesh, self.current_mesh_data.subdomains, self.c_per)
            objective_function_value = self.get_target_function_value(self.current_mesh_data.mesh,
                                                                      self.current_mesh_data.subdomains,
                                                                      self.current_mesh_data.interface,
                                                                      u_bar, self.nu, self.c_per)
            self.results.save_objective_function_value(objective_function_value)
            self.results.save_u(u)
            self.results.plot_and_save_additional_results(end_time)

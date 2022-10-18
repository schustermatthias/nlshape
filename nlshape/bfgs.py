import numpy as np
import fenics


class BFGS:
    def __init__(self, num_vertices, memory_length):
        self.gradients = np.zeros((memory_length + 1, 2 * num_vertices))  # maybe without + 1
        self.deformations = np.zeros((memory_length, 2 * num_vertices))
        self.step_number = 0
        self.curvature_condition_was_violated = 0
        self.deformation_sum = np.zeros(2 * num_vertices)
        self.num_vertices = num_vertices
        self.memory_length = memory_length

    def restart(self):
        self.gradients = np.zeros((self.memory_length + 1, 2 * self.num_vertices))
        self.deformations = np.zeros((self.memory_length, 2 * self.num_vertices))
        self.step_number = 0
        self.curvature_condition_was_violated = 0
        self.deformation_sum = np.zeros(2 * self.num_vertices)

    def bfgs_step(self, mesh, shape_gradient, mu):
        W = fenics.VectorFunctionSpace(mesh, "CG", 1)
        q = fenics.Function(W)
        q.vector()[:] = shape_gradient.vector()

        alpha = np.zeros(self.memory_length)
        rho = np.zeros(self.memory_length)

        difference_gradients = fenics.Function(W)
        latest_difference_gradients = fenics.Function(W)
        latest_deformation = fenics.Function(W)

        min_index = min(self.step_number, self.memory_length)
        for index in range(min_index)[::-1]:
            difference_gradients.vector()[:] = self.initialize_gradient(mesh, index + 1).vector() \
                                               - self.initialize_gradient(mesh, index).vector()
            deformation = self.initialize_deformation(mesh, index)
            rho[index] = 1.0 / self.linear_elasticity_value(mesh, difference_gradients, deformation, mu)
            alpha[index] = self.linear_elasticity_value(mesh, deformation, q, mu) * rho[index]
            q.vector()[:] = q.vector() - alpha[index] * difference_gradients.vector()

        latest_difference_gradients.vector()[:] = self.initialize_gradient(mesh, min_index).vector() \
                                                  - self.initialize_gradient(mesh, min_index - 1).vector()
        latest_deformation.vector()[:] = self.initialize_deformation(mesh, min_index - 1).vector()
        gamma = self.linear_elasticity_value(mesh, latest_difference_gradients, latest_deformation, mu) \
                / self.linear_elasticity_value(mesh, latest_difference_gradients, latest_difference_gradients, mu)
        q.vector()[:] = gamma * q.vector()

        for index in range(min_index):
            difference_gradients.vector()[:] = self.initialize_gradient(mesh, index + 1).vector() \
                                               - self.initialize_gradient(mesh, index).vector()
            beta = self.linear_elasticity_value(mesh, difference_gradients, q, mu) * rho[index]
            q.vector()[:] = q.vector() + (alpha[index] - beta) * self.initialize_deformation(mesh, index).vector()
        return q

    def initialize_gradient(self, mesh, index):
        W = fenics.VectorFunctionSpace(mesh, "CG", 1)
        gradient = fenics.Function(W)
        gradient.vector()[:] = self.gradients[index]
        return gradient

    def initialize_deformation(self, mesh, index):
        W = fenics.VectorFunctionSpace(mesh, "CG", 1)
        mesh_coordinates = fenics.Function(W)
        mesh_coordinates.vector()[:] = self.deformations[index]
        return mesh_coordinates

    def linear_elasticity_value(self, mesh, U, V, mu):
        dx = fenics.Measure("dx", domain=mesh)
        a = fenics.inner(2.0 * mu * fenics.sym(fenics.grad(U)), fenics.sym(fenics.grad(V))) * dx
        return fenics.assemble(a)

    def update_gradient(self, gradient):
        if self.step_number <= self.memory_length:
            self.gradients[self.step_number, ] = gradient[:]
        else:
            self.gradients[0:self.memory_length, ] = self.gradients[1:self.memory_length + 1, ]

            self.gradients[self.memory_length, ] = gradient[:]

    def update_deformation(self, deformation):
        if self.step_number <= self.memory_length:
            self.deformations[self.step_number - 1, ] = deformation[:]
        else:
            self.deformations[0:self.memory_length - 1, ] = self.deformations[1:self.memory_length, ]

            self.deformations[self.memory_length - 1, ] = deformation[:]

    def get_bfgs_update(self, mesh, last_deformation, shape_gradient, mu):
        if self.step_number == 0:
            self.update_gradient(shape_gradient.vector())
            self.step_number += 1
            return shape_gradient
        if self.step_number > 0:
            # Check, if curvature condition holds.
            min_index = min(self.memory_length, self.step_number)
            W = fenics.VectorFunctionSpace(mesh, "CG", 1)

            new_deformation = fenics.Function(W)
            if self.curvature_condition_was_violated:
                new_deformation.vector()[:] = last_deformation + self.deformation_sum
            else:
                new_deformation.vector()[:] = last_deformation

            new_difference_gradients = fenics.Function(W)
            new_difference_gradients.vector()[:] = shape_gradient.vector() \
                                                   - self.initialize_gradient(mesh, min_index - 1).vector()
            value = self.linear_elasticity_value(mesh, new_difference_gradients, new_deformation, mu)
            if value > 0 and np.allclose(value, 0.) is False:
                self.update_gradient(shape_gradient.vector())
                self.update_deformation(new_deformation.vector())

                self.curvature_condition_was_violated = 0
                self.deformation_sum = np.zeros(2 * mesh.num_vertices())

                bfgs_update = self.bfgs_step(mesh, shape_gradient, mu)
                self.step_number += 1
                return bfgs_update
            else:
                print("curvature condition violated")
                self.curvature_condition_was_violated += 1
                if self.curvature_condition_was_violated == 1:
                    self.deformation_sum[:] = new_deformation.vector()
                else:  # Restart bfgs, if curvature condition is violated the second time in a row
                    print("restart bfgs")
                    self.restart()
                    self.update_gradient(shape_gradient.vector())
                    self.step_number += 1
                return shape_gradient

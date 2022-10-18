import matplotlib.pyplot as plt
import datetime
import os

import fenics


def save_dict(dictionary, name):
    file = open(name + ".txt", "w")
    for key in dictionary.keys():
        file.write(str(key) + ": " + str(dictionary[key]) + "\n")
    file.close()


class Results:
    def __init__(self, conf, interface, subdomains, target_interface):
        timestamp = datetime.datetime.now().strftime("%Y_%m_%d") + "/" + datetime.datetime.now().strftime("%H_%M")
        self.results_folder = "results/" + timestamp + "/"

        self.file_u = fenics.File(self.results_folder + "u.pvd")
        self.file_v = fenics.File(self.results_folder + "v.pvd")

        self.file_target_interface = fenics.File(self.results_folder + "target_interface.pvd")
        self.file_interface = fenics.File(self.results_folder + "interface.pvd")
        # self.file_shape_gradient = fenics.File(self.results_folder + "shape_gradient.pvd")
        # self.file_shape_derivative = fenics.File(self.results_folder + "shape_derivative.pvd")
        self.file_subdomains = fenics.File(self.results_folder + "subdomains.pvd")

        interface.rename("interface", "")
        subdomains.rename("subdomains", "")
        target_interface.rename("target_interface", "")

        self.file_interface << interface
        self.file_subdomains << subdomains
        self.file_target_interface << target_interface

        save_dict(conf, self.results_folder + "configuration")

        self.shape_gradient_norm_history = []
        self.shape_derivative_norm_history = []
        self.deformation_norm_history = []
        self.objective_function_value_history = []

        self.counter = 1

    def save_target_data(self, u_bar, subdomains):
        file_u_bar = fenics.File(self.results_folder + "u_bar.pvd")
        file_u_bar << u_bar
        file_subdomains = fenics.File(self.results_folder + "target_subdomains.pvd")
        file_subdomains << subdomains

    def save_u(self, u):
        u.rename("u", "")
        self.file_u << u

    def save_v(self, v):
        v.rename("v", "")
        self.file_v << v

    def store_additional_results(self, mesh, shape_gradient, deformation, objective_function_value,
                                 shape_derivative):
        shape_gradient_norm = fenics.norm(shape_gradient, "L2", mesh)
        shape_derivative_norm = fenics.norm(shape_derivative, "L2", mesh)
        deformation_norm = fenics.norm(deformation, "L2", mesh)
        self.shape_gradient_norm_history.append(shape_gradient_norm)
        self.shape_derivative_norm_history.append(shape_derivative_norm)
        self.deformation_norm_history.append(deformation_norm)
        self.objective_function_value_history.append(objective_function_value)

        fenics.plot(deformation)
        plt.savefig(self.results_folder + 'deformation_' + str(self.counter))
        plt.clf()
        plt.close()

        # fenics.plot(shape_gradient)
        # plt.savefig('shape_gradient_' + str(self.counter))
        # plt.clf()
        # plt.close()
        self.counter += 1

    def save_objective_function_value(self, value):
        self.objective_function_value_history.append(value)

    def save_interface_and_subdomains(self, interface, subdomains):
        interface.rename("interface", "")
        subdomains.rename("subdomains", "")

        self.file_interface << interface
        self.file_subdomains << subdomains

    def plot_and_save_additional_results(self, time):
        fig = []
        axs = []
        for index in range(4):
            fig.append(plt.figure())
            axs.append(fig[index].add_subplot(111))
        axs[0].plot(self.shape_gradient_norm_history, label="shape gradient norm")
        axs[1].plot(self.shape_derivative_norm_history, label="shape derivative norm")
        axs[2].plot(self.deformation_norm_history, label="deformation norm")
        axs[3].plot(self.objective_function_value_history, label="objective function value")

        os.system("mkdir " + self.results_folder + "additional_results")
        names = ["shape_gradient_norm", "shape_derivative_norm", "deformation_norm", "objective_function_value"]
        for index in range(4):
            axs[index].legend()
            fig[index].savefig(self.results_folder + "additional_results/" + names[index] + ".png")
        # plt.show()

        # Write a csv file
        # file_additional_results = open(self.results_folder + "additional_results.csv", "w")
        # file_additional_results.write("objective_function_value" + "; " + "shape_gradient_norm" + "; "
        #                              + "deformation_norm" + "; " + "step_size" + "; " + "mu_max" + "\n")
        # for k in range(self.counter-1):
        #    file_additional_results.write(str(self.objective_function_value_history[k]) + "; "
        #                                  + str(self.shape_gradient_norm_history[k]) + "; "
        #                                  + str(self.deformation_norm_history[k]) + "; "
        #                                  + str(self.step_size_history[k]) + "; " + str(self.mu_max_history[k]) + "\n")

        # file_additional_results.write(str(self.objective_function_value_history[-1]))

        file_additional_results = open(self.results_folder + "additional_results.txt", "a")
        file_additional_results.write("shape gradient norm: " + str(self.shape_gradient_norm_history) + "\n")
        file_additional_results.write("shape derivative norm: " + str(self.shape_derivative_norm_history) + "\n")
        file_additional_results.write("deformation norm: " + str(self.deformation_norm_history) + "\n")
        file_additional_results.write("objective function value: " + str(self.objective_function_value_history) + "\n")
        file_additional_results.write("Algorithm ended after " + str(time) + "s")

        file_additional_results.close()

import warnings
import mesh_data
from ex2_sym import configuration
import nonlocal_shape_problem

warnings.filterwarnings("ignore", message="Changing the sparsity structure of a csr_matrix is expensive. "
                                          "lil_matrix is more efficient.")

if __name__ == '__main__':
    #  mesh_data.convert_msh_to_xdmf("example_1")
    conf = configuration
    interface_problem = nonlocal_shape_problem.NonlocalShapeProblem(conf)
    interface_problem.solve_shape_problem(conf["number_iterations"])

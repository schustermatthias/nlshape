import numpy as np
import os

import fenics
import meshio


def get_mesh(name):
    mesh = fenics.Mesh()
    with fenics.XDMFFile("mesh/" + name + ".xdmf") as infile:
        infile.read(mesh)
    return mesh


def get_mesh_and_mesh_data(name, interface_label):
    mesh = fenics.Mesh()
    with fenics.XDMFFile("mesh/" + name + ".xdmf") as infile:
        infile.read(mesh)

    subdomain_data = fenics.MeshValueCollection("size_t", mesh, 2)
    with fenics.XDMFFile("mesh/" + name + "_subdomains.xdmf") as infile:
        infile.read(subdomain_data)
    subdomain_function = fenics.cpp.mesh.MeshFunctionSizet(mesh, subdomain_data)

    boundary_data = fenics.MeshValueCollection("size_t", mesh, 1)
    with fenics.XDMFFile("mesh/" + name + "_boundaries.xdmf") as infile:
        infile.read(boundary_data)
    interface_function = fenics.cpp.mesh.MeshFunctionSizet(mesh, boundary_data)
    interface_array = interface_function.array()
    for l in range(interface_array.size):
        if interface_array[l] != interface_label:
            interface_function.set_value(l, 0)

    return mesh, subdomain_function, interface_function


def get_interface_indices(mesh, interface, interface_label):
    # Find facets on interior boundary
    indices_interface_facets = []
    for facet in range(len(interface)):
        if interface[facet] == interface_label:
            indices_interface_facets.append(facet)

    # Find vertices on interior boundary
    interface_vertices = []
    interface_elements = []
    for cell in fenics.cells(mesh):
        for facet in fenics.facets(cell):
            if facet.index() in indices_interface_facets:
                interface_elements.append(cell.index())
                for vertex in fenics.vertices(facet):
                    interface_vertices.append(vertex.index())

    return list(set(interface_vertices)), list(set(interface_elements))


def convert_mesh(mesh, subdomains):
    elements = np.array(mesh.cells(), dtype=int)
    vertices = mesh.coordinates()
    num_elements = elements.shape[0]
    elementLabels = np.array(subdomains.array(), dtype='long')
    for triangle in range(num_elements):
        if elementLabels[triangle] == 3:
            elementLabels[triangle] = -1.0
    return elements, vertices, elementLabels


def get_indices_boundary_vertices(mesh, subdomains, boundary_label):
    elements = np.array(mesh.cells(), dtype=int)
    num_elements = elements.shape[0]
    elementLabels = np.array(subdomains.array(), dtype='long')
    indices_boundary_vertices = []
    for triangle in range(num_elements):
        if elementLabels[triangle] == boundary_label:
            for vertex in elements[triangle]:
                indices_boundary_vertices.append(vertex)
    return list(set(indices_boundary_vertices))


def convert_msh_to_xdmf(name):
    msh = meshio.read("mesh/" + name + ".msh")
    for cell in msh.cells:
        if cell.type == "triangle":
            triangle_cells = cell.data
        elif cell.type == "line":
            line_cells = cell.data

    for key in msh.cell_data_dict["gmsh:physical"].keys():
        if key == "triangle":
            triangle_data = msh.cell_data_dict["gmsh:physical"][key]
        elif key == "line":
            line_data = msh.cell_data_dict["gmsh:physical"][key]

    mesh = meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells})
    mesh_function_subdomains = meshio.Mesh(points=msh.points,
                                cells=[("triangle", triangle_cells)],
                                cell_data={"name_to_read": [triangle_data]})
    mesh_function_boundaries = meshio.Mesh(points=msh.points,
                                           cells=[("line", line_cells)],
                                           cell_data={"name_to_read": [line_data]})
    meshio.write("mesh/" + name + ".xdmf", mesh)
    meshio.write("mesh/" + name + "_subdomains.xdmf", mesh_function_subdomains)
    meshio.write("mesh/" + name + "_boundaries.xdmf", mesh_function_boundaries)
    os.system('meshio-convert mesh/' + name + ".xdmf " + "mesh/" + name + ".xdmf -p -z")


def remesh(mesh, interface_vertices, name, num, element_size):
    new_name = name + "_remeshed_" + str(num)
    create_geo_file(mesh, interface_vertices, new_name)
    os.system('gmsh mesh/' + new_name + '.geo -2 -clscale ' + str(element_size)
              + ' -format msh22 -o mesh/' + new_name + '.msh')
    convert_msh_to_xdmf(new_name)
    return new_name


def create_geo_file(mesh, interface_vertices, name):
    dummy_file = open("mesh/dummy_shape.geo", "r")
    data = dummy_file.readlines()

    new_file = open("mesh/" + name + ".geo", "w")

    for line in data[0:10]:
        new_file.write(line)

    # Prepare text, that is added to the .geo-file.
    vertices = mesh.coordinates()
    counter = 0
    interface = "Spline(10) = {"
    for index in interface_vertices:
        coordinates = vertices[index]
        new_file.write("Point(" + str(10 + counter) + ") = {" + str(coordinates[0]) + ", " + str(coordinates[1])
                        + ", 0, lc};\n")
        interface += str(10 + counter) + ", "
        counter += 1
    interface += "10};"

    for line in data[10:40]:
        new_file.write(line)

    new_file.write(interface)

    for line in data[40:]:
        new_file.write(line)

    dummy_file.close()
    new_file.close()


class MeshData:

    def __init__(self, file_name, interface_label, boundary_label):
        self.mesh, self.subdomains, self.interface = get_mesh_and_mesh_data(file_name, interface_label)
        self.interface_vertices, self.interface_elements = get_interface_indices(self.mesh, self.interface,
                                                                                 interface_label)
        self.boundary_vertices = get_indices_boundary_vertices(self.mesh, self.subdomains, boundary_label)
        vertices = np.array([i for i in range(self.mesh.num_vertices())])
        self.indices_nodes_not_on_interface = list(set(vertices) - set(self.interface_vertices))

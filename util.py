import numpy as np
from collections import OrderedDict as oD
import scipy.spatial as sp


def border_temperature_boundary_conditions(mesh):
    """
    Function that returns three vectors for the standard boundary condition for the Poisson temperature problem.
    :param mesh: Mesh object to be used to obtain the points information.
    """
    # Acquiring borders
    #     _d
    #  a |_| c
    #     b
    vertex_a = np.where(mesh.x == np.min(mesh.x))[0]
    vertex_b = np.where(mesh.y == np.min(mesh.y))[0]
    vertex_c = np.where(mesh.x == np.max(mesh.x))[0]
    vertex_d = np.where(mesh.y == np.max(mesh.y))[0]

    # Defining indices, types and values
    indices = list(oD.fromkeys(np.append(vertex_a, vertex_b)))
    values = np.zeros(len(indices))
    types = np.zeros(len(values)) + 1
    indices = list(oD.fromkeys(np.append(indices, list(oD.fromkeys(np.append(vertex_c, vertex_d))))))
    types = np.hstack((types, np.zeros(len(indices) - len(values))))
    values = np.append(values, np.zeros(len(indices) - len(values)) + 1)

    return indices, values, types


def hagen_poiseuille_boundary_conditions(mesh):
    """
    Function that returns three vectors for the standard boundary condition for the Poisson temperature problem.
    :param mesh: Mesh object to be used to obtain the points information.
    """
    # Acquiring borders
    #     _d
    #  a |_| c
    #     b
    vertex_a = np.where(mesh.x == np.min(mesh.x))[0]
    vertex_b = np.where(mesh.y == np.min(mesh.y))[0]
    vertex_c = np.where(mesh.x == np.max(mesh.x))[0]
    vertex_d = np.where(mesh.y == np.max(mesh.y))[0]

    # Defining psi
    indices = np.append(vertex_a, vertex_c)
    values = list(map(lambda x: mesh.y[x], indices))
    indices = np.array(list(oD.fromkeys(np.append(indices, vertex_b))))
    values = np.append(values, np.zeros(len(indices) - len(values)))
    indices = np.array(list(oD.fromkeys(np.append(indices, vertex_d))))
    values = np.append(values, np.zeros(len(indices) - len(values)) + 1.0)
    vector = [{"name": "psi", "indices": np.copy(indices), "values": np.copy(values), "type": True}]

    # Defining velocity (x axis component)
    indices = np.copy(vertex_a)
    values = np.zeros(len(indices)) + 1e-6
    indices = np.array(list(oD.fromkeys(np.append(indices, np.append(vertex_b, vertex_d)))))
    values = np.append(values, np.zeros(len(indices) - len(values)))
    vector.append({"name": "vel_x", "indices": np.copy(indices), "values": np.copy(values), "type": True})

    # Defining velocity (y axis component)
    vector.append({"name": "vel_y", "indices": np.copy(indices), "values": np.copy(values * 0.), "type": True})

    return vector


def check_method_call(*args):
    """
    Tests if the arguments are valid
    :param args: Any argument group that required for the method
    """
    for arg in args:
        if arg is not None:
            try:
                if len(arg):
                    return
            except TypeError:
                if arg:
                    return
    raise ValueError("Method called incorrectly, please read the documentation and try changing the arguments.")


def check_list_or_object(_list: list, _class: type):
    """
    Checks if _list parameter is a list of the same type as _class, or if it is a single value of that type.
    :param _list: List or single value to be checked.
    :param _class: Type of object desired.
    :return: List of values of the type _class.
    """
    if _list:
        if isinstance(_list, (list, tuple)):
            if all([isinstance(obj, _class) for obj in _list]):
                return _list
            else:
                raise TypeError("Object in argument list is not a {0}!".format(_class.__name__))
        else:
            if isinstance(_list, _class):
                return [_list]
            else:
                raise TypeError("Argument used is not a {0}!".format(_class.__name__))
    else:
        return []


def check_if_instance(obj, clazz):
    """
    Checks if the object is an instance of the class type.
    :param obj: Object instance.
    :param clazz: Class type, or list of types.
    """
    if not isinstance(obj, clazz):
        if isinstance(clazz, (list, tuple)):
            names = ", ".join([_clazz.__name__ for _clazz in clazz])
            raise TypeError("Incorrect method usage, parameter of type --{0}-- used,"
                            " please use a ({1}) object.".format(obj.__class__.__name__, names))
        raise TypeError("Incorrect method usage, parameter of type --{0}-- used,"
                        " please use a ({1}) object.".format(obj.__class__.__name__, clazz.__name__))


def style_plot(param_x, param_y):
    """
    Alter the plot styling.
    :param param_x: List of x coordinates.
    :param param_y: List of y coordinates.
    """
    import matplotlib.pyplot as plt

    def max_amplitude(_list):
        return max(_list) - min(_list)

    default_size = (6.4, 4.8)

    fig = plt.gcf()
    fig.set_size_inches((default_size[0] * max_amplitude(param_x) / max_amplitude(param_y), 4.8))
    fig.subplots_adjust(left=0.1 - 0.01 * max_amplitude(param_x) / max_amplitude(param_y), right=0.95)


def create_new_surface(*imported_points, lt_version=True):
    """
    Create new surface
    :return: Element information
    """

    x, y, ien, delauney_surfaces = 0, 0, 0, 0

    if imported_points:
        # Custom geometry.
        imported_points = np.array(imported_points[0])
        delauney_surfaces = sp.Delaunay(imported_points[:, :2])

        if lt_version:
            x = delauney_surfaces.points[:, 0]
            y = delauney_surfaces.points[:, 1]
            ien = delauney_surfaces.simplices
        else:
            import pygmsh
            geom = pygmsh.built_in.Geometry()
            for tri in delauney_surfaces.simplices:
                geom.add_polygon([[delauney_surfaces.points[tri[0]][0], delauney_surfaces.points[tri[0]][1], 0.0],
                                  [delauney_surfaces.points[tri[1]][0], delauney_surfaces.points[tri[1]][1], 0.0],
                                  [delauney_surfaces.points[tri[2]][0], delauney_surfaces.points[tri[2]][1], 0.0]])
            x, y, ien, delauney_surfaces = use_meshio(None, geom)

    else:
        if not lt_version:
            import pygmsh
            geom = pygmsh.built_in.Geometry()

            # Default surface.
            geom.add_polygon([
                              [0.0, 0.0, 0.0],
                              [1.0, 0.0, 0.0],
                              [0.5, 0.5, 0.0],
                              [1.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                             ],
                             lcar=0.05)

            x, y, ien, delauney_surfaces = use_meshio(None, geom)

    return x, y, ien, delauney_surfaces


def use_meshio(filename, geometry):
    """
    Use MeshIO library for creating or importing the mesh point structure from Gmsh.
    :param filename: Name of ".msh" file to be imported.
    :param geometry: A Geometry object from the PyGmsh library
    :return: x, y, ien - The mesh point structure.
    """
    import fnmatch
    import pygmsh
    import meshio

    # Saving mesh as .vtk exportable file
    # points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geometry)
    # meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data) # TODO: FIX LIBRARY

    if fnmatch.fnmatch(filename, "*.msh"):
        points = meshio.read(filename).points
    else:
        points = meshio.read(filename + ".msh").points

    x_ = points[:, 0]
    y_ = points[:, 1]

    delauney_surfaces = sp.Delaunay(points[:, :2])
    ien_ = delauney_surfaces.simplices

    return x_, y_, ien_, delauney_surfaces


def sparse_to_vector(vector):
    """
    Converts one dimensional sparse matrix to a vector array to allow more features.
    :param vector: Vector as a sparse matrix (x,1) or (1,x).
    :return: Vector as an one dimensional array (x,).
    """
    return np.ravel(vector.toarray())


def get_area(x_coord, y_coord):
    """
    Calculate area of a triangle, given it's vertices.
    :param x_coord: The x coordinates of the vertices, in order.
    :param y_coord: The y coordinates of the vertices, in order.
    :return: The area of the triangle.
    """
    check_method_call(x_coord, y_coord)
    return ((x_coord[0] * y_coord[1] - x_coord[1] * y_coord[0]) +
            (x_coord[1] * y_coord[2] - x_coord[2] * y_coord[1]) +
            (x_coord[2] * y_coord[0] - x_coord[0] * y_coord[2])) / 2.0


def get_dt(mesh):
    """
    Calculates optimal dt based on mesh.
    :param mesh: Mesh object that contains element and point informations.
    :return: Optimal dt.
    """
    import itertools as it

    def _h(elem):
        _a = ((mesh.x[elem[0]] * mesh.y[elem[1]] - mesh.x[elem[1]] * mesh.y[elem[0]]) +
              (mesh.x[elem[1]] * mesh.y[elem[2]] - mesh.x[elem[2]] * mesh.y[elem[1]]) +
              (mesh.x[elem[2]] * mesh.y[elem[0]] - mesh.x[elem[0]] * mesh.y[elem[2]])) / 2.0
        _aux = np.array(list(it.combinations(elem, 2)))
        _l = min(list(map(lambda p_1, p_2: np.sqrt((mesh.x[p_1] - mesh.x[p_2]) ** 2 + (mesh.y[p_1] - mesh.y[p_2]) ** 2),
                          _aux[:, 0], _aux[:, 1])))
        return _a / _l

    return min(list(map(lambda x: _h(x), mesh.ien)))

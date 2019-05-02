import numpy as np
from collections import OrderedDict as oD
import scipy.spatial as sp
import scipy.sparse
import pickle
import matplotlib.pyplot as plt
from cycler import cycler


def get_dict_value(values: list, default: any, expression_input: np.ndarray = None):
    """
    Returns first value in list that is not None. To prevent False values when facing zero.
    :param values: List of values.
    :param default: Default value if no other is found.
    :param expression_input: Input of the expression, if provided.
    :return: First not None value or default value.
    """
    for value in values:
        if value is not None:
            if callable(value):
                return value(expression_input) if expression_input is not None else default
            return value
    return default


def build_boundary_conditions(mesh, values_dict: dict = None, types_dict: dict = None, default_value: float = 0.0,
                              default_type: bool = False):
    """
    Function that returns three vectors for the standard boundary condition for the Poisson temperature problem.
    Function used to create the vectors for use as boundary conditions builders. Creates the information for points
    distributed along the borders of the mesh. Uses a dict with each side name as a key and its value or type as values.
    :param mesh: Mesh object to be used to obtain the points information.
    :param values_dict: Dictionary of values for each side of the contours.
                        Each key corresponds to a side, acceptable keys are (in order of priority):
                        "n", "north", "t", "top";
                        "e", "east", "r", "right";
                        "s", "south", "b", "bottom";
                        "w", "west", "l", "left".
                        If more than one key is used for a side the one with more priority is used.
                        If no keys are found for a side the default_value is used.
                        There is also the key "a", for a repeatable condition for all sides.
    :param types_dict: Dictionary of types for each side of the contours.
                       True for Dirichlet Condition,
                       False for Neumann Condition.
    :param default_value: Default value for unset boundaries, will be used for any boundary not defined in the dict.
                          Value defaults to 0.
    :param default_type: Default type for unset boundaries, will be used for any boundary not defined in the dict.
                          Type defaults to Neumann (False).
    :return: [indices, values, types]. Each represents the boundary conditions information vectors of each parameter.
             Note: The boundaries are defined in clockwise order, mutual points, such as origin, are defined in order:
             north, east, south and west.
    """
    # Acquiring borders
    #     _N_
    #  W |   |
    #    |___| E
    #      S

    if not values_dict:
        values_dict = {}
    if not types_dict:
        types_dict = {}
    if "all" in values_dict:
        values_dict["n"] = values_dict.get("all")
        values_dict["e"] = values_dict.get("all")
        values_dict["s"] = values_dict.get("all")
        values_dict["w"] = values_dict.get("all")
    if "all" in types_dict:
        types_dict["n"] = types_dict.get("all")
        types_dict["e"] = types_dict.get("all")
        types_dict["s"] = types_dict.get("all")
        types_dict["w"] = types_dict.get("all")

    # Finding indices of nodes in contours
    vertex_n = np.where(mesh.y == mesh.y.max())[0]
    vertex_e = np.where(mesh.x == mesh.x.max())[0]
    vertex_s = np.where(mesh.y == mesh.y.min())[0]
    vertex_w = np.where(mesh.x == mesh.x.min())[0]

    # Defining indices, types and values of each side, moving clockwise
    # Setting information of north side
    indices = vertex_n
    values = np.zeros(len(indices)) + get_dict_value([values_dict.get("n"), values_dict.get("north"),
                                                      values_dict.get("t"), values_dict.get("top")],
                                                     default_value, indices)
    types = np.zeros(len(indices)) + get_dict_value([types_dict.get("n"), types_dict.get("north"),
                                                     types_dict.get("t"), types_dict.get("top")],
                                                    default_type)

    # Setting information of east side
    indices = np.array(list(oD.fromkeys(np.append(indices, vertex_e))))
    values = np.append(values, np.zeros(len(indices) - len(values)) +
                       get_dict_value([values_dict.get("e"), values_dict.get("east"),
                                       values_dict.get("r"), values_dict.get("right")],
                                      default_value, indices[-(len(indices) - len(values))::]))
    types = np.append(types, np.zeros(len(indices) - len(types)) +
                      get_dict_value([types_dict.get("e"), types_dict.get("east"),
                                      types_dict.get("r"), types_dict.get("right")],
                                     default_type))

    # Setting information of south side
    indices = np.array(list(oD.fromkeys(np.append(indices, vertex_s))))
    values = np.append(values, np.zeros(len(indices) - len(values)) +
                       get_dict_value([values_dict.get("s"), values_dict.get("south"),
                                       values_dict.get("b"), values_dict.get("bottom")],
                                      default_value, indices[-(len(indices) - len(values))::]))
    types = np.append(types, np.zeros(len(indices) - len(types)) +
                      get_dict_value([types_dict.get("s"), types_dict.get("south"),
                                      types_dict.get("b"), types_dict.get("bottom")],
                                     default_type))

    # Setting information of west side
    indices = np.array(list(oD.fromkeys(np.append(indices, vertex_w))))
    values = np.append(values, np.zeros(len(indices) - len(values)) +
                       get_dict_value([values_dict.get("w"), values_dict.get("west"),
                                       values_dict.get("l"), values_dict.get("left")],
                                      default_value, indices[-(len(indices) - len(values))::]))
    types = np.append(types, np.zeros(len(indices) - len(types)) +
                      get_dict_value([types_dict.get("w"), types_dict.get("west"),
                                      types_dict.get("l"), types_dict.get("left")],
                                     default_type))

    return indices, values, types


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


def style_plot(param_x: list, param_y: list):
    """
    Alter the plot styling.
    :param param_x: List of x coordinates.
    :param param_y: List of y coordinates.
    """
    def max_amplitude(_list):
        return max(_list) - min(_list)

    default_size = (6.4, 4.8)

    fig = plt.gcf()
    fig.set_size_inches((default_size[0] * max_amplitude(param_x) / max_amplitude(param_y), 4.8))
    fig.subplots_adjust(left=0.1 - 0.01 * max_amplitude(param_x) / max_amplitude(param_y), right=0.95)


def show_comparison(x_coordinates: np.ndarray, analytic_expression: callable, numeric_solution: [dict, np.ndarray],
                    numeric_label: str = "Numeric Solution", analytic_label: str = "Analytic Solution",
                    title: str = None, x_label: str = None, y_label: str = None, save_file_as: str = None):
    """
    Method that shows the comparison between the analytic and numeric solutions.
    :param x_coordinates: Array of input values for function.
    :param numeric_solution: Array of values for the numeric solution.
    :param analytic_expression: Function that describes the analytic solution.
    :param numeric_label: Label for numeric solution on graph.
    :param analytic_label: Label for analytic solution on graph.
    :param title: Title of plot figure.
    :param x_label: Label for the x axis.
    :param y_label: Label for the y axis.
    :param save_file_as: Filename used to save generated figure. If not defined figure is not saved.
    :return: Displays the graphical comparison.
    """
    check_method_call(x_coordinates)
    check_method_call(analytic_expression)
    check_method_call(numeric_solution)

    analytic_solution = analytic_expression(x_coordinates)

    default_cycler = cycler('color', ['b', 'g', 'k']) * cycler('linestyle', ['--', '-.', ':'])
    plt.rc('axes', prop_cycle=default_cycler)

    plt.plot(x_coordinates, analytic_solution, "r-", label=analytic_label)
    if isinstance(numeric_solution, dict):
        [plt.plot(x_coordinates, numeric_solution[key], label=("{:.4f}s".format(key) if isinstance(key, (float, int))
                                                               else key)) for key in sorted(numeric_solution)]
    else:
        plt.plot(x_coordinates, numeric_solution, "b--", label=numeric_label)

    axes = plt.gca()
    if x_label:
        axes.set_xlabel(x_label)
    if y_label:
        axes.set_ylabel(y_label)
    if title:
        axes.set_title(title)

    plt.grid()
    plt.legend()

    # Calculate errors
    numeric_solution = np.array(numeric_solution if not isinstance(numeric_solution, dict) else
                                numeric_solution[max(numeric_solution.keys())])

    error_array = np.nan_to_num(abs(numeric_solution - analytic_solution)/analytic_solution)
    print("Mean Error: {0}\nStandard Error: {1}".format(np.mean(error_array), np.std(error_array)))

    if save_file_as is not None and isinstance(save_file_as, str):
        plt.savefig("{0}".format(save_file_as))

    return plt.show()


def create_new_surface(imported_points: list = None, lt_version: bool = True):
    """
    Create new surface
    :param imported_points: List of tuples containing each point coordinate.
                            Eg: [(x_1, y_1), (x_2, y_2), ...]
    :param lt_version: Light version does not utilise meshio library.
    :return: Element information
    """
    x, y, ien, delauney_surfaces = 0, 0, 0, 0

    if imported_points:
        # Custom geometry.
        imported_points = np.array(imported_points)
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
            x, y, ien, delauney_surfaces = use_meshio("", geom)

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

            x, y, ien, delauney_surfaces = use_meshio("", geom)

    return x, y, ien, delauney_surfaces


def use_meshio(filename: str, geometry):
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

    file_data = meshio.read(filename + ("" if fnmatch.fnmatch(filename, "*.msh") else ".msh"))

    x_ = file_data.points[:, 0]
    y_ = file_data.points[:, 1]
    ien_ = file_data.cells["triangle"]

    delauney_surfaces = sp.Delaunay(file_data.points[:, :2])
    # TODO: REFINE BY SPLITTING IS STILL NOT WORKING
    delauney_surfaces.simplices = ien_.astype("int32")
    # ien_ = delauney_surfaces.simplices

    return x_, y_, ien_, delauney_surfaces


def sparse_to_vector(vector: scipy.sparse.spmatrix):
    """
    Converts one dimensional sparse matrix to a vector array to allow more features.
    :param vector: Vector as a sparse matrix (x,1) or (1,x).
    :return: Vector as an one dimensional array (x,).
    """
    return np.ravel(vector.toarray())


def get_area(x_coord: list, y_coord: list):
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


def get_side_length(x_coord: list, y_coord: list):
    """
    Calculate average side length of a triangle, given it's vertices.
    """
    return (np.sqrt((x_coord[0] - x_coord[1])**2 + (y_coord[0] - y_coord[1])**2) +
            np.sqrt((x_coord[1] - x_coord[2])**2 + (y_coord[1] - y_coord[2])**2) +
            np.sqrt((x_coord[2] - x_coord[0])**2 + (y_coord[2] - y_coord[0])**2)) / 3.0


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


def save(obj, filename: str = "temp"):
    """
    Saves a temporary Python obj to prevent repeated work.
    :param obj: Any object.
    :param filename: Optional file name.
    """
    with open(filename + ".dat", "wb") as file:
        pickle.dump(obj, file)


def load(filename: str = "temp"):
    """
    Loads a previously created Python obj to prevent repeated work.
    :param filename: Optional file name.
    :return: Original object.
    """
    with open(filename + ".dat", "rb") as file:
        return pickle.load(file)

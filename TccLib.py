# -*- coding: utf-8 -*-
# -------------------------------------- by: LUCAS CARVALHO DE SOUSA ---------------------------------------------------
# Esta biblioteca foi criada para o a solução de sistemas diferenciais através do Método de Elementos Finitos
# This library was created for educational purposes, it is meant to be used for the solution of differential equation
# systems
__author__ = "Lucas Carvalho de Sousa"
__author_email__ = "encarvlucas@gmail.com"
__website__ = "https://github.com/encarvlucas/LucasCarvalhoTCC"


# -- Functions ---------------------------------------------------------------------------------------------------------
def border_temperature_boundary_conditions(mesh):
    """
    Function that returns three vectors for the standard boundary condition for the Poisson temperature problem.
    :param mesh: Mesh object to be used to obtain the points information.
    """
    import numpy as np
    from collections import OrderedDict as oD

    # Acquiring borders
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


def check_method_call(*args):
    """
    Tests if the arguments are valid
    :param args: Any argument group that required for the method
    """
    for arg in args:
        try:
            if len(arg):
                return
        except TypeError:
            if arg:
                return
    raise ValueError("Method called incorrectly, please read the documentation and try changing the arguments.")


def create_new_surface(*imported_points, lt_version=True):
    """
    Create new surface
    :return: Element information
    """
    import numpy as np
    import scipy.spatial as dl

    x, y, ien = 0, 0, 0

    if imported_points:
        # Custom geometry.
        imported_points = np.array(imported_points[0])
        delauney_surfaces = dl.Delaunay(imported_points[:, :2])

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
            x, y, ien = use_meshio(None, geom)

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

            x, y, ien = use_meshio(None, geom)

    return x, y, ien


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
    import scipy.spatial as sp

    # Saving mesh as .vtk exportable file
    # points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geometry)
    # meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data) # TODO: FIX LIBRARY

    if fnmatch.fnmatch(filename, "*.msh"):
        points = meshio.read(filename).points
    else:
        points = meshio.read(filename + ".msh").points

    x_ = points[:, 0]
    y_ = points[:, 1]
    ien_ = sp.Delaunay(points[:, :2]).simplices

    return x_, y_, ien_


def get_dt(mesh):
    """
    Calculates optimal dt based on mesh.
    :param mesh: Mesh object that contains element and point informations.
    :return: Optimal dt.
    """
    import numpy as np
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


def get_matrices(mesh):
    """
    Function that generates the algebraic components of the solution method
    :param mesh: Mesh object to be used to generate the matrices
    :return:
    """
    import numpy as np
    from scipy import sparse

    # Default parameters:
    thickness = 1.0

    # Matrices:
    gx_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix (x component)
    gy_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix (y component)
    m_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Mass matrix

    for elem in mesh.ien:
        x = mesh.x[elem]
        y = mesh.y[elem]

        area = ((x[0] * y[1] - x[1] * y[0]) +
                (x[1] * y[2] - x[2] * y[1]) +
                (x[2] * y[0] - x[0] * y[2])) / 2.0

        # B = [b_i, b_j, b_k]
        b = np.array([y[1] - y[2],
                      y[2] - y[0],
                      y[0] - y[1]])

        # C = [c_i, c_j, c_k]
        c = np.array([x[2] - x[1],
                      x[0] - x[2],
                      x[1] - x[0]])

        #       [b_i*b_i, b_i*b_j, b_i*b_k],
        # K_x = [b_j*b_i, b_j*b_j, b_j*b_k],
        #       [b_k*b_i, b_k*b_j, b_k*b_k]
        k_x = np.array([
                        [b[0] * b[0], b[0] * b[1], b[0] * b[2]],
                        [b[0] * b[1], b[1] * b[1], b[1] * b[2]],
                        [b[0] * b[2], b[1] * b[2], b[2] * b[2]]
                       ])

        #       [c_i*c_i, c_i*c_j, c_i*c_k],
        # K_y = [c_j*c_i, c_j*c_j, c_j*c_k],
        #       [c_k*c_i, c_k*c_j, c_k*c_k]
        k_y = np.array([
                        [c[0] * c[0], c[0] * c[1], c[0] * c[2]],
                        [c[0] * c[1], c[1] * c[1], c[1] * c[2]],
                        [c[0] * c[2], c[1] * c[2], c[2] * c[2]]
                       ])

        #     [2, 1, 1],
        # M = [1, 2, 1],
        #     [1, 1, 2]
        m = (area / 12.0) * np.array([[2, 1, 1],
                                      [1, 2, 1],
                                      [1, 1, 2]])

        for i in range(3):  # Used so because of the triangular elements
            for j in range(3):
                gx_sparse[elem[i], elem[j]] += (thickness / (4.0 * area)) * k_x[i][j]
                gy_sparse[elem[i], elem[j]] += (thickness / (4.0 * area)) * k_y[i][j]
                m_sparse[elem[i], elem[j]] += m[i][j]

    return gx_sparse, gy_sparse, m_sparse


def apply_space_boundary_conditions(mesh, matrix_a, vector_b):
    """
    Performs the evaluation of boundary conditions and applies the changes to the main matrix and vector.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param matrix_a: The coefficients matrix A in the linear solution method [A*x = b].
    :param vector_b: The results vector b in the linear solution method [A*x = b].
    """
    for _relative_index, _column_index in enumerate(mesh.space_boundary_conditions.point_index_vector):
        if mesh.space_boundary_conditions.type_of_condition_vector[_relative_index]:
            # Dirichlet Treatment
            if matrix_a is not None:
                # Skip applying matrix a values for loops
                for _line_index in matrix_a.tocsc()[:, _column_index].indices:
                    vector_b[_line_index, 0] -= (matrix_a[_line_index, _column_index] *
                                                 mesh.space_boundary_conditions.values_vector[_relative_index])
                    matrix_a[_line_index, _column_index] = 0.
                    matrix_a[_column_index, _line_index] = 0.

                matrix_a[_column_index, _column_index] = 1.

            vector_b[_column_index, 0] = mesh.space_boundary_conditions.values_vector[_relative_index]
        else:
            # Neumann Treatment
            vector_b[_column_index, 0] += mesh.space_boundary_conditions.values_vector[_relative_index]


def apply_initial_boundary_conditions(mesh, vector_v):
    """
    Performs the evaluation of boundary conditions and applies initial values to the vector.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param vector_v: The vector v of initial values for the transient solution.
    """
    for _relative_index, _point in enumerate(mesh.time_boundary_conditions.point_index_vector):
        if mesh.time_boundary_conditions.type_of_condition_vector[_relative_index]:
            vector_v[0, _point] = mesh.time_boundary_conditions.values_vector[_relative_index]


def solve_poisson(mesh, permanent_solution=True, k_coef=0., k_coef_x=1.0, k_coef_y=1.0, q=0, total_time=1.):
    # TODO: REMOVE TOTAL TIME, CALCULATE BY DIFFERENCE BETWEEN FRAMES
    """
    Solves the mesh defined 2D Poisson equation problem:
        DT = -∇(k.∇T) + Q   ->   (M + K)*T_i^n = M*T_i^n-1 + M*Q_i
        dt                       dt              dt
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param permanent_solution: Parameter that defines if the solution will be calculated for the transient (True) or
                                permanent (False) problem.
    :param k_coef: Thermal conductivity coefficient both axis.
    :param k_coef_x: Thermal conductivity coefficient for x axis.
    :param k_coef_y: Thermal conductivity coefficient for y axis.
    :param q: Heat generation for each point.
    :return: Temperature value for each point in the mesh.
    """
    from scipy import sparse
    import scipy.sparse.linalg as linalg
    import numpy as np

    k_coef_x = k_coef or k_coef_x
    k_coef_y = k_coef or k_coef_y

    # ------------------ Contingency -----------------------------------------------------------------------------------
    if not (len(mesh.x) and len(mesh.y) and len(mesh.ien)):
        raise ValueError("The mesh is empty. Try using import_point_structure() before solving.")

    try:
        for boundary_condition in [mesh.space_boundary_conditions, mesh.time_boundary_conditions]:
            if (not isinstance(boundary_condition, Mesh.BoundaryConditions) or
                    not (len(boundary_condition.point_index_vector) and
                         len(boundary_condition.values_vector) and
                         len(boundary_condition.type_of_condition_vector))):
                raise ValueError("There are no boundary conditions defined. "
                                 "Try using mesh.boundary_condition.set_new_boundary_conditions() before solving.")
    except TypeError:
        raise ValueError("There are no boundary conditions defined. "
                         "Try using mesh.boundary_conditions.set_new_boundary_conditions() before solving.")

    # --- Defining the Matrices ----------------------------------------------------------------------------------------
    gx_matrix, gy_matrix, m_matrix = get_matrices(mesh)  # Stiffness (G_x, G_y) and Mass matrices (M)
    k_matrix = (k_coef_x * gx_matrix) + (k_coef_y * gy_matrix)  # K_xy = k_coef_x.G_x + k_coef_y.G_y
    q_matrix = sparse.lil_matrix((mesh.size, 1))  # Heat generation
    if isinstance(q, ComplexPointList):
        for _relative_index, _q in enumerate(q.indexes):
            q_matrix[_q] = q.values[_relative_index]

    if permanent_solution:
        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        b_vector = sparse.lil_matrix(m_matrix.dot(q_matrix))
        apply_space_boundary_conditions(mesh, k_matrix, b_vector)

        # --------------------------------- Solver ---------------------------------------------------------------------
        return linalg.spsolve(k_matrix.tocsc(), b_vector)

    else:
        # Optimal dt based on size of elements
        dt = get_dt(mesh)

        #      A = M/dt + K
        a_matrix = m_matrix / dt + k_matrix

        # First frame of the solution (time = 0)
        t_vector = sparse.lil_matrix((1, mesh.size))
        apply_initial_boundary_conditions(mesh, t_vector)

        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        # TODO: REFACTOR BOUNDARY CONDITION APPLY TO ALLOW DIRICHLET SOLUTION
        for _relative_index, _point in enumerate(mesh.space_boundary_conditions.point_index_vector):
            if mesh.space_boundary_conditions.type_of_condition_vector[_relative_index]:
                for _column_index in a_matrix.tocsr()[_point, :].indices:
                    a_matrix[_point, _column_index] = 0.
                a_matrix[_point, _point] = 1.
            else:
                q_matrix[_point, 0] += mesh.space_boundary_conditions.values_vector[_relative_index]

        frames = [np.ravel(t_vector.toarray())]
        for _frame_index in range(int(total_time / dt)):
            #      b = M * Q_i + M/dt * T_i^n-1
            b_vector = sparse.lil_matrix(m_matrix.dot(q_matrix) + m_matrix.dot(t_vector.reshape(-1, 1)) / dt)

            # --------------------------------- Boundary conditions treatment ------------------------------------------
            #     b += C.C. Dirichlet/Neumann
            apply_space_boundary_conditions(mesh, None, b_vector)

            #  A * x = b   ->   x = solve(A, b)
            t_vector = linalg.spsolve(a_matrix, b_vector)
            frames.append(t_vector)

        return frames


def solve_poiseuille(mesh, nu_coef=1.0):
    """
    Solves the mesh defined 2D Poiseuille equation problem:
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :return: Velocity vectors and pressure values for each point in the mesh.
    """
    import numpy as np
    from scipy import sparse
    import scipy.sparse.linalg as linalg
    # TODO: CONTINUE SOLUTION

    dt = get_dt(mesh)

    # --- Defining the Matrices ----------------------------------------------------------------------------------------
    gx_matrix, gy_matrix, m_matrix = get_matrices(mesh)

    k_matrix = (gx_matrix + gy_matrix)  # K_xy = G_x + G_y

    velocity_vector = sparse.csc_matrix((2, mesh.size))
    omega_vector = sparse.lil_matrix((0, mesh.size))

    # --------------------------------- Boundary conditions treatment --------------------------------------------------
    apply_initial_boundary_conditions(mesh, velocity_vector)

    # --------------------------------- Solve Loop ---------------------------------------------------------------------
    for i in range(10):
        a_matrix = (m_matrix / dt + nu_coef * k_matrix)
        b_vector = -velocity_vector.dot(gx_matrix.dot(omega_vector)) + (m_matrix / dt).dot(omega_vector)

        omega_vector = linalg.spsolve(a_matrix, b_vector)

        phi_vector = linalg.spsolve(k_matrix.tocsc(), m_matrix.dot(omega_vector))

        velocity_vector[0, :] = gy_matrix.dot(phi_vector)
        velocity_vector[1, :] = -gx_matrix.dot(phi_vector)


# -- Classes -----------------------------------------------------------------------------------------------------------

class Mesh:
    """
    Mesh element to be used in the calculations
    """
    x, y, ien = [], [], []
    name = "default"
    size = 0

    class BoundaryConditions:
        """
        Boundary conditions of the simulation
        """
        point_index_vector = []
        values_vector = []
        type_of_condition_vector = []

        def set_new_boundary_conditions(self, *vect_argm, point_index=range(0), values=0., type_of_boundary=True):
            """
            Sets the boundary conditions for the mesh
            :param vect_argm: Vector(s) of boundary conditions
            :param point_index: Vector of oreder of points
            :param values: Values of the condition in each point
            :param type_of_boundary: Type True for Dirichlet and False for Neumann
            """
            import numpy as np

            if vect_argm:
                check_method_call(vect_argm)
            else:
                check_method_call(point_index)

            try:
                if isinstance(vect_argm[0], list) or isinstance(values, np.ndarray):
                    array = np.array(vect_argm[0])

                    if len(vect_argm[0][0]) == 3:
                        self.point_index_vector = array[:, 0]
                        self.values_vector = array[:, 1]
                        self.type_of_condition_vector = list(map(bool, array[:, 2]))
                        return

                    if len(vect_argm[0][0]) == 2:
                        self.point_index_vector = array[:, 0]
                        self.values_vector = array[:, 1]
                        self.type_of_condition_vector = [True] * len(self.point_index_vector)
                        return

            except (TypeError, IndexError):
                self.point_index_vector = point_index
                if isinstance(values, list) or isinstance(values, np.ndarray):
                    if len(values) == len(point_index):
                        self.values_vector = np.array(values)
                    else:
                        raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                         "types and definitions")
                else:
                    self.values_vector = np.array([values] * len(point_index))

                if isinstance(type_of_boundary, list) or isinstance(type_of_boundary, np.ndarray):
                    if len(type_of_boundary) == len(point_index):
                        self.type_of_condition_vector = list(map(bool, np.array(type_of_boundary)))

                    else:
                        raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                         "types and definitions")
                else:
                    self.type_of_condition_vector = np.array([type_of_boundary] * len(point_index))

    def __init__(self, name="untitled", default_boundary_conditions=True):
        """
        Class constructor, initializes geometry.
        :param name: Mesh's main name.
        :param default_boundary_conditions: Determine if the default conditions are to be applied
        """
        import os
        self.import_point_structure(import_mesh_file=self.name)

        try:
            os.chdir("./results/{0}/".format(self.name))
        except FileNotFoundError:
            os.mkdir("./results/{0}".format(self.name))
            os.chdir("./results/{0}/".format(self.name))
        # TODO: COPY ORIGINAL .msh FILE TO NEW DIRECTORY

        self.space_boundary_conditions = self.BoundaryConditions()
        self.time_boundary_conditions = self.BoundaryConditions()

        # Default Boundary coditions declaration
        if default_boundary_conditions:
            max_x = max(self.x)
            max_y = max(self.y)
            min_x = min(self.x)
            min_y = min(self.y)
            vector = []

            for index, coord_x in enumerate(self.x):
                if (coord_x == max_x) or (self.y[index] == max_y):
                    # Boundaries are defined to be the upper and right sides with value 100
                    vector.append([index, 1., True])

            for index, coord_x in enumerate(self.x):
                if (coord_x == min_x) or (self.y[index] == min_y):
                    # Boundaries are defined to be the lower and left sides with value 0
                    vector.append([index, 0., True])

            self.space_boundary_conditions.set_new_boundary_conditions(vector)
            self.time_boundary_conditions.set_new_boundary_conditions(point_index=range(self.size), values=0.)

    def import_point_structure(self, *args, points=False, light_version=True, import_mesh_file=""):
        """
        Imports points position to create mesh from source file.
        :param args: Name of source file, defaults to "points.txt".
        :param points: Custom set of defined points, in form of list of shape (x, 2).
        :param light_version: Defines if script will use Gmsh to obtain elements.
        :param import_mesh_file: Use public library meshio to import already generated ".msh" file.
        """
        if not args:
            filename = "points.txt"
        else:
            filename = args[0]

        if isinstance(points, list):
            surface = create_new_surface(points, lt_version=light_version)
            self.space_boundary_conditions = self.BoundaryConditions()
            self.time_boundary_conditions = self.BoundaryConditions()

        else:
            try:
                if import_mesh_file:
                    surface = use_meshio("results/{0}".format(import_mesh_file), None)

                else:
                    with open(filename, "r") as arq:
                        points = []
                        for line in arq:
                            points.append([float(i) for i in line.split(";")] + [0.0])
                        surface = create_new_surface(points, lt_version=light_version)

            except FileNotFoundError:
                print("File not found, generating new default mesh.")
                # surface = create_new_surface(lt_version=False)  # TODO: FIX LIBRARY
                surface = use_meshio("results/untitled", None)

        self.x, self.y, self.ien = surface
        self.size = len(self.x)

    def output(self, result_vector, extension="VTK", dt=0., data_name = "Temperature"):
        """
        Export result to .vtk or .csv file
        :param result_vector: The vector of the value in each point
        :param extension: File extension
        :param dt: Value of time between frames
        :param data_name: Data name for display.
        """
        import os
        import fnmatch

        try:
            os.chdir("./paraview_results/")
        except FileNotFoundError:
            pass

        number_elements = len(self.ien)

        if extension == "CSV":
            # ------------------------- Saving results to CSV file ----------------------------------------------------
            with open("{0}_results.csv".format(self.name), "w") as arq:
                arq.write("{0}, Points:0, Points:1, Points:2\n".format(data_name))
                for i in range(self.size):
                    arq.write("{0},{1},{2},{3}\n".format(result_vector[i], self.x[i], self.y[i], 0))

        if extension == "VTK":
            try:
                size = len(result_vector[0])
                # -------- Deleting previous results -------------------------------------------------------------------
                list(map(os.remove, [file for file in os.listdir('.') if fnmatch.fnmatch(file, 'results_*.vtk')]))

                # --------- Saving multiple results to VTK files -------------------------------------------------------
                for j in range(len(result_vector)):
                    with open("{0}_results_{1}.vtk".format(self.name, j), "w") as arq:
                        # ------------------------------------ Header --------------------------------------------------
                        arq.write(
                            "# vtk DataFile Version 3.0\n{0}\n{1}\n\nDATASET {2}\n".format("LucasCarvalhoTCC Results",
                                                                                           "ASCII", "POLYDATA"))
                        arq.write("FIELD FieldData 1\nTIME 1 1 double\n{}\n".format(dt))
                        # ------------------------------------ Points coordinates --------------------------------------
                        arq.write("\nPOINTS {0} {1}\n".format(size, "float"))
                        for i in range(size):
                            arq.write("{0} {1} 0.0\n".format(self.x[i], self.y[i]))
                        # --------------------------------------- Cells ------------------------------------------------
                        arq.write("\nPOLYGONS {0} {1}\n".format(number_elements, number_elements * 4))
                        for i in range(number_elements):
                            arq.write("{0} {1} {2} {3}\n".format(3, self.ien[i][0], self.ien[i][1], self.ien[i][2]))
                        # ------------------------------------ Data in each point --------------------------------------
                        arq.write("\nPOINT_DATA {0}\n\nSCALARS {1} float 1\n".format(size, data_name))
                        arq.write("\nLOOKUP_TABLE {0}\n".format(data_name))
                        for i in range(size):
                            arq.write("{}\n".format(result_vector[j][i]))

            except TypeError:
                size = self.size
                # ------------------------- Saving results to VTK file -------------------------------------------------
                with open("{0}_results.vtk".format(self.name), "w") as arq:
                    # ------------------------------------ Header ------------------------------------------------------
                    arq.write(
                        "# vtk DataFile Version 3.0\n{0}\n{1}\n\nDATASET {2}\n".format("LucasCarvalhoTCC Results",
                                                                                       "ASCII", "POLYDATA"))
                    # ------------------------------------ Points coordinates ------------------------------------------
                    arq.write("\nPOINTS {0} {1}\n".format(size, "float"))
                    for i in range(size):
                        arq.write("{0} {1} 0.0\n".format(self.x[i], self.y[i]))
                    # --------------------------------------- Cells ----------------------------------------------------
                    arq.write("\nPOLYGONS {0} {1}\n".format(number_elements, number_elements * 4))
                    for i in range(number_elements):
                        arq.write("{0} {1} {2} {3}\n".format(3, self.ien[i][0], self.ien[i][1], self.ien[i][2]))
                    # ----------------------------------- Data in each point--------------------------------------------
                    arq.write("\nPOINT_DATA {0}\n\nSCALARS {1} float 1\n".format(size, data_name))
                    arq.write("\nLOOKUP_TABLE {0}\n".format(data_name))
                    for i in range(size):
                        arq.write("{}\n".format(result_vector[i]))

    def show_geometry(self, names=False, rainbow=False, save=True):
        """
        Display mesh geometry on screen using matplotlib.
        :param names: Show the index of each point next to it.
        :param rainbow: Color in the edge of each element in a different color.
        :param save: Save generated image.
        :return:
        """
        import numpy as np
        import matplotlib.pyplot as plt

        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        if names:
            for _index in range(self.size):
                plt.gca().annotate(_index, (self.x[_index], self.y[_index]))

        if rainbow:
            plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.hsv(np.linspace(0.0, 1.0, len(self.ien)))))

            for element in self.ien:
                plot_coordinates = (self.x[element[0]], self.x[element[1]], self.x[element[2]], self.x[element[0]]), \
                                   (self.y[element[0]], self.y[element[1]], self.y[element[2]], self.y[element[0]])

                plt.plot(plot_coordinates[0], plot_coordinates[1])

        else:
            plt.triplot(self.x, self.y, triangles=self.ien[0])

        if save:
            plt.savefig("{0}_mesh".format(self.name))

        plt.show()

    def show_solution(self, solution_vector):
        """
        Display 3D solution of the mesh geometry.
        :param solution_vector: Vector that contains the value of the solution for each point in the mesh.
        :return: Display image.
        """
        from matplotlib import pyplot
        from mpl_toolkits.mplot3d import Axes3D

        check_method_call(solution_vector)
        try:
            if len(solution_vector) != self.size:
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Solution must be a vector")

        self.output(solution_vector)

        fig = pyplot.gcf()
        axes = Axes3D(fig)
        surf = axes.plot_trisurf(self.x, self.y, solution_vector, cmap="jet")
        axes.view_init(90, 270)
        fig.colorbar(surf, shrink=0.4, aspect=9)

        pyplot.savefig("{0}_permanent_results".format(self.name))

        return pyplot.show()

    def show_animated_solution(self, frames_vector, dt=0.):
        """
        Display animated version of the 3D solution.
        :param frames_vector: Vector which each element contains a vector with the value of the solution for each point
                              in the mesh.
        :param dt: Time between each frame, if not specified the animation won't be saved.
        :return: Display image.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib.animation import FuncAnimation

        check_method_call(frames_vector)
        try:
            if np.any([len(i) != self.size for i in frames_vector]):
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                             "{0}".format(self.size))

        fig = plt.gcf()
        axes = Axes3D(fig)
        _min_value = np.min(frames_vector)
        _max_value = np.max(frames_vector)
        surf = axes.plot_trisurf(self.x, self.y, frames_vector[0], cmap="jet", vmin=_min_value, vmax=_max_value)
        fig.colorbar(surf, shrink=0.4, aspect=9)
        frame_index = 0

        def update(_current_frame):
            plt.cla()
            axes.plot_trisurf(self.x, self.y, _current_frame, cmap="jet", vmin=_min_value, vmax=_max_value)
            axes.set_zlim3d([_min_value, _max_value])
            if dt:
                # global frame_index
                axes.text2D(0., 0.9, "Frame: {0}\nTime: {1}".format(frame_index, frame_index * dt),
                            transform=axes.transAxes)
                # frame_index += 1
                # TODO: ADD TIME INFORMATION TO FRAME
            return

        animation = FuncAnimation(fig, update, frames=frames_vector, interval=100, save_count=False)

        if dt:
            self.output(frames_vector, dt=dt)
            animation.save("{0}_transient_results.gif".format(self.name), dpi=80, writer='imagemagick')

        return plt.show()


class ComplexPointList:
    """
    Class that defines a list of points with an index value and a property value each.
    """
    indexes = []
    values = []

    def __init__(self, _indexes, _values):
        """
        Class constructor that generates one list for each property
        :param _indexes: The point index in the ordered list of the mesh point coordinates
        :param _values: The value of the desired property in each point.
        """
        check_method_call(_indexes)

        self.indexes = _indexes

        try:
            if isinstance(list(_values), list):
                self.values = _values
            else:
                raise TypeError

        except TypeError:
            for i in range(len(_indexes)):
                self.values.append(_values)

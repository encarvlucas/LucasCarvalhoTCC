# -*- coding: utf-8 -*-
# -------------------------------------- by: LUCAS CARVALHO DE SOUSA ---------------------------------------------------
# Esta biblioteca foi criada para o a solução de sistemas diferenciais através do Método de Elementos Finitos
# This library was created for educational purposes, it is meant to be used for the solution of differential equation
# systems
__author__ = "Lucas Carvalho de Sousa"
__author_email__ = "encarvlucas@gmail.com"
__website__ = "https://github.com/encarvlucas/LucasCarvalhoTCC"

# -- Functions ---------------------------------------------------------------------------------------------------------

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
            x, y, ien = use_meshio(geom)

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

            x, y, ien = use_meshio(geom)

    return x, y, ien


def use_meshio(geometry):
    """
    Use MeshIO library for creating the mesh point structure from Gmsh
    :param geometry: A Geometry object from the PyGmsh library
    :return: x, y, ien - The mesh point structure
    """
    import pygmsh
    import meshio
    import scipy.spatial as sp

    # Saving mesh as .vtk exportable file
    # points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geometry)
    # meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data) # TODO: FIX LIBRARY
    points = meshio.read("untitled.msh").points

    x_ = points[:, 0]
    y_ = points[:, 1]
    ien_ = sp.Delaunay(points[:, :2]).simplices

    return x_, y_, ien_


def get_matrices(mesh):
    """
    Function that generates the algebraic components of the solution method
    :param mesh: Mesh object to be used to generate the matrices
    :return:
    """
    import numpy as np
    from scipy import sparse

    # Default parameters:
    k_coef_x = 1.0
    k_coef_y = 1.0
    thickness = 1.0

    # Matrices:
    k_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix
    m_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Mass matrix

    for elem in mesh.ien:
        x = mesh.x[elem]
        y = mesh.y[elem]

        area = ((x[0] * y[1] - x[1] * y[0]) +
                (x[1] * y[2] - x[2] * y[1]) +
                (x[2] * y[0] - x[0] * y[2])) / 2.0

        b = np.array([y[1] - y[2],
                      y[2] - y[0],
                      y[0] - y[1]])

        c = np.array([x[2] - x[1],
                      x[0] - x[2],
                      x[1] - x[0]])

        k = -(thickness / (4.0 * area)) * (k_coef_x * np.array([
                                                                [b[0] * b[0], b[0] * b[1], b[0] * b[2]],
                                                                [b[0] * b[1], b[1] * b[1], b[1] * b[2]],
                                                                [b[0] * b[2], b[1] * b[2], b[2] * b[2]]
                                                               ]) +
                                           k_coef_y * np.array([
                                                                [c[0] * c[0], c[0] * c[1], c[0] * c[2]],
                                                                [c[0] * c[1], c[1] * c[1], c[1] * c[2]],
                                                                [c[0] * c[2], c[1] * c[2], c[2] * c[2]]
                                                               ]))
        m = (area / 12.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])

        for i in range(3):  # Used so because of the triangular elements
            for j in range(3):
                k_sparse[elem[i], elem[j]] += k[i][j]
                m_sparse[elem[i], elem[j]] += m[i][j]

    return k_sparse, m_sparse


def solve(mesh, permanent_solution=True, dt=0.01, total_time=1.):
    """
    Solves the mesh defined 2D problem
    :return: The solution for the permanent problem
    """
    import numpy as np
    from scipy import sparse
    import scipy.sparse.linalg as linalg

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

    if permanent_solution:
        # --- Defining the Matrices ------------------------------------------------------------------------------------
        q_matrix = sparse.lil_matrix((mesh.size, 1))  # Heat generation
        k_matrix, m_matrix = get_matrices(mesh)  # Stiffness and Mass matrices

        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        for relative_index, column_index in enumerate(mesh.space_boundary_conditions.point_index_vector):
            if mesh.space_boundary_conditions.type_of_condition_vector[relative_index]:
                # Dirichlet Treatment
                for line_index in k_matrix.tocsc()[:, column_index].indices:
                    q_matrix[line_index, 0] -= (k_matrix[line_index, column_index] *
                                                mesh.space_boundary_conditions.values_vector[relative_index])
                    k_matrix[line_index, column_index] = 0.
                    k_matrix[column_index, line_index] = 0.

                k_matrix[column_index, column_index] = 1.
                q_matrix[column_index, 0] = mesh.space_boundary_conditions.values_vector[relative_index]
            else:
                # Neumann Treatment
                q_matrix[column_index, 0] -= mesh.space_boundary_conditions.values_vector[relative_index]

        # --------------------------------- Solver ---------------------------------------------------------------------
        return linalg.spsolve(k_matrix.tocsc(), q_matrix)

    else:
        # --- Defining the Matrices-------------------------------------------------------------------------------------
        q_matrix = sparse.lil_matrix((mesh.size, 1))  # Heat generation
        k_matrix, m_matrix = get_matrices(mesh)  # Stiffness and Mass matrices

        a_matrix = m_matrix / dt - k_matrix

        # First frame of the solution (time = 0)
        initial = sparse.lil_matrix((1, mesh.size))
        for relative_index, point in enumerate(mesh.time_boundary_conditions.point_index_vector):
            if mesh.time_boundary_conditions.type_of_condition_vector[point]:
                initial[0, point] = mesh.time_boundary_conditions.values_vector[relative_index]

        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        def boundary_treatment(vec):
            for _relative_index, _point in enumerate(mesh.space_boundary_conditions.point_index_vector):
                vec[_point] = mesh.time_boundary_conditions.values_vector[_relative_index]
            return vec

        for relative_index, point in enumerate(mesh.space_boundary_conditions.point_index_vector):
            if mesh.space_boundary_conditions.type_of_condition_vector[relative_index]:
                a_matrix[point, :] = 0.
                a_matrix[point, point] = 1.
            else:
                q_matrix[point, 0] -= mesh.space_boundary_conditions.values_vector[relative_index]
                # TODO: ASK GUSTAVO WHY IT WAS 1/dt INSTEAD OF dt

        # --------------------------------- Solver ---------------------------------------------------------------------
        t_matrix = initial
        frames = [np.ravel(initial.toarray())]
        for frame_index in range(int(total_time / dt)):
            b_matrix = m_matrix.dot(q_matrix) / dt + sparse.lil_matrix(m_matrix.dot(t_matrix.reshape(-1, 1)))
            b_matrix = boundary_treatment(b_matrix)
            t_matrix = linalg.spsolve(a_matrix, b_matrix)
            frames.append(t_matrix)

        return frames


# -- Classes -----------------------------------------------------------------------------------------------------------

class Mesh:
    """
    Mesh element to be used in the calculations
    """
    x, y, ien = [], [], []
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

            # TODO: INCREASE CHECKS
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
                        self.type_of_condition_vector = array[:, 2]
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
                        self.type_of_condition_vector = np.array(type_of_boundary)

                    else:
                        raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                         "types and definitions")
                else:
                    self.type_of_condition_vector = np.array([type_of_boundary] * len(point_index))

    def __init__(self, default_boudary_conditions=True):
        """
        Class constructor,
        initializes geometry
        :param default_boudary_conditions: Determine if the default conditions are to be applied
        """
        self.import_point_structure()

        self.space_boundary_conditions = self.BoundaryConditions()
        self.time_boundary_conditions = self.BoundaryConditions()

        # Default Boundary coditions declaration
        if default_boudary_conditions:
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

    def import_point_structure(self, *args, points=False, light_version=True):
        """
        Imports points position to create mesh from source file
        :param args: Name of source file, defaults to points.txt
        :param points: Custom set of defined points, in form of list of shape (x, 2)
        :param light_version: Defines if script will use Gmsh to obtain elements
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
                with open(filename, "r") as arq:
                    points = []
                    for line in arq:
                        points.append([float(i) for i in line.split(";")] + [0.0])
                    surface = create_new_surface(points)

            except FileNotFoundError:
                surface = create_new_surface(lt_version=False)

        self.x, self.y, self.ien = surface
        self.size = len(self.x)

    def output(self, result_vector, extension="VTK", dt=0.):
        """
        Export result to .vtk or .csv file
        :param result_vector: The vector of the value in each point
        :param extension: File extension
        :param dt: Value of time between frames
        """
        import os
        import fnmatch

        os.chdir("results/")
        number_elements = len(self.ien)
        data_name = "Temperature"

        if extension == "CSV":
            # ------------------------- Saving results to CSV file ----------------------------------------------------
            with open("results.csv", "w") as arq:
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
                    with open("results_{}.vtk".format(j), "w") as arq:
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
                with open("results.vtk", "w") as arq:
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

    def show_geometry(self, rainbow=False):
        """
        Display mesh geometry on screen using matplotlib
        """
        import numpy as np
        import matplotlib.pyplot as plt

        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        if rainbow:
            plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.hsv(np.linspace(0.0, 1.0, len(self.ien)))))

            for element in self.ien:
                plot_coordinates = (self.x[element[0]], self.x[element[1]], self.x[element[2]], self.x[element[0]]), \
                                   (self.y[element[0]], self.y[element[1]], self.y[element[2]], self.y[element[0]])

                plt.plot(plot_coordinates[0], plot_coordinates[1])

        else:
            plt.triplot(self.x, self.y, triangles=self.ien[0])

        plt.show()

    def show_solution(self, solution_vector):
        """
        Display 3D solution of the mesh geometry
        :param solution_vector: Vector that contains the value of the solution for each point in the mesh
        :return: Display image
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

        return pyplot.show()

    def show_animated_solution(self, frames_vector, dt=0.):
        """
        Display animated version of the 3D solution
        :param frames_vector: Vector which each element contains a vector with the value of the solution for each point
            in the mesh
        :param dt: Time between each frame
        :return:
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

        self.output(frames_vector, dt=dt)

        fig = plt.figure()
        axes = plt.gca()
        axes = Axes3D(fig)
        _min_value = np.min(frames_vector)
        _max_value = np.max(frames_vector)
        surf = axes.plot_trisurf(self.x, self.y, frames_vector[0], cmap="jet", vmin=_min_value, vmax=_max_value)
        fig.colorbar(surf, shrink=0.4, aspect=9)

        def update(current_frame):
            plt.cla()
            axes.plot_trisurf(self.x, self.y, current_frame, cmap="jet", vmin=_min_value, vmax=_max_value)
            axes.set_zlim3d([_min_value, _max_value])
            return

        anim = FuncAnimation(fig, update, frames=frames_vector, interval=100, save_count=False)

        return plt.show()

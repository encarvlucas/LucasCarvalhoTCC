# -*- coding: utf-8 -*-
# -------------------------------------- by: LUCAS CARVALHO DE SOUSA ---------------------------------------------------
# Esta biblioteca foi criada para o a solução de sistemas diferenciais através do Método de Elementos Finitos
# This library was created for educational purposes, it is meant to be used for the solution of differential equation
# systems
__author__ = "Lucas Carvalho de Sousa"


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
            ])

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

    # Saving mesh as .vtk exportable file
    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geometry)
    meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data)

    x_ = points[:, 0]
    y_ = points[:, 1]
    ien_ = cells["triangle"]

    return x_, y_, ien_


def solve(mesh, permanent_solution=True):
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
            if not (len(boundary_condition.point_index_vector) and
                    len(boundary_condition.values_vector) and
                    len(boundary_condition.type_of_condition_vector)):
                raise ValueError("There are no boundary conditions defined. "
                                 "Try using mesh.boundary_condition.set_new_boundary_conditions() before solving.")
    except TypeError:
        raise ValueError("There are no boundary conditions defined. "
                         "Try using mesh.boundary_conditions.set_new_boundary_conditions() before solving.")

    k_coef_x = 1.0
    k_coef_y = 1.0
    thickness = 1.0
    dt = 0.1
    total_time = 10.0

    def get_matrix():
        """
        Function that generates the algebraic components of the solution method
        :return:
        """
        k_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix
        m_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Mass matrix

        for elem in mesh.ien:
            x = mesh.x[elem]
            y = mesh.y[elem]

            a = ((x[0] * y[1] - x[1] * y[0]) +
                 (x[1] * y[2] - x[2] * y[1]) +
                 (x[2] * y[0] - x[0] * y[2])) / 2.0

            b = np.array([y[1] - y[2],
                          y[2] - y[0],
                          y[0] - y[1]])

            c = np.array([x[2] - x[1],
                          x[0] - x[2],
                          x[1] - x[0]])

            k = -(thickness / (4.0 * a)) * (k_coef_x * np.array([
                                                                 [b[0] * b[0], b[0] * b[1], b[0] * b[2]],
                                                                 [b[0] * b[1], b[1] * b[1], b[1] * b[2]],
                                                                 [b[0] * b[2], b[1] * b[2], b[2] * b[2]]
                                                                ]) +
                                            k_coef_y * np.array([
                                                                 [c[0] * c[0], c[0] * c[1], c[0] * c[2]],
                                                                 [c[0] * c[1], c[1] * c[1], c[1] * c[2]],
                                                                 [c[0] * c[2], c[1] * c[2], c[2] * c[2]]
                                                                ]))
            if permanent_solution:
                m = (a / 12.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])

            for i in range(3):  # Used so because of the triangular elements
                for j in range(3):
                    k_sparse[elem[i], elem[j]] += k[i][j]
                    if permanent_solution:
                        m_sparse[elem[i], elem[j]] += m[i][j]

        return k_sparse, m_sparse

    if permanent_solution:

        # --- Defining the Matrices-------------------------------------------------------------------------------------
        q_matrix = sparse.lil_matrix((mesh.size, 1))  # Heat generation
        k_matrix, m_matrix = get_matrix()  # Stiffness and Mass matrices

        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        k_matrix = k_matrix.tocsc()
        for index, point_index in enumerate(mesh.space_boundary_conditions.point_index_vector):
            for i in k_matrix[:, point_index].indices:
                if mesh.space_boundary_conditions.type_of_condition_vector[index]:
                    # Dirichlet Treatment
                    q_matrix[i, 0] -= k_matrix[i, point_index] * mesh.space_boundary_conditions.values_vector[index]
                k_matrix[i, point_index] = 0
                k_matrix[point_index, i] = 0
            k_matrix[point_index, point_index] = 1
            q_matrix[point_index, 0] = mesh.x[point_index] ** 2 + 1

        # --------------------------------- Solver ---------------------------------------------------------------------
        return linalg.spsolve(k_matrix, q_matrix)

    else:
        # --- Defining the Matrices-------------------------------------------------------------------------------------
        q_matrix = sparse.lil_matrix((mesh.size, 1))  # Heat generation
        k_matrix, m_matrix = get_matrix()  # Stiffness and Mass matrices

        a_matrix = m_matrix - k_matrix * dt

        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        def cc1(vec):
            for i in mesh.space_boundary_conditions.point_index_vector:
                vec[i] = mesh.time_boundary_conditions.values_vector[i]
            return vec

        def cc2():
            for i in mesh.space_boundary_conditions.point_index_vector:
                a_matrix[i, :] = 0
                a_matrix[i, i] = 1
            return

        # --------------------------------- Solver ---------------------------------------------------------------------
        cc2()
        frames = [mesh.time_boundary_conditions.values_vector]
        t_matrix = frames[-1]
        for t in range(int(total_time / dt)):
            b_matrix = dt * m_matrix.dot(q_matrix) + sparse.lil_matrix(m_matrix.dot(t_matrix)).transpose()
            b_matrix = cc1(b_matrix)
            t_matrix = linalg.spsolve(a_matrix, b_matrix)
            frames.append(t_matrix)

        return frames

    return 0


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

        def set_new_boundary_conditions(self, *vect_argm, point_index=range(0), values=0, type_of_boundary=True):
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
                if isinstance(vect_argm[0], list):
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
                if isinstance(values, list):
                    if len(values) == len(point_index):
                        self.values_vector = values
                    else:
                        raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                         "types and definitions")
                else:
                    self.values_vector = [values] * len(point_index)
                if isinstance(type_of_boundary, list):
                    if len(type_of_boundary) == len(point_index):
                        self.type_of_condition_vector = type_of_boundary
                    else:
                        raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                         "types and definitions")
                else:
                    self.type_of_condition_vector = [type_of_boundary] * len(point_index)

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
                    vector.append([index, 100, True])

            for index, coord_x in enumerate(self.x):
                if (coord_x == min_x) or (self.y[index] == min_y):
                    # Boundaries are defined to be the lower and left sides with value 0
                    vector.append([index, 0, True])

            self.space_boundary_conditions.set_new_boundary_conditions(vector)
            self.time_boundary_conditions.set_new_boundary_conditions(point_index=range(self.size), values=0)

    def import_point_structure(self, *args):
        """
        Imports points position to create mesh from source file
        :param args: Name of source file, defaults to points.txt
        """
        if not args:
            filename = "points.txt"
        else:
            filename = args[0]

        try:
            with open(filename, "r") as arq:
                points = []
                for line in arq:
                    points.append([float(i) for i in line.split(";")] + [0.0])
                surface = create_new_surface(points)

        except FileNotFoundError:
            surface = create_new_surface()

        self.x, self.y, self.ien = surface
        self.size = len(self.x)

    def show(self, rainbow=False):
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

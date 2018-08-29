# -*- coding: utf-8 -*-
# -------------------------------------- by: LUCAS CARVALHO DE SOUSA ---------------------------------------------------
# Esta biblioteca foi criada para o a solução de sistemas diferenciais através do Método de Elementos Finitos
# This library was created for educational purposes, it is meant to be used for the solution of diferential equation
# systems
__author__ = "Lucas Carvalho de Sousa"


# -- Functions ---------------------------------------------------------------------------------------------------------

def create_new_surface(*imported_points, lt_version=False):
    """
    Create new surface
    :return: Element information
    """
    import numpy as np
    import scipy.spatial as dl
    if not lt_version:
        import pygmsh
        geom = pygmsh.built_in.Geometry()

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
            for tri in delauney_surfaces.simplices:
                geom.add_polygon([[delauney_surfaces.points[tri[0]][0], delauney_surfaces.points[tri[0]][1], 0.0],
                                  [delauney_surfaces.points[tri[1]][0], delauney_surfaces.points[tri[1]][1], 0.0],
                                  [delauney_surfaces.points[tri[2]][0], delauney_surfaces.points[tri[2]][1], 0.0]])

    else:
        if not lt_version:

            # Default surface.
            geom.add_polygon([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ])

    if not lt_version:
        import meshio

        # Saving mesh as .vtk exportable file
        points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
        meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data)

        x = points[:, 0]
        y = points[:, 1]
        ien = cells["triangle"]

    return x, y, ien


def solve(mesh, permanent_solution=True):
    """
    Solves the mesh defined 2D problem
    :return: The solution for the permanent problem
    """
    import numpy as np
    from scipy import sparse

    if not (len(mesh.x) and len(mesh.y) and len(mesh.ien)):
        raise ValueError("The mesh is empty. Try using import_point_structure() before solving.")

    # if not len(mesh.boundary_conditions):
    #     raise ValueError("There are no boundary conditions defined. "
    #                      "Try using set_boundary_conditions() before solving.")

    k_coef_x = 1.0  # TODO: DEFINE OUTSIDE FUNCTION
    k_coef_y = 1.0
    thickness = 1.0

    if permanent_solution:
        # --- Defining the Matrices--------------------------------------------------------------------------------
        q_matrix = sparse.lil_matrix((1, mesh.size))  # Heat generation
        k_matrix = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix

        for elem in mesh.ien:
            x = mesh.x[elem]
            y = mesh.y[elem]

            A = ((x[0] * y[1] - x[1] * y[0]) +
                 (x[1] * y[2] - x[2] * y[1]) +
                 (x[2] * y[0] - x[0] * y[2]))

            b = np.array([y[1] - y[2],
                          y[2] - y[0],
                          y[0] - y[1]])

            c = np.array([x[2] - x[1],
                          x[0] - x[2],
                          x[1] - x[0]])

            k = -(thickness / (4 * A)) * (k_coef_x * np.array([[b[0] ** 2, b[0] * b[1], b[0] * b[2]],
                                                        [b[0] * b[1], b[1] ** 2, b[1] * b[2]],
                                                        [b[0] * b[2], b[1] * b[2], b[2] ** 2]]) +
                                          k_coef_y * np.array([[c[0] ** 2, c[0] * c[1], c[0] * c[2]],
                                                              [c[0] * c[1], c[1] ** 2, c[1] * c[2]],
                                                              [c[0] * c[2], c[1] * c[2], c[2] ** 2]]))

            for i in range(3):  # Used so because of the triangular elements
                for j in range(3):
                    k_matrix[elem[i], elem[j]] += k[i][j]

        # --------------------------------- Boundary conditions treatment Dirichlet ------------------------------------
        for point_index in []:  # TODO: FIND BOUNDARIE POINTS points where there are boundary conditions
            for i in k_matrix.nonzero():
                q_matrix[i] -= k_matrix[i, point_index] * ((mesh.x[point_index]) ** 2 + 1)
                #				/\ - valor da função no ponto
                k_matrix[i, point_index] = 0
                k_matrix[point_index, i] = 0
            k_matrix[point_index, point_index] = 1
            q_matrix[point_index] = mesh.x[point_index] ** 2 + 1

        # --------------------------------- Solver ---------------------------------------------------------------------
        T = np.linalg.solve(k_matrix, q_matrix)

    # TODO: APPLY TRANSIENT SOLUTION
    """
    def elemfintrans(Lx, Ly, nx, ny, k_condx, k_condy, esp, xy, IEN, dt, nt, T_0, cc_id):
        numpnts_x = nx + 1
        numpnts_y = ny + 1
        # --------------------------------Geração das matrizes---------------------------------------
        Q = np.zeros(len(xy))  # Geração de calor
        K = np.zeros((len(xy), len(xy)))
        M = np.copy(K)
        T = np.copy(T_0)

        for elem in IEN:
            x = xy[elem, 0]
            y = xy[elem, 1]
            A = ((x[0] * y[1] - x[1] * y[0]) +
                 (x[1] * y[2] - x[2] * y[1]) +
                 (x[2] * y[0] - x[0] * y[2])) / 2
            b = np.array([y[1] - y[2],
                          y[2] - y[0],
                          y[0] - y[1]])
            c = np.array([x[2] - x[1],
                          x[0] - x[2],
                          x[1] - x[0]])
            k = -(esp / (4.0 * A)) * (k_condx * np.array([[b[0] ** 2, b[0] * b[1], b[0] * b[2]],
                                                          [b[0] * b[1], b[1] ** 2, b[1] * b[2]],
                                                          [b[0] * b[2], b[1] * b[2], b[2] ** 2]])
                                      + k_condy * np.array([[c[0] ** 2, c[0] * c[1], c[0] * c[2]],
                                                            [c[0] * c[1], c[1] ** 2, c[1] * c[2]],
                                                            [c[0] * c[2], c[1] * c[2], c[2] ** 2]]))
            m = (A / 12.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    K[elem[i]][elem[j]] += k[i][j]
                    M[elem[i]][elem[j]] += m[i][j]

        A = M - K * dt

        # ---------------------------------Condição de contorno-----------------------------------------------
        def cc1(vec):
            for i in cc_id:
                vec[i] = T_0[i]
            return vec

        def cc2():
            for i in cc_id:
                A[i, :] = 0
                A[i][i] = 1
                # vec[i] = T_0[i]
            return

        # -----------------------------------Solução no tempo-------------------------------------------
        # PROCURAR método dos gradientes conjugados
        # Gmesh
        # Delauney
        cc2()
        frames = [T_0]
        for t in range(nt):
            B = dt * np.dot(M, Q) + np.dot(M, T)
            B = cc1(B)
            T = np.linalg.solve(A, B)
            frames.append(T)
    """

    return 0


# -- Classes -----------------------------------------------------------------------------------------------------------

class Mesh:
    """
    Mesh element to be used in the calculations
    """
    x, y, ien = [], [], []
    size = 0
    boundary_conditions = []

    def __init__(self):
        """
        Class constructor,
        initializes geometry
        """
        self.import_point_structure()
        self.set_boundary_conditions()

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

            for triangl in self.ien:
                plot_coordinates = (self.x[triangl[0]], self.x[triangl[1]], self.x[triangl[2]], self.x[triangl[0]]), \
                                   (self.y[triangl[0]], self.y[triangl[1]], self.y[triangl[2]], self.y[triangl[0]])

                plt.plot(plot_coordinates[0], plot_coordinates[1])

        else:
            plt.triplot(self.x, self.y, triangles=self.ien[0])

        plt.show()

    def set_boundary_conditions(self):
        """
        Sets boundary conditions,
        they are determined as a tuple of (point, Dirichlet_value)
        """
        # TODO: SET BOUNDARY CONDITIONS

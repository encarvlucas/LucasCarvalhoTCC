# -*- coding: utf-8 -*-
# -------------------------------------- by: LUCAS CARVALHO DE SOUSA ---------------------------------------------------
# Esta biblioteca foi criada para o uso academico para a solucao de sistemas diferenciais atraves do Metodo de Elementos
#  Finitos
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
    import numpy as np
    from collections import OrderedDict as oD

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
    values = np.zeros(len(indices)) + 1.0
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


def check_list_or_object(_list, _class):
    # TODO: CHECK IF ITS BEING USED CORRECTLY
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
    import numpy as np
    import scipy.spatial as sp

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

    delauney_surfaces = sp.Delaunay(points[:, :2])
    ien_ = delauney_surfaces.simplices

    return x_, y_, ien_, delauney_surfaces


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


def sparse_to_vector(vector):
    """
    Converts one dimensional sparse matrix to a vector array to allow more features.
    :param vector: Vector as a sparse matrix (x,1) or (1,x).
    :return: Vector as an one dimensional array (x,).
    """
    import numpy as np
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
    kx_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix (x component)
    ky_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Stiffness matrix (y component)
    gx_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Convection matrix (x component)
    gy_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Convection matrix (y component)
    m_sparse = sparse.lil_matrix((mesh.size, mesh.size))  # Mass matrix

    for elem in mesh.ien:
        x = mesh.x[elem]
        y = mesh.y[elem]

        area = get_area(x, y)

        # B = [b_i, b_j, b_k]
        b = np.array([y[1] - y[2],
                      y[2] - y[0],
                      y[0] - y[1]])

        # C = [c_i, c_j, c_k]
        c = np.array([x[2] - x[1],
                      x[0] - x[2],
                      x[1] - x[0]])

        #          [b_i*b_i, b_i*b_j, b_i*b_k],
        # K_x =  1 [b_j*b_i, b_j*b_j, b_j*b_k],
        #       4A [b_k*b_i, b_k*b_j, b_k*b_k]
        k_x = (thickness / (4.0 * area)) * np.array([
                                                     [b[0] * b[0], b[0] * b[1], b[0] * b[2]],
                                                     [b[0] * b[1], b[1] * b[1], b[1] * b[2]],
                                                     [b[0] * b[2], b[1] * b[2], b[2] * b[2]]
                                                    ])

        #          [c_i*c_i, c_i*c_j, c_i*c_k],
        # K_y =  1 [c_j*c_i, c_j*c_j, c_j*c_k],
        #       4A [c_k*c_i, c_k*c_j, c_k*c_k]
        k_y = (thickness / (4.0 * area)) * np.array([
                                                     [c[0] * c[0], c[0] * c[1], c[0] * c[2]],
                                                     [c[0] * c[1], c[1] * c[1], c[1] * c[2]],
                                                     [c[0] * c[2], c[1] * c[2], c[2] * c[2]]
                                                    ])

        #         [b_i, b_j, b_k],
        # G_x = 1 [b_i, b_j, b_k],
        #       6 [b_i, b_j, b_k]
        g_x = (1./6.) * np.array([
                                  [b[0], b[1], b[2]],
                                  [b[0], b[1], b[2]],
                                  [b[0], b[1], b[2]]
                                 ])

        #         [c_i, c_j, c_k],
        # G_y = 1 [c_i, c_j, c_k],
        #       6 [c_i, c_j, c_k]
        g_y = (1./6.) * np.array([
                                  [c[0], c[1], c[2]],
                                  [c[0], c[1], c[2]],
                                  [c[0], c[1], c[2]]
                                 ])

        #     [2, 1, 1],
        # M = [1, 2, 1],
        #     [1, 1, 2]
        m = (area / 12.0) * np.array([[2, 1, 1],
                                      [1, 2, 1],
                                      [1, 1, 2]])

        for i in range(3):  # Used so because of the triangular elements
            for j in range(3):
                kx_sparse[elem[i], elem[j]] += k_x[i][j]
                ky_sparse[elem[i], elem[j]] += k_y[i][j]
                m_sparse[elem[i], elem[j]] += m[i][j]
                gx_sparse[elem[i], elem[j]] += g_x[i][j]
                gy_sparse[elem[i], elem[j]] += g_y[i][j]

    return kx_sparse, ky_sparse, m_sparse, gx_sparse, gy_sparse


def apply_boundary_conditions(mesh, boundary_name, matrix_a, vector_b):
    """
    Performs the evaluation of boundary conditions and applies the changes to the main matrix and vector.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param boundary_name: The boundary name as defined for each problem.
    :param matrix_a: The coefficients matrix A in the linear solution method [A*x = b].
    :param vector_b: The results vector b in the linear solution method [A*x = b].
    """
    for _relative_index, _column_index in enumerate(mesh.boundary_conditions[boundary_name].point_index_vector):
        if mesh.boundary_conditions[boundary_name].type_of_condition_vector[_relative_index]:
            # Dirichlet Treatment
            if matrix_a is not None:
                # Skip applying matrix a values for loops
                for _line_index in matrix_a.tocsc()[:, _column_index].indices:
                    vector_b[_line_index, 0] -= (matrix_a[_line_index, _column_index] *
                                                 mesh.boundary_conditions[boundary_name].values_vector[_relative_index])
                    matrix_a[_line_index, _column_index] = 0.
                    matrix_a[_column_index, _line_index] = 0.

                matrix_a[_column_index, _column_index] = 1.

            vector_b[_column_index, 0] = mesh.boundary_conditions[boundary_name].values_vector[_relative_index]
        else:
            # Neumann Treatment
            vector_b[_column_index, 0] += mesh.boundary_conditions[boundary_name].values_vector[_relative_index]


def apply_initial_boundary_conditions(mesh, boundary_name, vector_v):
    """
    Performs the evaluation of boundary conditions and applies initial values to the vector.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param vector_v: The vector v of initial values for the transient solution.
    :param boundary_name: The boundary name as defined for each problem.
    """
    for _relative_index, _point in enumerate(mesh.boundary_conditions[boundary_name].point_index_vector):
        if mesh.boundary_conditions[boundary_name].type_of_condition_vector[_relative_index]:
            vector_v[_point] = mesh.boundary_conditions[boundary_name].values_vector[_relative_index]


def solve_poisson(mesh, permanent_solution=True, k_coef=0., k_coef_x=1.0, k_coef_y=1.0, q=None, dt=None, total_time=1.):
    """
    Solves the mesh defined 2D Poisson equation problem:
        DT = -∇(k*∇T) + Q   ->   (M + K).T_i^n =  M.T_i^n-1 + M.Q_i
        dt                       dt              dt
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param permanent_solution: Parameter that defines if the solution will be calculated for the transient (True) or
                                permanent (False) problem.
    :param k_coef: Thermal conductivity coefficient both axis.
    :param k_coef_x: Thermal conductivity coefficient for x axis.
    :param k_coef_y: Thermal conductivity coefficient for y axis.
    :param q: Heat generation for each point.
    :param dt: Value of time between frames.
    :param total_time: Length of time the calculation takes place (only necessary for transient solutions).
    :return: Temperature value for each point in the mesh.
    """
    from scipy import sparse
    import scipy.sparse.linalg as linalg

    k_coef_x = k_coef or k_coef_x
    k_coef_y = k_coef or k_coef_y

    # ------------------ Contingency -----------------------------------------------------------------------------------
    if not (len(mesh.x) and len(mesh.y) and len(mesh.ien)):
        raise ValueError("The mesh is empty. Try using import_point_structure() before solving.")

    try:
        for boundary_condition in [mesh.boundary_conditions["space"], mesh.boundary_conditions["time"]]:
            if (not isinstance(boundary_condition, BoundaryConditions) or
                    not (len(boundary_condition.point_index_vector) and
                         len(boundary_condition.values_vector) and
                         len(boundary_condition.type_of_condition_vector))):
                raise ValueError("There are no boundary conditions defined. "
                                 "Try using mesh.new_boundary_condition() before solving.")
    except (TypeError, KeyError):
        raise ValueError("Incorrect or no boundary conditions are defined. "
                         "Try using mesh.new_boundary_condition() before solving.")

    # --- Defining the Matrices ----------------------------------------------------------------------------------------
    # Stiffness (K_x, K_y) and Mass matrices (M)
    kx_matrix, ky_matrix, m_matrix, *_ = get_matrices(mesh)
    del _

    # K_xy = k_coef_x*G_x + k_coef_y*G_y
    k_matrix = (k_coef_x * kx_matrix) + (k_coef_y * ky_matrix)

    # Heat generation
    q_matrix = sparse.lil_matrix((mesh.size, 1))
    if isinstance(q, ComplexPointList):
        for _relative_index, _q in enumerate(q.indexes):
            q_matrix[_q] = q.values[_relative_index]

    if permanent_solution:
        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        #      A = K
        #      b = M / dt
        b_vector = sparse.lil_matrix(m_matrix.dot(q_matrix))
        apply_boundary_conditions(mesh, "space", k_matrix, b_vector)

        # --------------------------------- Solver ---------------------------------------------------------------------
        #    A.x = b   ->   x = solve(A, b)
        return linalg.spsolve(k_matrix.tocsc(), b_vector)

    else:
        # Use custom dt or obtain optimal dt based on size of elements
        dt = dt or mesh.default_dt

        #      A = M/dt + K
        a_matrix = m_matrix / dt + k_matrix

        # First frame of the solution (time = 0)
        t_vector = sparse.lil_matrix((mesh.size, 1))
        apply_initial_boundary_conditions(mesh, "time", t_vector)

        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        for _relative_index, _point in enumerate(mesh.boundary_conditions["space"].point_index_vector):
            if mesh.boundary_conditions["space"].type_of_condition_vector[_relative_index]:
                for _column_index in a_matrix.tocsr()[_point, :].indices:
                    a_matrix[_point, _column_index] = 0.
                a_matrix[_point, _point] = 1.
            else:
                q_matrix[_point, 0] += mesh.boundary_conditions["space"].values_vector[_relative_index]

        frames = [sparse_to_vector(t_vector)]
        for _frame_index in range(int(total_time / dt)):
            #      b = M.Q_i + (M/dt).(T_i^n-1)
            b_vector = sparse.lil_matrix(m_matrix.dot(q_matrix) + m_matrix.dot(t_vector.reshape(-1, 1)) / dt)

            # --------------------------------- Boundary conditions treatment ------------------------------------------
            #     b += C.C. Dirichlet/Neumann
            apply_boundary_conditions(mesh, "space", None, b_vector)

            #    A.x = b   ->   x = solve(A, b)
            t_vector = linalg.spsolve(a_matrix, b_vector)
            frames.append(t_vector)

        return frames


def solve_poiseuille(mesh, nu_coef=1.0, dt=None, total_time=1.0, save_each_frame=True):
    """
    Solves the mesh defined 2D Poiseuille equation problem:
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param dt: Value of time between frames.
    :param total_time: Length of time the calculation takes place.
    :param save_each_frame: True if every loop saves the current velocity values.
    :return: Velocity vectors and pressure values for each point in the mesh.
    """
    from scipy import sparse
    import scipy.sparse.linalg as linalg

    dt = dt or mesh.default_dt
    num_frames = int(total_time/dt)

    # --- Defining the Matrices ----------------------------------------------------------------------------------------
    #    K_x,       K_y,        M,       G_x,       G_y
    k_matrix, ky_matrix, m_matrix, gx_matrix, gy_matrix, = get_matrices(mesh)

    k_matrix = (k_matrix + ky_matrix)  # K_xy = K_x + K_y

    velocity_x_vector = sparse.lil_matrix((mesh.size, 1))
    velocity_y_vector = sparse.lil_matrix((mesh.size, 1))

    # --------------------------------- Set initial boundary conditions ------------------------------------------------
    apply_initial_boundary_conditions(mesh, "vel_x", velocity_x_vector)
    apply_initial_boundary_conditions(mesh, "vel_y", velocity_y_vector)

    # Wipe previous result set
    mesh.remove_previous_results(True)
    mesh.output_results(result_dictionary={"Velocity_X": sparse_to_vector(velocity_x_vector),
                                           "Velocity_Y": sparse_to_vector(velocity_y_vector)},
                        dt=dt, frame_num=0)

    # --------------------------------- Adding particles ---------------------------------------------------------------
    mesh.add_particle("A", (0.1 * (max(mesh.x) - min(mesh.x)), 0.5 * (max(mesh.y) - min(mesh.y))))
    particles = [Particle("B", (0.11 * (max(mesh.x) - min(mesh.x)), 0.8 * (max(mesh.y) - min(mesh.y))), color="b")]
    mesh.add_particle(list_of_particles=particles)
    mesh.add_particle("B", (0.11 * (max(mesh.x) - min(mesh.x)), 0.8 * (max(mesh.y) - min(mesh.y))), color="b")

    # --------------------------------- Solve Loop ---------------------------------------------------------------------
    for frame_num in range(1, num_frames + 1):
        print("Performing loop {0}".format(frame_num))

        # ------------------------ Acquire omega boundary condition ----------------------------------------------------
        #        M.w = (G_x.v_y) - (G_y.v_x)
        omega_vector = linalg.spsolve(m_matrix.tocsc(),
                                      (gx_matrix.dot(velocity_y_vector) - gy_matrix.dot(velocity_x_vector)))
        mesh.new_boundary_condition("omega", point_index=range(mesh.size), values=omega_vector, type_of_boundary=True)

        # ------------------------ Solve omega -------------------------------------------------------------------------
        #      A = M/dt + nu*K + (v.G)  -> v.G = G_x.(diagonal(v_x)) + G_y.(diagonal(v_y))
        a_matrix = sparse.lil_matrix(m_matrix / dt + nu_coef * k_matrix +
                                     (sparse_to_vector(velocity_x_vector) * gx_matrix +
                                      sparse_to_vector(velocity_y_vector) * gy_matrix))

        #      b = (M/dt).(omega^n-1)
        b_vector = sparse.lil_matrix((m_matrix / dt).dot(omega_vector)).T

        # Applying b.c.
        apply_boundary_conditions(mesh, "omega", a_matrix, b_vector)

        #        A.x = b   ->   x = solve(A, b)
        omega_vector = sparse.lil_matrix(linalg.spsolve(a_matrix.tocsc(), b_vector)).T

        # ------------------------ Solve psi ---------------------------------------------------------------------------
        #      A = K
        a_matrix = sparse.lil_matrix(k_matrix, copy=True)

        #      b = M.w
        b_vector = sparse.lil_matrix(m_matrix.dot(omega_vector), copy=True)

        # Applying b.c.
        apply_boundary_conditions(mesh, "psi", a_matrix, b_vector)

        #    A.x = b   ->   x = solve(A, b)
        psi_vector = sparse.lil_matrix(linalg.spsolve(a_matrix.tocsr(), b_vector)).T

        # ------------------------ Solve velocities --------------------------------------------------------------------
        #           M.v_x = G_y.psi
        velocity_x_vector = sparse.lil_matrix(linalg.spsolve(m_matrix.tocsc(), gy_matrix.dot(psi_vector))).T
        #           M.v_y = -G_x.psi
        velocity_y_vector = sparse.lil_matrix(linalg.spsolve(m_matrix.tocsc(), -gx_matrix.dot(psi_vector))).T

        # Applying b.c.
        apply_initial_boundary_conditions(mesh, "vel_x", velocity_x_vector)
        apply_initial_boundary_conditions(mesh, "vel_y", velocity_y_vector)

        # Show particle progress
        if False:  # TODO: REMOVE CHECK AND CREATE GIF
            if all([mesh.contains_particle(particle) for particle in mesh.particles]):
                mesh.show_geometry()
            else:
                print("There are no visible particles in the mesh domain.")

        # Move particles
        mesh.move_particles((sparse_to_vector(velocity_x_vector), sparse_to_vector(velocity_y_vector)), dt=dt)

        # Saving frames
        mesh.output_results(result_dictionary={"Velocity_X": sparse_to_vector(velocity_x_vector),
                                               "Velocity_Y": sparse_to_vector(velocity_y_vector)},
                            dt=dt, frame_num=frame_num)

    return sparse_to_vector(velocity_x_vector), sparse_to_vector(velocity_y_vector)


# -- Classes -----------------------------------------------------------------------------------------------------------

class Particle:
    """
    Defines a moving particle.
    Particles are defined as having a spherical shape of constant diameter.
    """
    position_history = None
    velocity_x = 0.
    velocity_y = 0.

    def __init__(self, name, position=(0., 0.), density=1., diameter=0.1, color="r"):
        """
        Particle class constructor.
        :param name: Name of the particle.
        :param position: Positional argument, (x, y) coordinates.
        :param density: Particle density.
        :param diameter: Diameter of particle (approximated as a sphere).
        :param color: Custom color for particle, default is "r" (red).
        """
        check_method_call(name)

        self.name = name
        self.pos_x = position[0]
        self.pos_y = position[1]
        self.position_history = [(self.pos_x, self.pos_y)]

        self.density = density
        self.diameter = diameter

        self.color = color

    def set_velocities(self, velocities, vel_x=0., vel_y=0.):
        """
        Determines the velocity of the particle for each axis.
        :param velocities: The list that contains each velocity value in the form (vel_x,vel_y)
        :param vel_x: The particle velocity in the x axis.
        :param vel_y: The particle velocity in the y axis.
        """
        if isinstance(velocities, (list, tuple)):
            self.velocity_x = velocities[0]
            self.velocity_y = velocities[1]
        else:
            self.velocity_x = vel_x
            self.velocity_y = vel_y

    def mark_new_position(self, dt):
        """
        Assigns a new position for the particle and tracks it's location history.
        :param dt: The time difference between frames.
        """
        self.pos_x += self.velocity_x * dt
        self.pos_y += self.velocity_y * dt
        self.position_history.append((self.pos_x, self.pos_y))


class BoundaryConditions:
    """
    Boundary conditions of the simulation.
    """
    point_index_vector = None
    values_vector = None
    type_of_condition_vector = None

    def set_new_boundary_condition(self, *vect_argm, point_index=range(0), values=0., type_of_boundary=True):
        """
        Sets the boundary condition for the mesh.
        :param vect_argm: Vector(s) of boundary conditions.
        :param point_index: Vector of order of points.
        :param values: Value or vector of values of the condition in each point.
        :param type_of_boundary: Value or vector of values, defined: True for Dirichlet and False for Neumann.
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


class Mesh:
    """
    Mesh element to be used in the calculations
    """
    x, y, ien, delauney_surfaces = None, None, None, None
    particles = None
    name = "default"
    size = 0
    default_dt = 0.

    def __init__(self, name="untitled", points=None):
        """
        Class constructor, initializes geometry.
        :param name: Mesh's main name.
        """
        import os
        from shutil import copy
        self.name = name
        self.particles = []

        if isinstance(points, list):
            self.import_point_structure(points=points)
        else:
            self.import_point_structure(import_mesh_file=self.name)

        self.boundary_conditions = {}

        try:
            os.chdir("./results/{0}/".format(self.name))
        except FileNotFoundError:
            os.mkdir("./results/{0}".format(self.name))
            os.chdir("./results/{0}/".format(self.name))

        if points is not None:
            copy("../{0}.msh".format(self.name), "./")

    def import_point_structure(self, *args, points=None, light_version=True, import_mesh_file=""):
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

        self.x, self.y, self.ien, self.delauney_surfaces = surface
        self.size = len(self.x)
        self.default_dt = get_dt(self)

    def new_boundary_condition(self, name, point_index=range(0), values=0., type_of_boundary=True):
        """
        Creates a new entry in the boundary conditions dictionary.
        :param name: Name of the property that has these boundary conditions.
        :param point_index: Vector of order of points.
        :param values: Value or vector of values of the condition in each point.
        :param type_of_boundary: Value or vector of values, defined: True for Dirichlet and False for Neumann.
        """
        self.boundary_conditions[name] = BoundaryConditions()
        self.boundary_conditions[name].set_new_boundary_condition(point_index=point_index, values=values,
                                                                  type_of_boundary=type_of_boundary)

    def add_particle(self, name=None, position=None, density=1., diameter=0.1, color="r", list_of_particles=None):
        """
        Associates a new particle with the mesh.
        :param name: Name of the particle. Each particle defined in a mesh must have different names.
        :param position: Positional argument, (x, y) coordinates.
        :param density: Particle density.
        :param diameter: Diameter of particle (approximated as a sphere).
        :param color: Custom color for particle, default is "r" (red).
        :param list_of_particles: Adds a list of predefined particle objects.
        """
        if list_of_particles:
            list_of_particles = check_list_or_object(list_of_particles, Particle)
            [self.particles.append(particle) for particle in list_of_particles if particle.name not in
             [obj.name for obj in self.particles]]

        else:
            check_method_call(name)
            if name in [obj.name for obj in self.particles]:
                print("\nThere is a particle already with that name. "
                      "Each particle defined inside a mesh must have a different name.\n")

            else:
                check_method_call(position)

                self.particles.append(Particle(name, position, density, diameter, color))

    def move_particles(self, velocity=None, velocity_vector_x=None, velocity_vector_y=None, dt=None):
        """
        Method that moves all particles currently inside the mesh domain.
        :param velocity: List or Tuple that contains the vector of velocity for each point in the mesh.
        :param velocity_vector_x: Vector of velocity in the x axis.
        :param velocity_vector_y: Vector of velocity in the y axis.
        :param dt: The time difference between frames.
        """
        check_method_call(dt)

        if isinstance(velocity, (list, tuple)) and len(velocity) == 2:
            velocity_vector_x = velocity[0]
            velocity_vector_y = velocity[1]

        check_method_call(velocity_vector_x, velocity_vector_y)

        for particle in [_particle for _particle in self.particles if self.contains_particle(_particle)]:
            # Determine which element contains the particle
            element = self.ien[self.delauney_surfaces.find_simplex((particle.pos_x, particle.pos_y))]

            total_area = get_area(self.x[element], self.y[element])
            # Interpolate velocity value for particle coordinates

            # Adds the velocity component from each point, interpolated by the proportion of the area opposite it.
            velocity_x = 0.
            velocity_y = 0.
            for point in element:
                _element = element.tolist()
                _element.remove(point)
                component = abs(get_area([particle.pos_x] + list(self.x[_element]),
                                         [particle.pos_y] + list(self.y[_element])) / total_area)
                velocity_x += component * velocity_vector_x[point]
                velocity_y += component * velocity_vector_y[point]

            particle.set_velocities((velocity_x, velocity_y))
            particle.mark_new_position(dt)

    def contains_particle(self, particle):
        """
        Checks if a particle is inside the mesh geometric domain.
        :param particle: Particle object that contains it's current position.
        :return: True if it is inside, otherwise False.
        """
        if not ((self.x.min() < particle.pos_x) and (particle.pos_x < self.x.max()) and
                (self.y.min() < particle.pos_y) and (particle.pos_y < self.y.max())):
            return False
        else:
            return True

    def remove_previous_results(self, return_to_previous_directory=False):
        """
        Removes previous results to prevent incorrect data analysis.
        :param return_to_previous_directory: Returns to previous directory.
        """
        import os
        import fnmatch

        # ---------------- Deleting previous results -------------------------------------------------------------------
        try:
            os.chdir("./paraview_results/")
        except FileNotFoundError:
            if not fnmatch.fnmatch(os.getcwd(), "*/paraview_results"):
                os.mkdir("paraview_results/")
                os.chdir("./paraview_results/")

        list(map(os.remove,
                 [file for file in os.listdir(".") if fnmatch.fnmatch(file, "{0}_results_*.vtk".format(self.name))]))

        if return_to_previous_directory:
            os.chdir("..")

    def output_results(self, result_vector=None, extension="VTK", dt=0., frame_num=None, data_names="Temperature",
                       result_dictionary=None):
        """
        Export result to .vtk or .csv file
        :param result_vector: The vector of the value in each point
        :param extension: File extension
        :param dt: Value of time between frames
        :param frame_num: Save one frame at a time (the number of that frame).
        :param data_names: Data name for display.
        :param result_dictionary: A dictionary that contains the names for each data and its values.
        """
        import os
        import fnmatch

        # ------------------ Contingency -------------------------------------------------------------------------------
        check_method_call(result_vector, result_dictionary)

        if not result_dictionary:
            result_dictionary = {}
            if isinstance(data_names, list):
                for index, data_name in enumerate(data_names):
                    result_dictionary[data_name] = result_vector[index]
            else:
                result_dictionary[data_names] = result_vector

        number_frames = 0
        for _name, _vector in result_dictionary.items():
            if " " in _name:
                raise AttributeError("There can't be any spaces in the property names!")

            if dt and frame_num is None:
                if len(_vector[0]) != self.size:
                    raise ValueError("Incorrect size for result _vector.")
                number_frames = len(_vector)

            elif len(_vector) != self.size:
                raise ValueError("Incorrect size for result vector.")

        # ----------------------- Change to result directory -----------------------------------------------------------
        try:
            os.chdir("./paraview_results/")
        except FileNotFoundError:
            if not fnmatch.fnmatch(os.getcwd(), "*/paraview_results"):
                os.mkdir("paraview_results/")
                os.chdir("./paraview_results/")

        number_elements = len(self.ien)

        if extension == "CSV":
            # ------------------------- Saving results to CSV file -----------------------------------------------------
            # TODO: FIX WITH NEW DICTIONARY APPLICATION
            with open("{0}_results.csv".format(self.name), "w") as arq:
                arq.write("Points:0, Points:1, Points:2, {0}\n".format(data_names))
                for i in range(self.size):
                    arq.write("{0},{1},{2},{3}\n".format(self.x[i], self.y[i], 0, result_vector[i]))
            return os.chdir("..")

        if extension == "VTK":
            size = self.size

            def write_scalars(file, property_name, vector):
                file.write("\nSCALARS {0} float 1\n".format(property_name))
                file.write("\nLOOKUP_TABLE {0}\n".format(property_name))
                for _i in range(size):
                    file.write("{}\n".format(vector[_i]))

            def write_header_and_cells(file):
                # ------------------------------------ Header ----------------------------------------------------------
                file.write("# vtk DataFile Version 3.0\n{0}\n{1}\n\nDATASET {2}\n".format("LucasCarvalhoTCC Results",
                                                                                          "ASCII", "POLYDATA"))
                if dt:
                    file.write("FIELD FieldData 1\nTIME 1 1 double\n{}\n".format(dt))
                # ------------------------------------ Points coordinates ----------------------------------------------
                file.write("\nPOINTS {0} {1}\n".format(size, "float"))
                for _i in range(size):
                    file.write("{0} {1} 0.0\n".format(self.x[_i], self.y[_i]))
                # --------------------------------------- Cells --------------------------------------------------------
                file.write("\nPOLYGONS {0} {1}\n".format(number_elements, number_elements * 4))
                for _i in range(number_elements):
                    file.write("{0} {1} {2} {3}\n".format(3, self.ien[_i][0], self.ien[_i][1], self.ien[_i][2]))

            if dt:
                if frame_num is not None:
                    # ----- Saving a single result (frame) to VTK file -------------------------------------------------
                    with open("{0}_results_{1}.vtk".format(self.name, frame_num), "w") as arq:
                        # --------------------------------- Header and cells -------------------------------------------
                        write_header_and_cells(arq)

                        # ------------------------------------ Data in each point --------------------------------------
                        arq.write("\nPOINT_DATA {0}\n".format(size))
                        [write_scalars(arq, name, vector) for name, vector in result_dictionary.items()]

                else:
                    self.remove_previous_results()
                    # ----- Saving multiple results to VTK files -------------------------------------------------------
                    for j in range(number_frames):
                        with open("{0}_results_{1}.vtk".format(self.name, j), "w") as arq:
                            # ----------------------------- Header and cells -------------------------------------------
                            write_header_and_cells(arq)

                            # -------------------------------- Data in each point --------------------------------------
                            arq.write("\nPOINT_DATA {0}\n".format(size))
                            [write_scalars(arq, name, vector[j]) for name, vector in result_dictionary.items()]

            else:
                # ------------------------- Saving results to VTK file -------------------------------------------------
                with open("{0}_results.vtk".format(self.name), "w") as arq:
                    # --------------------------------- Header and cells -----------------------------------------------
                    write_header_and_cells(arq)

                    # ----------------------------------- Data in each point--------------------------------------------
                    arq.write("\nPOINT_DATA {0}\n".format(size))
                    [write_scalars(arq, name, vector) for name, vector in result_dictionary.items()]

            return os.chdir("..")

        raise NameError("Format not available. Try VTK or CSV.")

    def show_geometry(self, particles=None, names=False, rainbow=False, save=False):
        """
        Display mesh geometry on screen using matplotlib.
        :param particles: Particle or list of particles to be displayed inside the geometry.
        :param names: Show the index of each point next to it.
        :param rainbow: Color in the edge of each element in a different color.
        :param save: Save generated image.
        :return: Display image.
        """
        import numpy as np
        import matplotlib.pyplot as plt

        particles = check_list_or_object(particles, Particle)
        particles.extend(self.particles)

        # Draw mesh points
        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        # Draw particles
        [plt.scatter(particle.pos_x, particle.pos_y, c=particle.color) for particle in particles
         if self.contains_particle(particle)]

        # Sets plot styling
        style_plot(self.x, self.y)

        # Display id's of points and names of objects
        if names:
            for _index in range(self.size):
                plt.gca().annotate(_index, (self.x[_index], self.y[_index]))

        [plt.gca().annotate(particle.name, (particle.pos_x, particle.pos_y), color=particle.color, fontsize=15)
         for particle in particles if self.contains_particle(particle)]

        # Displays elements as different colors to help distinguish each one
        if rainbow:
            plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.hsv(np.linspace(0.0, 1.0, len(self.ien)))))

            for element in self.ien:
                plot_coordinates = ((self.x[element[0]], self.x[element[1]], self.x[element[2]], self.x[element[0]]),
                                    (self.y[element[0]], self.y[element[1]], self.y[element[2]], self.y[element[0]]))

                plt.plot(plot_coordinates[0], plot_coordinates[1])

        else:
            plt.triplot(self.x, self.y, triangles=self.ien[0])

        if save:
            plt.savefig("{0}_mesh".format(self.name))

        return plt.show()

    def show_3d_solution(self, solution_vector):
        """
        Display 3D solution of the mesh geometry.
        :param solution_vector: Vector that contains the value of the solution for each point in the mesh.
        :return: Display image.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        check_method_call(solution_vector)
        try:
            if len(solution_vector) != self.size:
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Solution must be a vector")

        self.output_results(result_dictionary={"Temperature": solution_vector})

        fig = plt.gcf()
        axes = Axes3D(fig)
        surf = axes.plot_trisurf(self.x, self.y, solution_vector, cmap="jet")
        axes.view_init(90, 270)
        fig.colorbar(surf, shrink=0.4, aspect=9)

        plt.savefig("{0}_permanent_results".format(self.name))

        return plt.show()

    def show_animated_3d_solution(self, frames_vector, dt=0.):
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
            self.output_results(result_dictionary={"Temperature": frames_vector}, dt=dt)
            animation.save("{0}_transient_results.gif".format(self.name), dpi=80, writer='imagemagick')

        return plt.show()

    def show_velocity_quiver(self, velocity_x, velocity_y):
        """

        :param velocity_x:
        :param velocity_y:
        :return:
        """
        import matplotlib.pyplot as plt
        # TODO: CONTINUE METHOD IMPLEMENTATION

        check_method_call(velocity_x, velocity_y)

        try:
            if len(velocity_x) != self.size or len(velocity_y) != self.size:
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Solution must be a vector")

        self.output_results(result_dictionary={"Velocity_X": velocity_x, "Velocity_Y": velocity_y})

        fig, axes = plt.subplots()
        axes.quiver(self.x, self.y, velocity_x, velocity_y)

        # plt.savefig("{0}_permanent_results".format(self.name))

        return plt.show()

    def show_particle_movement(self, save=False):
        """
        Displays an animated image of the particles trajectories.
        :return: Display image.
        """
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        # Draw mesh points
        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        number_of_frames = max([len(self.particles[i].position_history) for i in range(len(self.particles))])
        list_of_dots = {}

        for particle in self.particles:
            list_of_dots[particle] = plt.scatter(0, 0, s=100, c=particle.color)
            # TODO: CHANGE SIZE TO BE DETERMINED BY DIAMETER OF PARTICLE

        def update(_frame):
            # Draw particles
            for _particle, dot in list_of_dots.items():
                if _frame < len(_particle.position_history):
                    dot.set_offsets((_particle.position_history[_frame][0], _particle.position_history[_frame][1]))

            return

        animation = FuncAnimation(plt.gcf(), update, frames=number_of_frames, interval=100, save_count=False)

        if save:
            animation.save("{0}_particle_movement.gif".format(self.name), dpi=80, writer='imagemagick')

        return plt.show()


class ComplexPointList:
    """
    Class that defines a list of points with an index value and a property value each.
    """
    indexes = None
    values = None

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
            self.values = []
            for i in range(len(_indexes)):
                self.values.append(_values)

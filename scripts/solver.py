import scipy.sparse.linalg as linalg
from scipy import sparse

from scripts.complexPointList import *
from scripts.mesh import *


def get_matrices(mesh: Mesh):
    """
    Function that generates the algebraic components of the solution method
    :param mesh: Mesh object to be used to generate the matrices
    :return:
    """
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

        area = util.get_area(x, y)

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


def apply_boundary_conditions(mesh: Mesh, boundary_name: str, matrix_a, vector_b, const=1.0):
    """
    Performs the evaluation of boundary conditions and applies the changes to the main matrix and vector.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param boundary_name: The boundary name as defined for each problem.
    :param matrix_a: The coefficients matrix A in the linear solution method [A*x = b].
    :param vector_b: The results vector b in the linear solution method [A*x = b].
    :param const: Optional correction constant.
    """
    for _rel_index, _column_index in enumerate(mesh.boundary_conditions[boundary_name].point_index_vector):
        if mesh.boundary_conditions[boundary_name].type_of_condition_vector[_rel_index]:
            # Dirichlet Treatment
            if matrix_a is not None:
                # Skip applying matrix a values for loops
                for _line_index in matrix_a.tocsc()[:, _column_index].indices:
                    vector_b[_line_index, 0] -= (matrix_a[_line_index, _column_index] *
                                                 mesh.boundary_conditions[boundary_name].values_vector[_rel_index])
                    matrix_a[_line_index, _column_index] = 0.
                    matrix_a[_column_index, _line_index] = 0.

                matrix_a[_column_index, _column_index] = 1.

            vector_b[_column_index, 0] = mesh.boundary_conditions[boundary_name].values_vector[_rel_index]
        else:
            # Neumann Treatment
            # TODO: FIX NEUMANN BOUNDARY TREATMENT
            vector_b[_column_index, 0] += 0. if mesh.boundary_conditions[boundary_name].values_vector[_rel_index] == 0\
                else (mesh.boundary_conditions[boundary_name].values_vector[_rel_index] * const * mesh.mean_side_length)


def apply_initial_boundary_conditions(mesh: Mesh, boundary_name, vector_v):
    """
    Performs the evaluation of boundary conditions and applies initial values to the vector.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param vector_v: The vector v of initial values for the transient solution.
    :param boundary_name: The boundary name as defined for each problem.
    """
    for _relative_index, _point in enumerate(mesh.boundary_conditions[boundary_name].point_index_vector):
        if mesh.boundary_conditions[boundary_name].type_of_condition_vector[_relative_index]:
            vector_v[_point] = mesh.boundary_conditions[boundary_name].values_vector[_relative_index]


def solve_poisson(mesh: Mesh, permanent_solution: bool = True, k_coef: float = None, k_coef_x: float = 1.0,
                  k_coef_y: float = 1.0, q: [float, ComplexPointList] = None, dt: float = None,
                  total_time: float = None, stop_rule: float = None, return_history: bool = False):
    """
    Solves the mesh defined 2D Poisson equation problem:
        DT = -∇(k*∇T) + Q   ->   (M + K).T_i^n =  M.T_i^n-1 + M.Q_i
        dt                       dt              dt
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param permanent_solution: Parameter that defines if the solution will be calculated for the transient (True) or
                                permanent (False) problem.
    :param k_coef: Thermal conductivity coefficient both axis [W/(m*K)].
    :param k_coef_x: Thermal conductivity coefficient for x axis [W/(m*K)].
    :param k_coef_y: Thermal conductivity coefficient for y axis [W/(m*K)].
    :param q: Heat generation for each point.
    :param dt: Value of time between frames [s].
    :param total_time: Length of time the calculation takes place (only necessary for transient solutions) [s].
    :param stop_rule: Precision used to detect if method can be stopped early.
    :param return_history: Flag used to check if return is resulting array of property values in the mesh,
                           or a MeshPropertyStates object that contains the information of values in various timestamps.
    :return: Temperature value for each point in the mesh [K or °C].
    """
    k_coef_x = k_coef or k_coef_x
    k_coef_y = k_coef or k_coef_y
    correction_coef = (k_coef_x + k_coef_y)/(2.*(0.909759 + 15.4696*mesh.mean_area))

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
    if isinstance(q, (float, int)):
        q_matrix = q_matrix.toarray() + q

    if permanent_solution:
        # --------------------------------- Boundary conditions treatment ----------------------------------------------
        #      A = K
        #      b = M / dt
        b_vector = sparse.lil_matrix(m_matrix.dot(q_matrix))
        apply_boundary_conditions(mesh, "space", k_matrix, b_vector, correction_coef)

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
            # Dirichlet Treatment
            if mesh.boundary_conditions["space"].type_of_condition_vector[_relative_index]:
                for _column_index in a_matrix.tocsr()[_point, :].indices:
                    a_matrix[_point, _column_index] = 0.
                a_matrix[_point, _point] = 1.
            else:
                # Neumann Treatment
                q_matrix[_point, 0] += 0. if mesh.boundary_conditions["space"].values_vector[_relative_index] == 0\
                    else (mesh.boundary_conditions["space"].values_vector[_relative_index] * correction_coef *
                          mesh.mean_side_length)

        t_vector = util.sparse_to_vector(t_vector)
        states = MeshPropertyStates(t_vector)
        for time in np.arange(dt, total_time or dt*1e3, dt):
            #      b = M.Q_i + (M/dt).(T_i^n-1)
            b_vector = sparse.lil_matrix(m_matrix.dot(q_matrix) + m_matrix.dot(t_vector.reshape(-1, 1)) / dt)

            # --------------------------------- Boundary conditions treatment ------------------------------------------
            #     b += C.C. Dirichlet/Neumann
            apply_boundary_conditions(mesh, "space", None, b_vector, correction_coef)

            #    A.x = b   ->   x = solve(A, b)
            t_vector = linalg.spsolve(a_matrix, b_vector)
            if stop_rule is not None and abs(np.mean(t_vector - states[-1])) < stop_rule:
                break
            states.append(t_vector, time)

        return states if return_history else states[-1]


def solve_velocity_field(mesh: Mesh, dt: float = None, total_time: float = 1.0, reynolds: float = None,
                         save_each_frame: bool = False, stop_rule: float = None, return_history: bool = True):
    """
    Solves the mesh defined 2D current-vorticity equation problem:
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param dt: Value of time between frames [s].
    :param total_time: Length of time the calculation takes place [s].
    :param reynolds: Option to provide the value of Reynolds Number [1].
    :param save_each_frame: True if every loop saves the current velocity values.
    :param stop_rule: Precision used to detect if method can be stopped early.
    :param return_history: Flag used to check if return is resulting array of property values in the mesh,
                           or a MeshPropertyStates object that contains the information of values in various timestamps.
    :return: Velocity vectors and pressure values for each point in the mesh.
    """
    dt = dt or mesh.default_dt

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
    mesh.output_results(result_dictionary={"Velocity_X": util.sparse_to_vector(velocity_x_vector),
                                           "Velocity_Y": util.sparse_to_vector(velocity_y_vector)},
                        dt=dt, frame_num=0)

    # Show initial particle position
    # if save_each_frame:
    #     mesh.show_geometry()
    #     mesh.save_frame(0)

    velocity_x_states = MeshPropertyStates(util.sparse_to_vector(velocity_x_vector))
    velocity_y_states = MeshPropertyStates(util.sparse_to_vector(velocity_y_vector))
    # --------------------------------- Solve Loop ---------------------------------------------------------------------
    for frame_index, time in enumerate(np.arange(dt, total_time, dt)):
        print("\rSolving velocity {0:.2f}%".format(100 * time / total_time), end="")

        # Defining Reynolds number
        re = reynolds or mesh.density * max(velocity_x_vector)[0, 0] * mesh.length_y / mesh.viscosity
        if re > 100.:
            print("Reynolds value is beyond defined maximum value: Re = {0} > 100".format(int(re)))

        # ------------------------ Acquire omega boundary condition ----------------------------------------------------
        #        M.w = (G_x.v_y) - (G_y.v_x)
        omega_vector = linalg.spsolve(m_matrix.tocsc(),
                                      (gx_matrix.dot(velocity_y_vector) - gy_matrix.dot(velocity_x_vector)))
        mesh.new_boundary_condition("omega", point_index=range(mesh.size), values=omega_vector, type_of_boundary=True)

        # ------------------------ Solve omega -------------------------------------------------------------------------
        #      A = M/dt + nu*K + (v.G)  -> v.G = G_x.(diagonal(v_x)) + G_y.(diagonal(v_y))
        a_matrix = sparse.lil_matrix(m_matrix / dt + (1./re) * k_matrix +
                                     (util.sparse_to_vector(velocity_x_vector) * gx_matrix +
                                      util.sparse_to_vector(velocity_y_vector) * gy_matrix))

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

        if stop_rule is not None and (abs(np.mean(velocity_x_vector.toarray() - velocity_x_states.last)) < stop_rule and
                                      abs(np.mean(velocity_y_vector.toarray() - velocity_y_states.last)) < stop_rule):
            break

        # Saving frames
        if return_history:
            velocity_x_states.append(util.sparse_to_vector(velocity_x_vector), time)
            velocity_y_states.append(util.sparse_to_vector(velocity_y_vector), time)
        if save_each_frame:
            # mesh.save_frame(time)
            mesh.output_results(result_dictionary={"Velocity": {"x": util.sparse_to_vector(velocity_x_vector),
                                                                "y": util.sparse_to_vector(velocity_y_vector)}},
                                dt=dt, frame_num=frame_index)

    print("\rSolving velocity done!")
    return [velocity_x_states, velocity_y_states] if return_history else [velocity_x_vector, velocity_y_vector]


def move_particles(mesh: Mesh, velocity: (list, tuple) = None, velocity_x: [list, np.ndarray] = None,
                   velocity_y: [list, np.ndarray] = None, dt: float = None, single_force: str = None,
                   acceleration_x: [list, np.ndarray] = None, acceleration_y: [list, np.ndarray] = None,
                   no_gravity: bool = False):
    """
    Method that moves all particles currently inside the mesh domain.
    :param mesh: The Mesh object that defines the geometry of the problem and the boundary conditions associated.
    :param velocity: List or Tuple that contains the vector of velocity for each point in the mesh [m/s].
    :param velocity_x: Vector of velocity in the x axis [m/s].
    :param velocity_y: Vector of velocity in the y axis [m/s].
    :param acceleration_x: Vector of acceleration in the fluid field in the x axis [m/s²].
    :param acceleration_y: Vector of acceleration in the fluid field in the y axis [m/s²].
    :param dt: The time difference between frames [s].
    :param single_force: Parameter used for testing single forces one at a time.
    :param no_gravity: Flag for disabling gravity influence on particles, used for simulations viewed from above.
    """
    # Contingency
    dt = dt or mesh.default_dt
    dt = min([dt] + [particle.max_dt(mesh.viscosity) for particle in mesh.particles])

    if isinstance(velocity, (list, tuple)) and len(velocity) == 2:
        velocity_x = velocity[0]
        velocity_y = velocity[1]

    util.check_method_call(velocity_x, velocity_y)

    if (acceleration_x is not None) or (acceleration_y is not None):
        if acceleration_x is None:
            acceleration_x = np.zeros(mesh.size)
        if acceleration_y is None:
            acceleration_y = np.zeros(mesh.size)

    # Applying forces to each particle if it is still able.
    for particle in [_particle for _particle in mesh.particles if mesh.contains_particle(_particle)]:
        forces = dict()

        # ----------------- Obtain velocity ----------------------------------------------------------------------------
        fluid_velocity = np.array(mesh.get_interpolated_value((particle.pos_x, particle.pos_y),
                                                              velocity_x, velocity_y))
        relative_vel = fluid_velocity - particle.velocity
        # relative_vel_norm = np.sqrt(relative_vel.dot(relative_vel))

        particle.reynolds = (mesh.density * max(relative_vel) * particle.diameter / mesh.viscosity)

        # ----------------- Gravitational Force ------------------------------------------------------------------------
        if not no_gravity:
            forces["gravitational"] = (0., -9.80665 * (particle.density - mesh.density)*particle.volume)

        # ----------------- Drag Force ---------------------------------------------------------------------------------
        # if particle.reynolds > 1:
        #     print("Reynolds number beyond usable definitions.")
        forces["drag"] = 3 * np.pi * mesh.viscosity * particle.diameter * relative_vel

        # -------------------- Lift Force ------------------------------------------------------------------------------
        dn = np.array([-fluid_velocity[1], fluid_velocity[0]])
        dn = dn/np.sqrt(dn.dot(dn)) * particle.radius
        dv_dr = (np.array(mesh.get_interpolated_value((particle.pos_x + dn[0], particle.pos_y + dn[1]),
                                                      velocity_y, velocity_x, True)) -
                 np.array(mesh.get_interpolated_value((particle.pos_x - dn[0], particle.pos_y - dn[1]),
                                                      velocity_y, velocity_x, True)))
        forces["lift"] = (np.nan_to_num(dn/abs(dn)) * 1.61 * mesh.viscosity * particle.diameter * relative_vel *
                          np.sqrt(particle.diameter * mesh.density / mesh.viscosity * abs(dv_dr)))

        # -------------------- Added Mass Force ------------------------------------------------------------------------
        acceleration = (0, 0) if (acceleration_x is None) or (acceleration_y is None) else \
            np.array(mesh.get_interpolated_value((particle.pos_x, particle.pos_y), acceleration_x, acceleration_y))
        forces["added_mass"] = mesh.density/2. * particle.volume * (acceleration -
                                                                    (particle.velocity - particle.last_velocity)/dt)

        particle.apply_forces(forces, mesh, dt, single_force)

import TccLib

import numpy as np

# Import gmsh created mesh
mesh = TccLib.Mesh("Poisson_Dirichlet")

# Show mesh geometry
mesh.show_geometry(names=True)

# Define parameters
T_0 = 0.
T_L = 0.
k = 5
Q = 40

# Define analytic comparison expression
analytic_expression = lambda x: Q/(2.*k) * (-x**2 + mesh.length_y*x) + (T_L-T_0)/mesh.length_y * x + T_0

# Define boundary conditions
boundary_conditions_values = {
    "north": T_L,
    "south": T_0,
}

boundary_conditions_types = {
    "north": True,
    "south": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values,
                                                                   boundary_conditions_types)

# Set Boundary Conditions
mesh.new_boundary_condition("space", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)
mesh.new_boundary_condition("time", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)

# Solve for permanent solution
temperature_perm = TccLib.solve_poisson(mesh, permanent_solution=True, k_coef=k, q=Q)

# Show results in 3D graph
mesh.show_3d_solution(temperature_perm, view_from_above=False, axis_labels={"z": "T(°C)"})

# Define x vector of positions
x_vector = np.linspace(min(mesh.y), max(mesh.y), 100)

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = mesh.length_x/2.
y_vector = [mesh.get_interpolated_value([x_position, x], temperature_perm) for x in x_vector]

# Show comparison graph
TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Equação de Laplace Permanente",
                            x_label="Posição no Eixo de Comparação(m)", y_label="Valor da Temperatura (°C)",
                            save_file_as="{0}_permanent_validation".format(mesh.name))

# Solve for transient solution
temperature_trans = TccLib.solve_poisson(mesh, permanent_solution=False, k_coef=k, stop_criteria=1e-5, q=Q,
                                         return_history=True)

# Show results in 3D graph
mesh.show_animated_3d_solution(temperature_trans, dt=mesh.default_dt, axis_labels={"z": "T(°C)"})

# Get a small dictionary with each value state with timestamps as keys
small_dict = temperature_trans.reduced_dict_log(5)

# Find values of property in the mesh at a determined set position for every value in the vector of x at each timestamp
z_vectors = {"{:.4f}s".format(key): [mesh.get_interpolated_value([x_position, x], small_dict[key]) for x in x_vector]
             for key in small_dict}

# Show comparison graph
TccLib.util.show_comparison(x_vector, analytic_expression, z_vectors, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Equação de Laplace Transiente",
                            x_label="Posição no Eixo de Comparação(m)", y_label="Valor da Temperatura (°C)",
                            save_file_as="{0}_transient_validation".format(mesh.name))

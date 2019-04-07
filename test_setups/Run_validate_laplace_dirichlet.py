import TccLib

import numpy as np
import os

os.chdir("..")

# Import gmsh created mesh
mesh = TccLib.Mesh("Poisson")

# Show mesh geometry
mesh.show_geometry(names=True)

# Define boundary conditions parameters
boundary_conditions_values = {
    "north": lambda i: 1.0,
}

boundary_conditions_types = {
    "north": True,
    "south": True,
}
k_coef = 5

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values,
                                                                   boundary_conditions_types)

# Set Boundary Conditions
mesh.new_boundary_condition("space", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)
mesh.new_boundary_condition("time", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)

# Solve for permanent solution
temperature_perm = TccLib.solve_poisson(mesh, permanent_solution=True, k_coef=k_coef)

# Show results in 3D graph
mesh.show_3d_solution(temperature_perm, view_from_above=False)

# Define x vector of positions
x_vector = np.linspace(min(mesh.y), max(mesh.y), 100)

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = (max(mesh.x) - min(mesh.x))/2.
y_vector = [mesh.get_interpolated_value([x_position, x], temperature_perm) for x in x_vector]

# Show comparison graph
TccLib.util.show_comparison(x_vector, lambda x: x, y_vector, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Equação de Laplace Permanente",
                            x_label="Posição no Eixo X(m)", y_label="Valor da Temperatura (°C)")

# Solve for transient solution
temperature_trans = TccLib.solve_poisson(mesh, permanent_solution=False, k_coef=k_coef, stop_criteria=1e-5,
                                         return_history=True)

# Show results in 3D graph
mesh.show_animated_3d_solution(temperature_trans)

# Get a small dictionary with each value state with timestamps as keys
small_dict = temperature_trans.reduced_dict_log(5)

# Find values of property in the mesh at a determined set position for every value in the vector of x at each timestamp
z_vectors = {key: [mesh.get_interpolated_value([x_position, x], small_dict[key]) for x in x_vector]
             for key in small_dict}

# Show comparison graph
TccLib.util.show_comparison(x_vector, lambda x: x, z_vectors, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Equação de Laplace Transiente",
                            x_label="Posição no Eixo X(m)", y_label="Valor da Temperatura (°C)")

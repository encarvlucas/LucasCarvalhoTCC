import TccLib

import numpy as np

# Define boundary conditions and parameters
vel_top = 1.
vel_bot = 0.
dt = 1.0

# Set liquid parameters or declare liquid
density = 1e3
viscosity = 0.89e-3
liquid = "oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Couette", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=True)

# Define analytic comparison expression
analytic_expression = lambda y: (vel_top - vel_bot) * y / mesh.length_y + vel_bot

# -------------------------------------------------- PSI ---------------------------------------------------------------
# Define boundary conditions for psi
boundary_conditions_values_psi = {
    "north": lambda y: mesh.y[y],
    "east": lambda y: mesh.y[y],
    "south": lambda y: mesh.y[y],
    "west": lambda y: mesh.y[y],
}

boundary_conditions_types_psi = {
    "north": True,
    "east": True,
    "south": True,
    "west": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_psi,
                                                                   boundary_conditions_types_psi)

# Set Boundary Conditions for psi
mesh.new_boundary_condition("psi", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)

# -------------------------------------------------- VEL X -------------------------------------------------------------
# Define boundary conditions for vel_x
boundary_conditions_values_x = {
    "north": vel_top,
    "east": 0.,
    "south": vel_bot,
    "west": 0.,
}

boundary_conditions_types_x = {
    "north": True,
    "east": False,
    "south": True,
    "west": False,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_x,
                                                                   boundary_conditions_types_x)

# Set Boundary Conditions for vel_x
mesh.new_boundary_condition("vel_x", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)

# -------------------------------------------------- VEL Y -------------------------------------------------------------
# Define boundary conditions for vel_y
boundary_conditions_values_y = {
    "north": 0.,
    "east": 0.,
    "south": 0.,
    "west": 0.,
}

boundary_conditions_types_y = {
    "north": True,
    "east": False,
    "south": True,
    "west": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_y,
                                                                   boundary_conditions_types_y)

# Set Boundary Conditions for vel_y
mesh.new_boundary_condition("vel_y", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)

# Solve for FEM velocity field solution
velocity_x, velocity_y = TccLib.solve_velocity_field(mesh, dt=dt, total_time=30.)
TccLib.util.save(velocity_x, "vel_x")
# velocity_x = TccLib.util.load("vel_x")

# Show results in quiver plot
# mesh.show_velocity_quiver(velocity_x, velocity_y)

# Define x vector of positions
x_vector = np.linspace(min(mesh.y), max(mesh.y), 100)

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = (max(mesh.x) - min(mesh.x))/2.
y_vector = [mesh.get_interpolated_value([x_position, x], velocity_x) for x in x_vector]

# Show comparison graph
TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Equação de Laplace Permanente",
                            x_label="Posição no Eixo de Comparação(m)", y_label="Valor da Temperatura (°C)",
                            save_file_as="{0}_permanent_validation".format(mesh.name))

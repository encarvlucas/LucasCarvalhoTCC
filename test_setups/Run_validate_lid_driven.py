import TccLib

import numpy as np

# Define boundary conditions and parameters
vel_top = 1.
dt = 01.
total_time = 50.

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Lid_driven", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=True)

# Define analytic comparison expression
analytic_expression = lambda y: y

# -------------------------------------------------- PSI ---------------------------------------------------------------
# Define boundary conditions for psi
boundary_conditions_values_psi = {
    "all": 0,
}

boundary_conditions_types_psi = {
    "all": True,
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
    "south": 0.,
    "west": 0.,
}

boundary_conditions_types_x = {
    "all": True,
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
    "all": 0.,
}

boundary_conditions_types_y = {
    "all": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_y,
                                                                   boundary_conditions_types_y)

# Set Boundary Conditions for vel_y
mesh.new_boundary_condition("vel_y", point_index=xy_indices, values=xy_values,
                            type_of_boundary=xy_types)

# Solve for FEM velocity field solution
# velocity_x, velocity_y = TccLib.solve_velocity_field(mesh, dt=dt, total_time=total_time, save_each_frame=True,
#                                                      stop_criteria=1e-5)
# TccLib.util.save(velocity_x, "vel_x")
# TccLib.util.save(velocity_y, "vel_y")
velocity_x = TccLib.util.load("vel_x")
velocity_y = TccLib.util.load("vel_y")

# Show results in quiver plot
mesh.show_velocity_quiver(velocity_x, velocity_y)

# Define x vector of positions
x_vector = np.linspace(min(mesh.y), max(mesh.y), 100)
small_dict = velocity_x.reduced_dict_log(5)

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = mesh.length_x*0.6
y_vector = {"{:.4f}s".format(key): [mesh.get_interpolated_value([x_position, x], small_dict[key]) for x in x_vector]
            for key in small_dict}

# Show comparison graph
# TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
#                             analytic_label="Solução Analítica", title="Escoamento de Cavidade",
#                             x_label="Posição no Eixo de Comparação(m)", y_label="Valor da Velocidade (m/s)",
#                             save_file_as="{0}_validation".format(mesh.name))

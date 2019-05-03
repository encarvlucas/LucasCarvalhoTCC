import TccLib

import numpy as np

# Define boundary conditions and parameters
vel = 1.
dt = 0.1
total_time = 20.

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "super_oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Step_ref", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=True, save=True)

# Concave sides boundary points
# sides_vertical = np.array([
#     19, 20, 21, 22, 23, 24, 25, 26, 27,
#     43, 44, 45, 46, 47, 48, 49, 50, 51,
# ])
# sides_horizontal_right = np.array([3, 28, 29, 30, 31])
# sides_horizontal_left = np.array([8, 52, 53, 54, 55])
# Refined mesh
sides_vertical = np.array([
    24, 25, 26, 27, 28, 29, 30, 31, 32,
    60, 61, 62, 63, 64, 65, 66, 67, 68
])
sides_horizontal_right = np.array([3, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 4])
sides_horizontal_left = np.array([8, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 9])

sides = np.append(sides_vertical, np.append(sides_horizontal_right, sides_horizontal_left))

# -------------------------------------------------- PSI ---------------------------------------------------------------
# Define boundary conditions for psi
boundary_conditions_values_psi = {
    "north": lambda y: mesh.y[y] - mesh.length_y/2.,
    "east": lambda y: mesh.y[y],
    "south": lambda y: mesh.y[y],
    "west": lambda y: mesh.y[y] - mesh.length_y/2.,
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
mesh.new_boundary_condition("psi", point_index=xy_indices, values=xy_values, type_of_boundary=xy_types)

# Set boundary conditions for obstacle
mesh.new_boundary_condition("psi", point_index=sides_horizontal_left, values=0., type_of_boundary=True)
mesh.new_boundary_condition("psi", point_index=sides_horizontal_right, values=1., type_of_boundary=True)
mesh.new_boundary_condition("psi", point_index=sides_vertical, values=0., type_of_boundary=False)

# -------------------------------------------------- VEL X -------------------------------------------------------------
# Define boundary conditions for vel_x
boundary_conditions_values_x = {
    "north": 0.,
    "east": 0.,
    "south": 0.,
    "west": vel,
}

boundary_conditions_types_x = {
    "north": True,
    "east": False,
    "south": True,
    "west": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_x,
                                                                   boundary_conditions_types_x)

# Set Boundary Conditions for vel_x
mesh.new_boundary_condition("vel_x", point_index=xy_indices, values=xy_values, type_of_boundary=xy_types)

# Set boundary conditions for obstacle
mesh.new_boundary_condition("vel_x", point_index=sides, values=0., type_of_boundary=True)

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
mesh.new_boundary_condition("vel_y", point_index=xy_indices, values=xy_values, type_of_boundary=xy_types)

# Set boundary conditions for obstacle
mesh.new_boundary_condition("vel_y", point_index=sides, values=0., type_of_boundary=True)

# Solve for FEM velocity field solution
velocity_x, velocity_y = TccLib.solve_velocity_field(mesh, dt=dt, total_time=total_time, save_each_frame=True,
                                                     stop_criteria=1e-5)
TccLib.util.save(velocity_x, "vel_x")
TccLib.util.save(velocity_y, "vel_y")

# velocity_x = TccLib.util.load("vel_x")
# velocity_y = TccLib.util.load("vel_y")

# Show results in quiver plot
mesh.show_velocity_quiver(velocity_x, velocity_y)

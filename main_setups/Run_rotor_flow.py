import TccLib

import numpy as np

# Define boundary conditions and parameters
vel = 1.
dt = 0.5
total_time = 5.

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "super_oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Rotor_ref", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=True, save=True)

# Rotor boundary points
# north = np.array([42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61])
# east = np.array([2, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 3])
# south = np.array([16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
# west = np.array([1, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 0])
# Refined mesh
north = np.array([55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                  80, 81, 82])
east = np.array([2, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 3])
south = np.array([20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38])
west = np.array([1, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 0])

# -------------------------------------------------- PSI ---------------------------------------------------------------
# Define boundary conditions for psi
boundary_conditions_values_psi = {
    "north": 1.,
    "east": lambda x: np.where(east == x)[0]/(len(east)-1),
    "south": 0.,
    "west": lambda x: np.where(west == x)[0]/(len(west)-1),
}

boundary_conditions_types_psi = {
    "north": True,
    "east": True,
    "south": True,
    "west": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_psi,
                                                                   boundary_conditions_types_psi, north=north,
                                                                   east=east, south=south, west=west)

# Set Boundary Conditions for psi
mesh.new_boundary_condition("psi", point_index=xy_indices, values=xy_values, type_of_boundary=xy_types)

# -------------------------------------------------- VEL X -------------------------------------------------------------
# Define boundary conditions for vel_x
boundary_conditions_values_x = {
    "north": 0.,
    "east": 0.,
    "south": 0.,
    "west": 0.,
}

boundary_conditions_types_x = {
    "north": True,
    "east": False,
    "south": True,
    "west": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_x,
                                                                   boundary_conditions_types_x, north=north,
                                                                   east=east, south=south, west=west)

# Set Boundary Conditions for vel_x
mesh.new_boundary_condition("vel_x", point_index=xy_indices, values=xy_values, type_of_boundary=xy_types)

# -------------------------------------------------- VEL Y -------------------------------------------------------------
# Define boundary conditions for vel_y
boundary_conditions_values_y = {
    "north": 0.,
    "east": 0.,
    "south": 0.,
    "west": vel,
}

boundary_conditions_types_y = {
    "north": True,
    "east": False,
    "south": True,
    "west": True,
}

# Get vectors of nodes values
xy_indices, xy_values, xy_types = TccLib.build_boundary_conditions(mesh, boundary_conditions_values_y,
                                                                   boundary_conditions_types_y, north=north,
                                                                   east=east, south=south, west=west)

# Set Boundary Conditions for vel_y
mesh.new_boundary_condition("vel_y", point_index=xy_indices, values=xy_values, type_of_boundary=xy_types)

# Solve for FEM velocity field solution
velocity_x, velocity_y = TccLib.solve_velocity_field(mesh, dt=dt, total_time=total_time, save_each_frame=True,
                                                     stop_criteria=1e-5)
TccLib.util.save(velocity_x, "vel_x")
TccLib.util.save(velocity_y, "vel_y")

# velocity_x = TccLib.util.load("vel_x")
# velocity_y = TccLib.util.load("vel_y")

# Show results in quiver plot
mesh.show_velocity_quiver(velocity_x, velocity_y)

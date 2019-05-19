import TccLib

import numpy as np

# Define boundary conditions and parameters
vel_top = 1.
dt = 01.
total_time = 50.
particle_density = 3e4
particle_diameter = 1e-3

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Lid_driven_ref", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=False, save=True)

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
#                                                      stop_rule=1e-5)
# TccLib.util.save(velocity_x, "vel_x")
# TccLib.util.save(velocity_y, "vel_y")
velocity_x = TccLib.util.load("vel_x")
velocity_y = TccLib.util.load("vel_y")

# Show results in quiver plot
mesh.show_velocity_quiver(velocity_x, velocity_y)

# Define Particles
particle_a = TccLib.Particle("A", (0.5, 0.95), density=particle_density, diameter=particle_diameter)
particle_b = TccLib.Particle("B", (0.5, 0.90), density=particle_density, diameter=particle_diameter)
particle_c = TccLib.Particle("C", (0.5, 0.85), density=particle_density, diameter=particle_diameter)
particle_d = TccLib.Particle("D", (0.5, 0.80), density=particle_density, diameter=particle_diameter)
particle_e = TccLib.Particle("E", (0.5, 0.75), density=particle_density, diameter=particle_diameter)
particles = [particle_a, particle_b, particle_c, particle_d, particle_e]

# particles = TccLib.util.load("particles")
mesh.add_particle(list_of_particles=particles)

# Define dt based on convergence limit
total_time = 5.
dt = min(particle_a.max_dt(mesh.viscosity), 1e-4)/2**1.

# Define x vector of positions
x_vector = np.arange(0, total_time, dt)

# Move Particles
for time in x_vector:
    print("\rMoving particles {0:.2f}%".format(100 * time / total_time), end="")
    TccLib.move_particles(mesh, velocity_x=velocity_x.last, velocity_y=velocity_y.last, dt=dt)
    TccLib.util.save(particles, "particles")

# Show particles trajectories
mesh.show_particle_course()

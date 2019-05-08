import TccLib

import numpy as np

# Define boundary conditions and parameters
vel = 1.
dt = 0.1
total_time = 20.
particle_density = 3e4
particle_diameter = 1e-3

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "super_oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Channel", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=False, save=True)

# Define analytic comparison expression
analytic_expression = lambda y: 4 * 1.5 * vel / (mesh.length_y**2) * (mesh.length_y - y) * y

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
# velocity_x, velocity_y = TccLib.solve_velocity_field(mesh, dt=dt, total_time=total_time, save_each_frame=True,
#                                                      stop_rule=1e-5)
# TccLib.util.save(velocity_x, "vel_x")
# TccLib.util.save(velocity_y, "vel_y")

velocity_x = TccLib.util.load("vel_x")
velocity_y = TccLib.util.load("vel_y")

# Show results in quiver plot
# mesh.show_velocity_quiver(velocity_x, velocity_y)

# Define x vector of positions
# x_vector = np.linspace(min(mesh.y), max(mesh.y), 100)
# small_dict = velocity_x.reduced_dict_log(5)
#
# # Find values of property in the mesh at a determined set position for every value in the vector of x
# x_position = mesh.length_x*0.5775
# y_vector = {key: [mesh.get_interpolated_value([x_position, x], small_dict[key]) for x in x_vector]
#             for key in small_dict}
#
# # Show comparison graph
# TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
#                             analytic_label="Solução Analítica", title="Equação de Hagen-Poiseuille",
#                             x_label="Posição no Eixo Y(m)", y_label="Valor da Velocidade (m/s)",
#                             save_file_as="{0}_validation".format(mesh.name))

# Define Particles
x_0 = 0.1 * mesh.length_x
y_0 = mesh.length_y / 6

particle_a = TccLib.Particle("A", (x_0, 5 * y_0), density=particle_density, diameter=particle_diameter)
particle_b = TccLib.Particle("B", (x_0, 4 * y_0), density=particle_density, diameter=particle_diameter)
particle_c = TccLib.Particle("C", (x_0, 3 * y_0), density=particle_density, diameter=particle_diameter)
particle_d = TccLib.Particle("D", (x_0, 2 * y_0), density=particle_density, diameter=particle_diameter)
particle_e = TccLib.Particle("E", (x_0, 1 * y_0), density=particle_density, diameter=particle_diameter)
particles = [particle_a, particle_b, particle_c, particle_d, particle_e]

particles = TccLib.util.load("particles")
mesh.add_particle(list_of_particles=particles)

# Define dt based on convergence limit
total_time = 6.
dt = min(particle_a.max_dt(mesh.viscosity), 1e-4)/2**1.

# Define x vector of positions
x_vector = np.arange(0, total_time, dt)

# Move Particles
# for time in x_vector:
#     print("\rMoving particles {0:.2f}%".format(100 * time / total_time), end="")
#     TccLib.move_particles(mesh, velocity_x=velocity_x.last, velocity_y=velocity_y.last, dt=dt)
#     TccLib.util.save(particles, "particles")

# Show particles trajectories
mesh.show_particle_course()

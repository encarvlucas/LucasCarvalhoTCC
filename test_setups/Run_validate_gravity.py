import TccLib

import numpy as np

# Define boundary conditions and parameters
total_time = 0.4
TccLib.Particle.frame_skips = 1
force = "gravitational"
particle_density = 3e4
particle_diameter = 1e-3

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "water"

# Import gmsh created mesh, and set velocity field
mesh = TccLib.Mesh("Forces", liquid=liquid)
vel_x = np.zeros(mesh.size) + 1.0
vel_y = np.zeros(mesh.size)

# Show mesh geometry
# mesh.show_geometry(names=True)

# Define analytic comparison expression
analytic_expression = lambda t: -9.80665/2. * t**2 + mesh.length_y

# Define Particles
particle_a = TccLib.Particle("A", (0.07 * mesh.length_x, mesh.length_y), density=particle_density,
                             diameter=particle_diameter)
mesh.add_particle(list_of_particles=[particle_a])

# Define dt based on convergence limit
dt = min(particle_a.max_dt(mesh.viscosity), 1e-4)/2**4.

# Define x vector of positions
x_vector = np.arange(0, total_time, dt)

# Move Particles
for time in x_vector:
    print("\rMoving particles {0:.2f}%".format(100 * time / total_time), end="")
    TccLib.move_particles(mesh, velocity_x=vel_x, velocity_y=vel_y, dt=dt, single_force=force)

print("\rFinished moving particles")

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = mesh.length_x*0.5
y_vector = np.array(particle_a.position_history)[:-1, 1]

# Show comparison graph
TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Força Gravitacional",
                            x_label="Tempo(s)", y_label="Posição no Eixo Y(m)",
                            save_file_as="{0}_{1}_validation".format(mesh.name, force))

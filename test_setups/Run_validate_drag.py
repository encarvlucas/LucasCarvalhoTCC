import TccLib

import numpy as np

# Define boundary conditions and parameters
total_time = 0.4
TccLib.Particle.frame_skips = 1
force = "drag"
particle_density = 3e4
particle_diameter = 1e-3
vel_const = 2.0
vel_x_0 = 0.
vel_y_0 = 0.

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "water"

# Import gmsh created mesh, and set velocity field
mesh = TccLib.Mesh("Forces", liquid=liquid)
vel_x = np.zeros(mesh.size) + vel_const
vel_y = np.zeros(mesh.size)

# Show mesh geometry
# mesh.show_geometry(names=True)

# Define Particles
x_0 = 0.5 * mesh.length_x
y_0 = 0.8 * mesh.length_y
particle_a = TccLib.Particle("A", (x_0, y_0), density=particle_density, diameter=particle_diameter,
                             velocity=(vel_x_0, vel_y_0))
mesh.add_particle(list_of_particles=[particle_a])

# Define analytic comparison expression
m = particle_a.mass
c = 3 * np.pi * mesh.viscosity * particle_a.diameter
if c/m > 1:
    print("Particle conditions might cause an unexpected behavior!")
analytic_expression = lambda t: m/c*(vel_x_0*(1 - np.exp(-c*t/m)) + vel_const*(np.exp(-c*t/m) - 1)) + vel_const*t + x_0

# Define dt based on convergence limit
dt = min(particle_a.max_dt(mesh.viscosity), 1e-4)/2**6.

# Define x vector of positions
x_vector = np.arange(0, total_time, dt)

# Move Particles
for time in x_vector:
    print("\rMoving particles {0:.2f}%".format(100 * time / total_time), end="")
    TccLib.move_particles(mesh, velocity_x=vel_x, velocity_y=vel_y, dt=dt, single_force=force)

print("\rFinished moving particles")

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = mesh.length_x*0.5
y_vector = np.array(particle_a.position_history)[:-1, 0]

# Show comparison graph
TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Força de Arrasto",
                            x_label="Tempo(s)", y_label="Posição no Eixo x(m)",
                            save_file_as="{0}_{1}_validation".format(mesh.name, force))

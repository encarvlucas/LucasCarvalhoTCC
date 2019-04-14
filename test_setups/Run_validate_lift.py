import TccLib

import numpy as np

# Define boundary conditions and parameters
total_time = 0.4
TccLib.Particle.frame_skips = 1
force = "lift"
particle_density = 3e4
particle_diameter = 1e-3
u_const = 2.0
v_const = -1.0
du_dy = 5.0
vel_y_0 = -0.10

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "water"

# Import gmsh created mesh, and set velocity field
mesh = TccLib.Mesh("Forces", liquid=liquid)
vel_x = np.array([du_dy*y for y in mesh.y]) + u_const
vel_y = np.zeros(mesh.size) + v_const

# Show mesh geometry
# mesh.show_geometry(names=True)

# Define Particles
y_0 = 0.95*mesh.length_y
particle_a = TccLib.Particle("A", (0.0, y_0), density=particle_density,
                             diameter=particle_diameter, velocity=(0., vel_y_0))
mesh.add_particle(list_of_particles=[particle_a])

# Define analytic comparison expression
re_g = np.sqrt(particle_a.diameter**2 * mesh.density / mesh.viscosity * abs(du_dy))
m = particle_a.mass
c = 1.61 * mesh.viscosity * particle_a.diameter * re_g
if c/m > 1:
    print("Particle conditions might cause an unexpected behavior!")
analytic_expression = lambda t: (m/c) * (vel_y_0-v_const) * (np.exp(c*t/m) - 1) + v_const*t + y_0

# Define dt based on convergence limit
dt = min(particle_a.max_dt(mesh.viscosity), 1e-4)/2**5.

# Define x vector of positions
x_vector = np.arange(0, total_time, dt)

# Move Particles
for time in x_vector:
    print("\rMoving particles {0:.2f}%".format(100 * time / total_time), end="")
    TccLib.move_particles(mesh, velocity_x=vel_x, velocity_y=vel_y, dt=dt, single_force=force)

print("\rFinished moving particles")

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = mesh.length_x*0.5
y_vector = {
    "Solução Numérica":        np.array(particle_a.position_history)[:-1, 1],
    "Queda sem Sustentação":   vel_y_0*x_vector + y_0,
}

# Show comparison graph
TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
                            analytic_label="Solução Analítica", title="Força de Sustentação",
                            x_label="Tempo(s)", y_label="Posição no Eixo Y(m)",
                            save_file_as="{0}_{1}_validation".format(mesh.name, force))

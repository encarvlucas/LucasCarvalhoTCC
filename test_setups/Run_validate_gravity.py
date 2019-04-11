import TccLib

import numpy as np

# Define boundary conditions and parameters
vel = 1.0
dt = 01.
total_time = 50.

# Set liquid parameters or declare liquid
# density = 1e3
# viscosity = 0.89e-3
liquid = "oil"

# Import gmsh created mesh
mesh = TccLib.Mesh("Forces", liquid=liquid)

# Show mesh geometry
# mesh.show_geometry(names=True)

# Define analytic comparison expression
analytic_expression = lambda t: -9.80665/2. * t**2

# Define Particles
mesh.add_particle("A", (0.07 * mesh.length_x, 0.8 * mesh.length_y), density=5.4e5, diameter=5e-5)

# Move Particles
TccLib.move_particles(mesh, velocity_x=(np.zeros(mesh.size) + 1.0), velocity_y=np.zeros(mesh.size))

# Define x vector of positions
x_vector = np.linspace(0, total_time, 100)

# Find values of property in the mesh at a determined set position for every value in the vector of x
x_position = mesh.length_x*0.6
y_vector = 1

# Show comparison graph
# TccLib.util.show_comparison(x_vector, analytic_expression, y_vector, numeric_label="Solução Numérica",
#                             analytic_label="Solução Analítica", title="Escoamento de Cavidade",
#                             x_label="Posição no Eixo de Comparação(m)", y_label="Valor da Velocidade (m/s)",
#                             save_file_as="{0}_validation".format(mesh.name))

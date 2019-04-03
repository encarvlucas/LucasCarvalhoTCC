import TccLib

import numpy as np

mesh = TccLib.Mesh("Poisson")

# mesh.show_geometry(names=True)

xy_indices, xy_values, xy_types = TccLib.border_temperature_boundary_conditions(mesh)

mesh.new_boundary_condition("space", point_index=xy_indices, values=xy_values,
                            type_of_boundary=True)
mesh.new_boundary_condition("time", point_index=xy_indices, values=xy_values,
                            type_of_boundary=True)

temperature = TccLib.solve_poisson(mesh, permanent_solution=True)
# mesh.show_3d_solution(temperature)

# temperature = TccLib.solve_poisson(mesh, permanent_solution=False)
# mesh.show_animated_3d_solution(temperature)

x_position = (max(mesh.x) - min(mesh.x))/2.

x_vector = np.arange(min(mesh.y), max(mesh.y), (max(mesh.x) - min(mesh.x))/100.)
y_vector = [mesh.get_interpolated_value([x_position, x], temperature) for x in x_vector]

TccLib.util.show_comparison(x_vector, y_vector, lambda x: x ** 2)

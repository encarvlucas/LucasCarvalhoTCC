from TccLib import *

malha = Mesh("Poiseuille")

malha.show_geometry(names=True)

xy_indices, xy_values, xy_types = border_temperature_boundary_conditions(malha)

malha.new_boundary_condition("space", point_index=xy_indices, values=xy_values,
                             type_of_boundary=True)
malha.new_boundary_condition("time", point_index=xy_indices, values=xy_values,
                             type_of_boundary=True)

vect = solve_poisson(malha, permanent_solution=True)
malha.show_3d_solution(vect)

Q = ComplexPointList([32, 39, 64, 67, 68, 70], 50.)

vect = solve_poisson(malha, permanent_solution=True, q=Q)
malha.show_3d_solution(vect)

vect = solve_poisson(malha, permanent_solution=False)
malha.show_animated_3d_solution(vect)

vect = solve_poisson(malha, permanent_solution=False, q=Q)
malha.show_animated_3d_solution(vect)

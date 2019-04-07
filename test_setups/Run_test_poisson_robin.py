from TccLib import *

malha = Mesh("Poisson")

malha.show_geometry(names=True, save=True)

xy_indices, xy_values, xy_types = build_boundary_conditions(malha)

malha.new_boundary_condition("space", point_index=xy_indices, values=xy_values,
                             type_of_boundary=xy_types)
malha.new_boundary_condition("time", point_index=xy_indices, values=xy_values,
                             type_of_boundary=xy_types)


vect = solve_poisson(malha, permanent_solution=True)
malha.show_3d_solution(vect)

Q = ComplexPointList([32, 39, 64, 67, 68, 70], 500.)

vect = solve_poisson(malha, permanent_solution=True, q=Q)
malha.show_3d_solution(vect)

vect = solve_poisson(malha, permanent_solution=False)
malha.show_animated_3d_solution(vect)

vect = solve_poisson(malha, permanent_solution=False, q=Q)
malha.show_animated_3d_solution(vect)
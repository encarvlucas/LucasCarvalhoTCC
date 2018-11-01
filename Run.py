from TccLib import *
import numpy as np
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

malha = Mesh("Poisson")

# xx, yy = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 15))
# xx = np.array(np.reshape(xx, (xx.size, 1)))
# yy = np.array(np.reshape(yy, (yy.size, 1)))
# xy = np.hstack((xx, yy))

# rand = np.random.rand(200)
# rand_2 = np.random.rand(200)
# xy = list(zip(rand, rand_2)).append([0,0])
# xy = np.array(([0, 0], [1, 0],
#                [1, 1], [0, 1]))
# xy = np.vstack((xy, np.array(list(zip(rand, rand_2)))))
# malha.show_geometry(names=True)

xy_indices, xy_values, xy_types = border_temperature_boundary_conditions(malha)

malha.new_boundary_condition("space", point_index=xy_indices, values=xy_values,
                             type_of_boundary=True)
malha.new_boundary_condition("time", point_index=xy_indices, values=xy_values,
                             type_of_boundary=True)
Q = ComplexPointList([32, 39, 64, 67, 68, 70], 50.)
permanent = False
# quit()
vect = solve_poisson(malha, permanent_solution=permanent)
if permanent:
    malha.show_3d_solution(vect)
else:
    # malha.show_solution(vect[-1])
    malha.show_animated_3d_solution(vect)

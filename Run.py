from TccLib import *
import numpy as np
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

malha = Mesh()

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
malha.import_point_structure(import_mesh_file="results/untitled")
# malha.show_geometry(names=True)

from collections import OrderedDict as od
vertex_a = np.where(malha.x == np.min(malha.x))[0]
vertex_b = np.where(malha.y == np.min(malha.y))[0]
vertex_c = np.where(malha.x == np.max(malha.x))[0]
vertex_d = np.where(malha.y == np.max(malha.y))[0]
xy_indices = list(od.fromkeys(np.append(vertex_a, vertex_b)))
xy_values = np.zeros(len(xy_indices))
xy_type = np.zeros(len(xy_values)) + 1
xy_indices = list(od.fromkeys(np.append(xy_indices, list(od.fromkeys(np.append(vertex_c, vertex_d))))))
xy_type = np.hstack((xy_type, np.zeros(len(xy_indices) - len(xy_values))))
xy_values = np.append(xy_values, np.zeros(len(xy_indices) - len(xy_values)) + 1)
malha.space_boundary_conditions.set_new_boundary_conditions(point_index=xy_indices, values=xy_values,
                                                            type_of_boundary=xy_type)
malha.time_boundary_conditions.set_new_boundary_conditions(point_index=xy_indices, values=xy_values,
                                                           type_of_boundary=xy_type)
Q = ComplexPointList([32, 39, 64, 67, 68, 70], 5.)
permanent = True
vect = solve_poisson(malha, permanent_solution=permanent)
if permanent:
    malha.show_solution(vect)
else:
    # malha.show_solution(vect[-1])
    malha.show_animated_solution(vect)

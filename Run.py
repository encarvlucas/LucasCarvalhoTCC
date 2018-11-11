from TccLib import *
import numpy as np
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D


xx, yy = np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 1, 20))
xx = np.array(np.reshape(xx, (xx.size, 1)))
yy = np.array(np.reshape(yy, (yy.size, 1)))
xy = np.hstack((xx, yy))

malha = Mesh("Poiseuille")
# , points=list(xy))

# rand = np.random.rand(200)
# rand_2 = np.random.rand(200)
# xy = list(zip(rand, rand_2)).append([0,0])
# xy = np.array(([0, 0], [1, 0],
#                [1, 1], [0, 1]))
# xy = np.vstack((xy, np.array(list(zip(rand, rand_2)))))
# malha.show_geometry(names=True, save=True)

list(map(lambda _vect: malha.new_boundary_condition(_vect["name"], point_index=_vect["indices"], values=_vect["values"],
                                                    type_of_boundary=_vect["type"]),
         hagen_poiseuille_boundary_conditions(malha)))

poiseuille = True

if poiseuille:
    vel_x, vel_y = solve_poiseuille(malha, total_time=5., dt=0.9)
    # malha.show_velocity_solution(vel_x, vel_y)

else:
    Q = ComplexPointList([32, 39, 64, 67, 68, 70], 50.)
    permanent = False
    vect = solve_poisson(malha, permanent_solution=permanent)
    if permanent:
        malha.show_3d_solution(vect)
    else:
        # malha.show_solution(vect[-1])
        malha.show_animated_3d_solution(vect)

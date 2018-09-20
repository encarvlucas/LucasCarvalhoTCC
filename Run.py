from TccLib import *
import numpy as np

malha = Mesh()

xx, yy = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 15))
xx = np.array(np.reshape(xx, (xx.size, 1)))
yy = np.array(np.reshape(yy, (yy.size, 1)))
xy = np.hstack((xx, yy))
malha.import_point_structure(points=list(xy))
# malha.show(rainbow=True)

from collections import OrderedDict as od
vertex_a = np.where(xx == np.min(xx))[0]
vertex_b = np.where(yy == np.min(yy))[0]
vertex_c = np.where(xx == np.max(xx))[0]
vertex_d = np.where(yy == np.max(yy))[0]
xy_indices = list(od.fromkeys(np.append(vertex_a, vertex_b)))
xy_values = np.zeros(len(xy_indices))
xy_indices = list(od.fromkeys(np.append(xy_indices, list(od.fromkeys(np.append(vertex_c, vertex_d))))))
xy_values = np.append(xy_values, np.zeros(len(xy_indices) - len(xy_values)) + 1)
malha.space_boundary_conditions.set_new_boundary_conditions(point_index=xy_indices, values=xy_values)
malha.time_boundary_conditions.set_new_boundary_conditions(point_index=range(xy.size), values=0)
vect = solve(malha)

from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
# size = (malha.size, malha.size)
# pyplot.pcolormesh(malha.x, malha.y, vect, cmap='jet', vmin=min(vect), vmax=max(vect))
# pyplot.colorbar()
# pyplot.show()

fig = pyplot.gcf()
axes = Axes3D(fig)
surf = axes.plot_trisurf(malha.x, malha.y, vect, cmap="jet")
axes.view_init(90, 270)
fig.colorbar(surf, shrink=0.4, aspect=9)
pyplot.show()

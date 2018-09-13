from TccLib import *

malha = Mesh()
# malha.show(rainbow=True)

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

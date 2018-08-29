# -*- coding: utf-8 -*-
# -------------------------------------- by: LUCAS CARVALHO DE SOUSA ---------------------------------------------------
# Esta biblioteca foi criada para o a solução de sistemas diferenciais através do Método de Elementos Finitos
# This library was created for educational purposes, it is meant to be used for the solution of diferential equation
# systems
__author__ = "Lucas Carvalho de Sousa"


# -- Functions ---------------------------------------------------------------------------------------------------------

def import_point_structure(*args):
    """
    Imports points position to create mesh from source file
    :param args: Name of source file, defaults to points.txt
    :return:
    """
    if not args:
        filename = "points.txt"
    else:
        filename = args[0]

    try:
        with open(filename, "r") as arq:
            points = []
            for line in arq:
                points.append([float(i) for i in line.split(";")] + [0.0])
            surface = create_new_surface(points)

    except FileNotFoundError:
        surface = create_new_surface()

    return surface


def create_new_surface(*imported_points, lt_version=False):
    """
    Create new surface
    :return: Element information
    """
    import numpy as np
    import scipy.spatial as dl
    if not lt_version:
        import pygmsh
        geom = pygmsh.built_in.Geometry()

    x, y, ien = 0, 0, 0

    if imported_points:
        # Custom geometry.
        imported_points = np.array(imported_points[0])
        delauney_surfaces = dl.Delaunay(imported_points[:, :2])

        if lt_version:
            x = delauney_surfaces.points[:, 0]
            y = delauney_surfaces.points[:, 1]
            ien = delauney_surfaces.simplices
        else:
            for tri in delauney_surfaces.simplices:
                geom.add_polygon([[delauney_surfaces.points[tri[0]][0], delauney_surfaces.points[tri[0]][1], 0.0],
                                  [delauney_surfaces.points[tri[1]][0], delauney_surfaces.points[tri[1]][1], 0.0],
                                  [delauney_surfaces.points[tri[2]][0], delauney_surfaces.points[tri[2]][1], 0.0]])

    else:
        if not lt_version:

            # Default surface.
            geom.add_polygon([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ])

    if not lt_version:
        import meshio

        # Saving mesh as .vtk exportable file
        points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
        meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data)

        x = points[:, 0]
        y = points[:, 1]
        ien = cells["triangle"]

    return x, y, ien


# -- Classes -----------------------------------------------------------------------------------------------------------

class Mesh:
    """
    Mesh element to be used in the calculations
    """
    x, y, ien = 0, 0, 0  # TODO: USE SPARCE MATRIX

    """
    Class constructor,
    initializes geometry
    """
    def __init__(self):
        self.x, self.y, self.ien = import_point_structure()

    """
    Display mesh geometry on screen using matplotlib
    """
    def show(self, rainbow=False):
        import numpy as np
        import matplotlib.pyplot as plt

        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        if rainbow:
            plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.hsv(np.linspace(0.0, 1.0, len(self.ien)))))

            for triangl in self.ien:
                plot_coordinates = (self.x[triangl[0]], self.x[triangl[1]], self.x[triangl[2]], self.x[triangl[0]]), \
                                   (self.y[triangl[0]], self.y[triangl[1]], self.y[triangl[2]], self.y[triangl[0]])

                plt.plot(plot_coordinates[0], plot_coordinates[1])

        else:
            plt.triplot(self.x, self.y, triangles=self.ien[0])

        plt.show()

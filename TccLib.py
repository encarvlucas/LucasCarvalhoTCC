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


def create_new_surface(*points):
    """
    Create new surface
    :return: Element information
    """
    import pygmsh
    import meshio
    import numpy as np
    import scipy.spatial as dl

    geom = pygmsh.built_in.Geometry()

    if points:
        # Custom geometry.
        points = np.array(points[0])
        delauney_surfaces = dl.Delaunay(points[:, :2])
        for tri in delauney_surfaces.simplices:
            geom.add_polygon([[delauney_surfaces.points[tri[0]][0], delauney_surfaces.points[tri[0]][1], 0.0],
                              [delauney_surfaces.points[tri[1]][0], delauney_surfaces.points[tri[1]][1], 0.0],
                              [delauney_surfaces.points[tri[2]][0], delauney_surfaces.points[tri[2]][1], 0.0]])
    else:
        # Default surface.
        geom.add_polygon([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])

    # Saving mesh as .vtk exportable file
    points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
    meshio.write_points_cells('output.vtk', points, cells, point_data, cell_data, field_data)
    return points, cells, point_data, cell_data, field_data


# -- Classes -----------------------------------------------------------------------------------------------------------

class Mesh:
    """
    Mesh element to be used in the calculations
    """
    import numpy as np
    x = np.zeros(1)  # TODO: USE SPARCE MATRIX

    def __init__(self):
        self.mesh = import_point_structure()

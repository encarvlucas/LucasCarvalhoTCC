from TccLib import *

malha = Mesh()
# geom = pygmsh.built_in.Geometry()
#
# # Draw a cross.
# poly = geom.add_polygon([
#     [0.0,  0.0, 0.0],
#     [1.0,  0.0, 0.0],
#     [0.5,  0.5, 0.0],
#     [1.0,  1.0, 0.0],
#     [0.0,  1.0, 0.0],
#     ])
#
# test = geom.get_code()
#
# with open("teste.geo", "w") as arq:
#     arq.write(test)
#
# points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom)
#
# meshio.write_points_cells('test.vtk', points, cells, point_data, cell_data, field_data)

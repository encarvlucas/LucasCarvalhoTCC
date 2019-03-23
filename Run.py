import numpy as np
import TccLib

xx, yy = np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 1, 20))
xx = np.array(np.reshape(xx, (xx.size, 1)))
yy = np.array(np.reshape(yy, (yy.size, 1)))
xy = np.hstack((xx, yy))

malha = TccLib.Mesh("Poiseuille_refined")
# malha = Mesh(points=list(xy))

# rand = np.random.rand(200)
# rand_2 = np.random.rand(200)
# xy = list(zip(rand, rand_2)).append([0,0])
# xy = np.array(([0, 0], [1, 0],
#                [1, 1], [0, 1]))
# xy = np.vstack((xy, np.array(list(zip(rand, rand_2)))))
# malha.show_geometry(names=True, save=True)

list(map(lambda _vect: malha.new_boundary_condition(_vect["name"], point_index=_vect["indices"], values=_vect["values"],
                                                    type_of_boundary=_vect["type"]),
         TccLib.hagen_poiseuille_boundary_conditions(malha)))

# --------------------------------- Adding particles ---------------------------------------------------------------
malha.add_particle("A", (0.07 * (max(malha.x) - min(malha.x)) + min(malha.x), 0.5 * (max(malha.y) - min(malha.y)) +
                         min(malha.y)), density=430, diameter=1e-7, velocity=(1., 0.))
# particles = [TccLib.Particle("B", (0.11 * (max(malha.x) - min(malha.x)) + min(malha.x),
#                                    0.8 * (max(malha.y) - min(malha.y)) + min(malha.y)),
#                              color="b", density=15e2, diameter=1e-7, velocity=(0.8, 0.)),
#              TccLib.Particle("C", (0.21 * (max(malha.x) - min(malha.x)) + min(malha.x),
#                                    0.7 * (max(malha.y) - min(malha.y)) + min(malha.y)),
#                              color="g", density=250, diameter=5e-6)]
# malha.add_particle(list_of_particles=particles)

poiseuille = True

if poiseuille:
    vel_x, vel_y = TccLib.solve_poiseuille(malha, total_time=1.15, dt=.01, rho_coef=1e3, mu_coef=0.89e-3,
                                           save_each_frame=True)
    # malha.show_velocity_quiver(vel_x, vel_y)
    malha.show_particle_movement(save=True)

else:
    Q = TccLib.ComplexPointList([32, 39, 64, 67, 68, 70], 50.)
    permanent = False
    vect = TccLib.solve_poisson(malha, permanent_solution=permanent)
    if permanent:
        malha.show_3d_solution(vect)
    else:
        # malha.show_solution(vect[-1])
        malha.show_animated_3d_solution(vect)

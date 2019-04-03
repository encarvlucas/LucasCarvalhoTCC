import TccLib

import numpy as np

xx, yy = np.meshgrid(np.linspace(0, 5, 10), np.linspace(0, 1, 20))
xx = np.array(np.reshape(xx, (xx.size, 1)))
yy = np.array(np.reshape(yy, (yy.size, 1)))
xy = np.hstack((xx, yy))

malha = TccLib.Mesh("Poiseuille", density=1e3, viscosity=0.89e-3)
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
malha.add_particle("A", (0.07 * (max(malha.x) - min(malha.x)) + min(malha.x), 0.8 * (max(malha.y) - min(malha.y)) +
                         min(malha.y)), density=5.4e5, diameter=5e-5, velocity=(0.5, -0.5))
# particles = [TccLib.Particle("B", (0.11 * (max(malha.x) - min(malha.x)) + min(malha.x),
#                                    0.8 * (max(malha.y) - min(malha.y)) + min(malha.y)),
#                              color="b", density=7.5e6, diameter=1e-4, velocity=(0.8, 0.)),
#              TccLib.Particle("C", (0.21 * (max(malha.x) - min(malha.x)) + min(malha.x),
#                                    0.7 * (max(malha.y) - min(malha.y)) + min(malha.y)),
#                              color="g", density=2.5e6, diameter=3e-5)]
# malha.add_particle(list_of_particles=particles)

poiseuille = True

if poiseuille:
    # velocity_time = 5.
    #
    # vel_x, vel_y = TccLib.solve_velocity_field(malha, total_time=velocity_time, dt=.01, save_each_frame=False)
    # # malha.show_velocity_quiver(vel_x, vel_y)
    # malha.output_results(result_dictionary={"Velocity_X": vel_x, "Velocity_Y": vel_y})
    #
    # TccLib.util.save((vel_x, vel_y), "vel")
    vel_x, vel_y = TccLib.util.load("vel")

    particle_time = 5.
    particle_dt = 1e-5
    TccLib.Particle.frame_skips = int(particle_time / (60.*particle_dt))
    for t in np.arange(0, particle_time, particle_dt):
        print("\rMoving Particles {0:.2f}%".format(100 * t / particle_time), end="")

        # Move particles
        TccLib.move_particles(malha, velocity=(vel_x, vel_y), dt=particle_dt)

        # Show particle progress
        # mesh.show_geometry()

    print("\rMoving Particles done!")
    TccLib.util.save(malha, "last_mesh")
    print("Generating gif with {0} frames.".format(len(malha.particles[0].position_history)))
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

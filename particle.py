import util

import numpy as np


class Particle:
    """
    Defines a moving particle.
    Particles are defined as having a spherical shape of constant diameter.
    """
    position_history = None

    reynolds = 0.

    # Number of skipped frames between saved positions.
    frame_skips = 1000

    def __init__(self, name, position=(0., 0.), velocity=(0., 0.), density=1., diameter=0.1, color="r"):
        """
        Particle class constructor.
        :param name: Name of the particle.
        :param position: Positional argument, (x, y) coordinates [m].
        :param velocity: Initial particle velocity, (x, y) values [m/s].
        :param density: Particle density [kg/m³].
        :param diameter: Diameter of particle (approximated as a sphere) [m].
        :param color: Custom color for particle, default is "r" (red).
        """

        util.check_method_call(name)

        self.name = name
        self.pos_x = position[0]
        self.pos_y = position[1]
        self.position_history = [(self.pos_x, self.pos_y)]

        self.velocity_x = velocity[0]
        self.velocity_y = velocity[1]

        self.density = density
        self.diameter = diameter
        # Density * Volume of a sphere = rho * pi/6 * d^3
        self.mass = np.pi / 6. * self.density * self.diameter ** 3

        self.color = color

        # Number of times this particle's position has been calculated.
        self.time_count = 1

    @property
    def velocity(self):
        return np.array([self.velocity_x, self.velocity_y])

    def apply_forces(self, forces, mesh, dt):
        """
        Method that gather all the forces applied to a particle and calculates the movement of the particle.
        :param forces: Forces as a dictionary with their names as keys
                       and the tuple of force values in (x, y) [N || kg.m/s²].
        :param mesh: Mesh object that determines the boundaries for collision.
        :param dt: The time difference between frames [s].
        """
        util.check_if_instance(forces, dict)

        sum_forces_x = 0.
        sum_forces_y = 0.

        for force in forces.values():
            util.check_if_instance(force, (list, tuple, np.ndarray))
            sum_forces_x += force[0]
            sum_forces_y += force[1]

        self.velocity_x += sum_forces_x * dt / self.mass
        self.velocity_y += sum_forces_y * dt / self.mass
        self.mark_new_position(mesh, dt)

    def mark_new_position(self, mesh, dt):
        """
        Assigns a new position for the particle and tracks its location history.
        :param mesh: Mesh object that determines the boundaries for collision.
        :param dt: The time difference between frames [s].
        """
        # TODO: ADD GENERIC WALL COLLISIONS
        self.pos_x += self.velocity_x * dt
        if not mesh.contains_particle(self):
            self.pos_x = 0.

        self.pos_y += self.velocity_y * dt
        if not mesh.contains_particle(self):
            self.pos_y = 0.7*abs(self.pos_y)
            self.velocity_y = 0.7*abs(self.velocity_y)

        self.time_count += 1
        if self.time_count % self.frame_skips == 0:
            self.position_history.append((self.pos_x, self.pos_y))

    @property
    def pixel_size(self):
        """
        Returns the plot size os the particle.
        :return: Approximate plotting size.
        """
        return self.diameter * 1e7

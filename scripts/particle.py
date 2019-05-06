import scripts.util as util

import numpy as np


class Particle:
    """
    Defines a moving particle.
    Particles are defined as having a spherical shape of constant diameter.
    """
    # Number of skipped frames between saved positions.
    frame_skips = 1000

    # Defines size of particles when plotting
    resolution_multiplier = 1e4

    color_wheel = ["b", "g", "r", "c", "m", "y"]
    color_index = 0

    def __init__(self, name: str, position: (list, tuple) = (0., 0.), velocity: (list, tuple) = (0., 0.),
                 density: float = 1., diameter: float = 0.1, color: str = None):
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
        self.last_velocity = self.velocity

        self.density = density
        self.diameter = diameter
        # Density * Volume of a sphere = rho * pi/6 * d^3
        self.volume = np.pi / 6. * self.diameter ** 3
        self.mass = self.volume * self.density

        # Particle Reynolds
        self.reynolds = 0.

        self.color = color if color is not None else Particle.get_next_color()

        # Number of times this particle's position has been calculated.
        self.time_count = 1

    @property
    def velocity(self):
        """
        :return: Particle's current velocity as a ndarray.
        """
        return np.array([self.velocity_x, self.velocity_y])

    @property
    def radius(self):
        """
        Particle radius, defined as 1/2 of it's diameter.
        :return: Particle's radius.
        """
        return self.diameter / 2.

    @property
    def position(self):
        """
        Particle current position coordinates.
        :return: Particle's position.
        """
        return self.pos_x, self.pos_y

    def max_dt(self, fluid_viscosity):
        """
        Determines the max dt to still allow convergence.
        :param fluid_viscosity: Fluid's viscosity.
        :return: Maximum allowed dt.
        """
        return self.diameter**2 * self.density / (9. * fluid_viscosity)

    def apply_forces(self, forces: dict, mesh, dt: float, single_force: str = None):
        """
        Method that gather all the forces applied to a particle and calculates the movement of the particle.
        :param forces: Forces as a dictionary with their names as keys
                       and the tuple of force values in (x, y) [N || kg.m/s²].
        :param mesh: Mesh object that determines the boundaries for collision.
        :param dt: The time difference between frames [s].
        :param single_force: Parameter used for testing single forces one at a time.
        """
        util.check_if_instance(forces, dict)

        sum_forces_x = 0.
        sum_forces_y = 0.

        if single_force is not None:
            if single_force in forces:
                sum_forces_x += forces[single_force][0]
                sum_forces_y += forces[single_force][1]
            else:
                raise KeyError(single_force + " force is not available.")
        else:
            for force in forces.values():
                sum_forces_x += force[0]
                sum_forces_y += force[1]

        self.last_velocity = self.velocity
        self.velocity_x += sum_forces_x * dt / self.mass
        self.velocity_y += sum_forces_y * dt / self.mass
        self.mark_new_position(mesh, dt)

    def mark_new_position(self, mesh, dt):
        """
        Assigns a new position for the particle and tracks its location history.
        :param mesh: Mesh object that determines the boundaries for collision.
        :param dt: The time difference between frames [s].
        """
        last_position = self.position
        self.pos_x += self.velocity_x * dt
        self.pos_y += self.velocity_y * dt
        if not mesh.contains_particle(self):
            self.pos_x = last_position[0] - 0.7*self.velocity_y * dt
            self.pos_y = last_position[1] + 0.7*self.velocity_x * dt
            if not mesh.contains_particle(self):
                self.pos_x = last_position[0] + 0.7*self.velocity_y * dt
                self.pos_y = last_position[1] - 0.7*self.velocity_x * dt

        self.time_count += 1
        if self.time_count % self.frame_skips == 0:
            self.position_history.append(self.position)

    @property
    def pixel_size(self):
        """
        Returns the plot size os the particle.
        :return: Approximate plotting size.
        """
        return self.diameter * Particle.resolution_multiplier

    @classmethod
    def get_next_color(cls):
        current = cls.color_wheel[cls.color_index]
        cls.color_index = 0 if cls.color_index + 1 >= len(cls.color_wheel) else cls.color_index + 1
        return current

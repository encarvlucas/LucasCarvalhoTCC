from boundaryConditions import *
from particle import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

import os
import fnmatch


class Mesh:
    """
    Mesh element to be used in the calculations
    """
    x, y, ien, delauney_surfaces = None, None, None, None
    particles = None
    name = "default"
    size = 0
    default_dt = 0.
    density = None
    viscosity = None

    def __init__(self, name: str = "untitled", points=None, density: float = 1.0, viscosity: float = 1.0):
        """
        Class constructor, initializes geometry.
        :param name: Mesh's main name.
        :param density: Fluid density [kg/m³].
        :param viscosity: Fluid dynamic viscosity [Pa.s || kg/m.s].
        """
        from shutil import copy
        self.name = name
        self.particles = []
        self.density = density
        self.viscosity = viscosity

        if isinstance(points, list):
            self.import_point_structure(points=points)
        else:
            self.import_point_structure(import_mesh_file=self.name)

        self.boundary_conditions = {}

        try:
            os.chdir("./results/{0}/".format(self.name))
        except FileNotFoundError:
            os.mkdir("./results/{0}".format(self.name))
            os.chdir("./results/{0}/".format(self.name))

        if points is not None:
            copy("../{0}.msh".format(self.name), "./")

    def import_point_structure(self, *args, points=None, light_version=True, import_mesh_file=""):
        """
        Imports points position to create mesh from source file.
        :param args: Name of source file, defaults to "points.txt".
        :param points: Custom set of defined points, in form of list of shape (x, 2).
        :param light_version: Defines if script will use Gmsh to obtain elements.
        :param import_mesh_file: Use public library meshio to import already generated ".msh" file.
        """
        if not args:
            filename = "points.txt"
        else:
            filename = args[0]

        if isinstance(points, list):
            surface = util.create_new_surface(points, lt_version=light_version)

        else:
            try:
                if import_mesh_file:
                    surface = util.use_meshio("results/{0}".format(import_mesh_file), None)

                else:
                    with open(filename, "r") as arq:
                        points = []
                        for line in arq:
                            points.append([float(i) for i in line.split(";")] + [0.0])
                        surface = util.create_new_surface(points, lt_version=light_version)

            except FileNotFoundError:
                print("File not found, generating new default mesh.")
                # surface = create_new_surface(lt_version=False)  # TODO: FIX LIBRARY
                surface = util.use_meshio("results/untitled", None)

        self.x, self.y, self.ien, self.delauney_surfaces = surface
        self.size = len(self.x)
        self.default_dt = util.get_dt(self)

    def new_boundary_condition(self, name, point_index=range(0), values=0., type_of_boundary=True):
        """
        Creates a new entry in the boundary conditions dictionary.
        :param name: Name of the property that has these boundary conditions.
        :param point_index: Vector of order of points.
        :param values: Value or vector of values of the condition in each point.
        :param type_of_boundary: Value or vector of values, defined: True for Dirichlet and False for Neumann.
        """
        self.boundary_conditions[name] = BoundaryConditions()
        self.boundary_conditions[name].set_new_boundary_condition(point_index=point_index, values=values,
                                                                  type_of_boundary=type_of_boundary)

    def add_particle(self, name: str = None, position: (list, tuple) = None, velocity=(0., 0.), density=1.,
                     diameter=0.1, color="r", list_of_particles=None):
        """
        Associates a new particle with the mesh.
        :param name: Name of the particle. Each particle defined in a mesh must have different names.
        :param position: Positional argument, (x, y) coordinates [m].
        :param velocity: Initial particle velocity, (x, y) values [m/s].
        :param density: Particle density [kg/m³].
        :param diameter: Diameter of particle (approximated as a sphere) [m].
        :param color: Custom color for particle, default is "r" (red).
        :param list_of_particles: Adds a list of predefined particle objects.
        """
        if list_of_particles:
            list_of_particles = util.check_list_or_object(list_of_particles, Particle)
            [self.particles.append(particle) for particle in list_of_particles if particle.name not in
             [obj.name for obj in self.particles]]

        else:
            util.check_method_call(name)
            if name in [obj.name for obj in self.particles]:
                print("\nThere is a particle already with that name. "
                      "Each particle defined inside a mesh must have a different name.\n")

            else:
                util.check_method_call(position)

                self.particles.append(Particle(name=name, position=position, velocity=velocity, density=density,
                                               diameter=diameter, color=color))

    def contains_particle(self, particle: Particle):
        """
        Checks if a particle is inside the mesh geometric domain.
        :rtype: bool
        :param particle: Particle object that contains it's current position.
        :return: True if it is inside, otherwise False.
        """
        util.check_if_instance(particle, Particle)

        if self.delauney_surfaces.find_simplex((particle.pos_x, particle.pos_y)) < 0:
            return False
        else:
            return True

    def get_fluid_velocity(self, position, velocity_vector_x, velocity_vector_y):
        """
        Method that calculates the velocity of the fluid in the desired coordinates using interpolation of the closest
        points know velocities.
        :param position: Coordinates of the point (x, y) [m].
        :param velocity_vector_x: Vector of velocity in the x axis [m/s].
        :param velocity_vector_y: Vector of velocity in the y axis [m/s].
        :return: The interpolated velocity of the fluid at the coordinates.
        """
        util.check_if_instance(position, (list, tuple))

        # Determine which element contains the particle
        element = self.ien[self.delauney_surfaces.find_simplex(position)]

        total_area = util.get_area(self.x[element], self.y[element])
        # Interpolate velocity value for particle coordinates

        # Adds the velocity component from each point, interpolated by the proportion of the area opposite it.
        fluid_vel_x = 0.
        fluid_vel_y = 0.
        for point in element:
            _element = element.tolist()
            _element.remove(point)
            component = abs(util.get_area([position[0]] + list(self.x[_element]),
                                          [position[1]] + list(self.y[_element])) / total_area)
            fluid_vel_x += component * velocity_vector_x[point]
            fluid_vel_y += component * velocity_vector_y[point]

        return fluid_vel_x, fluid_vel_y

    def remove_previous_results(self, return_to_previous_directory=False):
        """
        Removes previous results to prevent incorrect data analysis.
        :param return_to_previous_directory: Returns to previous directory.
        """
        # ---------------- Deleting previous results -------------------------------------------------------------------
        try:
            os.chdir("./paraview_results/")
        except FileNotFoundError:
            if not fnmatch.fnmatch(os.getcwd(), "*/paraview_results"):
                os.mkdir("paraview_results/")
                os.chdir("./paraview_results/")

        list(map(os.remove,
                 [file for file in os.listdir(".") if fnmatch.fnmatch(file, "{0}_results_*.vtk".format(self.name))]))

        if return_to_previous_directory:
            os.chdir("..")

    def output_results(self, result_vector=None, extension="VTK", dt=0., frame_num=None, data_names="Temperature",
                       result_dictionary=None):
        """
        Export result to .vtk or .csv file
        :param result_vector: The vector of the value in each point
        :param extension: File extension
        :param dt: Value of time between frames
        :param frame_num: Save one frame at a time (the number of that frame).
        :param data_names: Data name for display.
        :param result_dictionary: A dictionary that contains the names for each data and its values.
        """
        # ------------------ Contingency -------------------------------------------------------------------------------
        util.check_method_call(result_vector, result_dictionary)

        if not result_dictionary:
            result_dictionary = {}
            if isinstance(data_names, list):
                for index, data_name in enumerate(data_names):
                    result_dictionary[data_name] = result_vector[index]
            else:
                result_dictionary[data_names] = result_vector

        number_frames = 0
        for _name, _vector in result_dictionary.items():
            if " " in _name:
                raise AttributeError("There can't be any spaces in the property names!")

            if dt and frame_num is None:
                if len(_vector[0]) != self.size:
                    raise ValueError("Incorrect size for result _vector.")
                number_frames = len(_vector)

            elif len(_vector) != self.size:
                raise ValueError("Incorrect size for result vector.")

        # ----------------------- Change to result directory -----------------------------------------------------------
        try:
            os.chdir("./paraview_results/")
        except FileNotFoundError:
            if not fnmatch.fnmatch(os.getcwd(), "*/paraview_results"):
                os.mkdir("paraview_results/")
                os.chdir("./paraview_results/")

        number_elements = len(self.ien)

        if extension == "CSV":
            # ------------------------- Saving results to CSV file -----------------------------------------------------
            with open("{0}_results.csv".format(self.name), "w") as arq:
                arq.write("Points:0, Points:1, Points:2, {0}\n".format(data_names))
                for i in range(self.size):
                    arq.write("{0},{1},{2},{3}\n".format(self.x[i], self.y[i], 0, result_vector[i]))
            return os.chdir("..")

        if extension == "VTK":
            size = self.size

            def write_scalars(file, property_name, vector):
                file.write("\nSCALARS {0} float 1\n".format(property_name))
                file.write("\nLOOKUP_TABLE {0}\n".format(property_name))
                for _i in range(size):
                    file.write("{}\n".format(vector[_i]))

            def write_header_and_cells(file):
                # ------------------------------------ Header ----------------------------------------------------------
                file.write("# vtk DataFile Version 3.0\n{0}\n{1}\n\nDATASET {2}\n".format("LucasCarvalhoTCC Results",
                                                                                          "ASCII", "POLYDATA"))
                if dt:
                    file.write("FIELD FieldData 1\nTIME 1 1 double\n{}\n".format(dt))
                # ------------------------------------ Points coordinates ----------------------------------------------
                file.write("\nPOINTS {0} {1}\n".format(size, "float"))
                for _i in range(size):
                    file.write("{0} {1} 0.0\n".format(self.x[_i], self.y[_i]))
                # --------------------------------------- Cells --------------------------------------------------------
                file.write("\nPOLYGONS {0} {1}\n".format(number_elements, number_elements * 4))
                for _i in range(number_elements):
                    file.write("{0} {1} {2} {3}\n".format(3, self.ien[_i][0], self.ien[_i][1], self.ien[_i][2]))

            if dt:
                if frame_num is not None:
                    # ----- Saving a single result (frame) to VTK file -------------------------------------------------
                    with open("{0}_results_{1}.vtk".format(self.name, frame_num), "w") as arq:
                        # --------------------------------- Header and cells -------------------------------------------
                        write_header_and_cells(arq)

                        # ------------------------------------ Data in each point --------------------------------------
                        arq.write("\nPOINT_DATA {0}\n".format(size))
                        [write_scalars(arq, name, vector) for name, vector in result_dictionary.items()]

                else:
                    self.remove_previous_results()
                    # ----- Saving multiple results to VTK files -------------------------------------------------------
                    for j in range(number_frames):
                        with open("{0}_results_{1}.vtk".format(self.name, j), "w") as arq:
                            # ----------------------------- Header and cells -------------------------------------------
                            write_header_and_cells(arq)

                            # -------------------------------- Data in each point --------------------------------------
                            arq.write("\nPOINT_DATA {0}\n".format(size))
                            [write_scalars(arq, name, vector[j]) for name, vector in result_dictionary.items()]

            else:
                # ------------------------- Saving results to VTK file -------------------------------------------------
                with open("{0}_results.vtk".format(self.name), "w") as arq:
                    # --------------------------------- Header and cells -----------------------------------------------
                    write_header_and_cells(arq)

                    # ----------------------------------- Data in each point--------------------------------------------
                    arq.write("\nPOINT_DATA {0}\n".format(size))
                    [write_scalars(arq, name, vector) for name, vector in result_dictionary.items()]

            return os.chdir("..")

        raise NameError("Format not available. Try VTK or CSV.")

    def show_geometry(self, show=True, particles: list = None, names=False, rainbow=False, save=False):
        """
        Display mesh geometry on screen using matplotlib.
        :param show: Display plot on screen when its generated.
        :param particles: Particle or list of particles to be displayed inside the geometry.
        :param names: Show the index of each point next to it.
        :param rainbow: Color in the edge of each element in a different color.
        :param save: Save generated image.
        :return: Display image.
        """
        particles = util.check_list_or_object(particles, Particle)
        particles.extend(self.particles)

        # Draw mesh points
        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        plt.gca().autoscale(False)

        # Draw particles
        [plt.scatter(particle.pos_x, particle.pos_y, s=particle.pixel_size, c=particle.color) for
         particle in particles if self.contains_particle(particle)]

        # Sets plot styling
        util.style_plot(self.x, self.y)

        # Display id's of points and names of objects
        if names:
            for _index in range(self.size):
                plt.gca().annotate(_index, (self.x[_index], self.y[_index]))

        [plt.gca().annotate(particle.name,
                            (particle.pos_x + particle.diameter / 10, particle.pos_y + particle.diameter / 10),
                            color=particle.color, fontsize=15)
         for particle in particles if self.contains_particle(particle)]

        # Displays elements as different colors to help distinguish each one
        if rainbow:
            plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.hsv(np.linspace(0.0, 1.0, len(self.ien)))))

            for element in self.ien:
                plot_coordinates = ((self.x[element[0]], self.x[element[1]], self.x[element[2]], self.x[element[0]]),
                                    (self.y[element[0]], self.y[element[1]], self.y[element[2]], self.y[element[0]]))

                plt.plot(plot_coordinates[0], plot_coordinates[1])

        else:
            plt.triplot(self.x, self.y, triangles=self.ien[0])

        if save:
            plt.savefig("{0}_mesh".format(self.name))

        if show:
            return plt.show()

    def save_frame(self, frame_num: int):
        """
        Saves frame image in the frames folder, to allow a better problem interpretation.
        :param frame_num: Current frame number
        """
        try:
            os.mkdir("./frames")
        except FileExistsError:
            pass

        self.show_geometry(show=False)

        plt.savefig("./frames/{0}_frame_{1}".format(self.name, frame_num))
        plt.close()

    def show_3d_solution(self, solution_vector):
        """
        Display 3D solution of the mesh geometry.
        :param solution_vector: Vector that contains the value of the solution for each point in the mesh.
        :return: Display image.
        """
        util.check_method_call(solution_vector)
        try:
            if len(solution_vector) != self.size:
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Solution must be a vector")

        self.output_results(result_dictionary={"Temperature": solution_vector})

        fig = plt.gcf()
        axes = Axes3D(fig)
        surf = axes.plot_trisurf(self.x, self.y, solution_vector, cmap="jet")
        axes.view_init(90, 270)
        fig.colorbar(surf, shrink=0.4, aspect=9)

        plt.savefig("{0}_permanent_results".format(self.name))

        return plt.show()

    def show_animated_3d_solution(self, frames_vector, dt=0.):
        """
        Display animated version of the 3D solution.
        :param frames_vector: Vector which each element contains a vector with the value of the solution for each point
                              in the mesh.
        :param dt: Time between each frame, if not specified the animation won't be saved.
        :return: Display image.
        """
        util.check_method_call(frames_vector)
        try:
            if np.any([len(i) != self.size for i in frames_vector]):
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                             "{0}".format(self.size))

        fig = plt.gcf()
        axes = Axes3D(fig)
        _min_value = np.min(frames_vector)
        _max_value = np.max(frames_vector)
        surf = axes.plot_trisurf(self.x, self.y, frames_vector[0], cmap="jet", vmin=_min_value, vmax=_max_value)
        fig.colorbar(surf, shrink=0.4, aspect=9)

        def update(_current_frame):
            plt.cla()
            axes.plot_trisurf(self.x, self.y, _current_frame, cmap="jet", vmin=_min_value, vmax=_max_value)
            axes.set_zlim3d([_min_value, _max_value])
            return

        animation = FuncAnimation(fig, update, frames=frames_vector, interval=100, save_count=False)

        if dt:
            self.output_results(result_dictionary={"Temperature": frames_vector}, dt=dt)
            animation.save("{0}_transient_results.gif".format(self.name), dpi=80, writer='imagemagick')

        return plt.show()

    def show_velocity_quiver(self, velocity_x, velocity_y):
        """
        Display velocity vector field as arrows that represent the intensity and the direction of the velocity of the
        fluid at that point.
        :param velocity_x: Vector of velocity in the x axis.
        :param velocity_y: Vector of velocity in the y axis.
        :return: Display image.
        """
        util.check_method_call(velocity_x, velocity_y)

        try:
            if len(velocity_x) != self.size or len(velocity_y) != self.size:
                raise ValueError("Incorrect size of solution vector, it must be the same size as the mesh: "
                                 "{0}".format(self.size))
        except TypeError:
            raise ValueError("Solution must be a vector")

        self.output_results(result_dictionary={"Velocity_X": velocity_x, "Velocity_Y": velocity_y})

        fig, axes = plt.subplots()
        axes.quiver(self.x, self.y, velocity_x, velocity_y)

        # Sets plot styling
        util.style_plot(self.x, self.y)

        plt.savefig("{0}_velocity_field".format(self.name))

        return plt.show()

    def show_particle_movement(self, save=False):
        """
        Displays an animated image of the particles trajectories.
        :param save: Determines if the generated image will be saved.
        :return: Display image.
        """
        # Draw mesh points
        plt.plot(self.x, self.y, marker=".", color="k", linestyle="none", ms=5)

        number_of_frames = max([len(self.particles[i].position_history) for i in range(len(self.particles))])
        list_of_dots = {}

        for particle in self.particles:
            list_of_dots[particle] = plt.scatter(0, 0, s=particle.pixel_size, c=particle.color)

        def update(_frame):
            # Draw particles
            [_dot.set_offsets((_particle.position_history[_frame][0], _particle.position_history[_frame][1]))
             for _particle, _dot in list_of_dots.items() if _frame < len(_particle.position_history)]

            return

        animation = FuncAnimation(plt.gcf(), update, frames=number_of_frames, interval=1, save_count=False)

        if save:
            animation.save("{0}_particle_movement.gif".format(self.name), dpi=80, writer='imagemagick')

        return plt.show()

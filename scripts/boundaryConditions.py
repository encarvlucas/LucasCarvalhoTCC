import scripts.util as util
import numpy as np


class BoundaryConditions:
    """
    Boundary conditions of the simulation.
    """

    def __init__(self):
        """
        Class constructor, initializes vector values.
        """
        self.point_index_vector = []
        self.values_vector = []
        self.type_of_condition_vector = []

    def set_new_boundary_condition(self, *vect_arg, point_index: (list, np.ndarray) = None,
                                   values: (int, float, list, np.ndarray) = 0.,
                                   type_of_boundary: (bool, list, np.ndarray) = True):
        """
        Sets the boundary condition for the mesh.
        :param vect_arg: Vector(s) of boundary conditions.
        :param point_index: Vector of order of points.
        :param values: Value or vector of values of the condition in each point.
        :param type_of_boundary: Value or vector of values, defined: True for Dirichlet and False for Neumann.
        """
        if vect_arg:
            util.check_method_call(vect_arg)
        else:
            util.check_method_call(point_index)

        try:
            if isinstance(vect_arg[0], list) or isinstance(values, np.ndarray):
                array = np.array(vect_arg[0])

                if len(vect_arg[0][0]) == 3:
                    self.point_index_vector = np.append(self.point_index_vector, array[:, 0])
                    self.values_vector = np.append(self.values_vector, array[:, 1])
                    self.type_of_condition_vector = np.append(self.type_of_condition_vector,
                                                              list(map(bool, array[:, 2])))
                    return

                if len(vect_arg[0][0]) == 2:
                    self.point_index_vector = np.append(self.point_index_vector, array[:, 0], )
                    self.values_vector = np.append(self.values_vector, array[:, 1], )
                    self.type_of_condition_vector = np.append(self.type_of_condition_vector,
                                                              [True] * len(self.point_index_vector))
                    return

        except (TypeError, IndexError):
            self.point_index_vector = np.append(self.point_index_vector, np.array(point_index)).astype("int")
            if isinstance(values, (list, np.ndarray)):
                if len(values) == len(point_index):
                    self.values_vector = np.append(self.values_vector, np.array(values))
                else:
                    raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                     "types and definitions")
            else:
                self.values_vector = np.append(self.values_vector, np.array([values] * len(point_index)))

            if isinstance(type_of_boundary, (list, np.ndarray)):
                if len(type_of_boundary) == len(point_index):
                    self.type_of_condition_vector = np.append(self.type_of_condition_vector,
                                                              list(map(bool, np.array(type_of_boundary))))

                else:
                    raise ValueError("Incorrect vector sizes, there must be an equal number of points and point "
                                     "types and definitions")
            else:
                self.type_of_condition_vector = np.append(self.type_of_condition_vector,
                                                          np.array([type_of_boundary] * len(point_index)))

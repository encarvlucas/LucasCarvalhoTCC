import numpy as np


class PropertyState:
    """
    Class that represents a state of a mesh property in a given time frame.
    """
    def __init__(self, values: np.ndarray, time: float):
        """
        Class constructor that creates a frame containing the information of the property in the mesh at a single time.
        :param values: Array of property values in the declared timestamp.
        :param time: Timestamp of frame [s].
        """
        self.values = values
        self.time = time
        self.min = min(self.values)
        self.max = max(self.values)

    def __len__(self):
        """
        Method that returns size of values array to allow list behavior.
        :return: Size of values array.
        """
        return len(self.values)


class MeshPropertyStates:
    """
    Class that represents the previous states of a mesh property in various time frames.
    """
    def __init__(self, states_list: [list, np.ndarray] = None, states_dict: dict = None):
        """
        Class constructor that
        :param states_list: Optional argument, list of arrays or single array of values of the property in order.
                            If a list of arrays is used, the default dt used is 1.0s between frames.
        :param states_dict: Optional argument, constructs the object using the dict keys as the timestamps and dict
                            values as the property values arrays.
        """
        self.dt = 1.0
        if states_dict is not None:
            self.states = [PropertyState(states_dict[state], state) for state in states_dict]
            self.dt = self.states[-1].time - self.states[-2].time

        elif states_list is not None:
            if isinstance(states_list, list):
                self.states = [PropertyState(state, float(index)) for index, state in enumerate(states_list)]
            else:
                self.states = [PropertyState(states_list, 0.)]

        else:
            self.states = []

    def append(self, state: np.ndarray, time: float):
        self.states.append(PropertyState(state, time))
        self.dt = self.states[-1].time - self.states[-2].time

    @property
    def min(self):
        """
        :return: Smallest value of the property of all timestamps.
        """
        return min([state.min for state in self.states])

    @property
    def max(self):
        """
        :return: Largest value of the property of all timestamps.
        """
        return max([state.max for state in self.states])

    @property
    def dict(self):
        """
        :return: Dictionary of object information, using each timestamp as key and the value array as dict value.
        """
        return dict([(state.time, state.values) for state in self.states])

    def reduced_dict(self, size: int):
        """
        Method that creates a dictionary containing the condensed information of this object.
        :param size: Number of elements in dictionary.
        :return: Dictionary of object information, using each timestamp as key and the value array as dict value.
        """
        small_dict = {self.states[0].time: self.states[0].values}
        for i in [int(len(self.states)*j/(size-1)) - 1 for j in range(size)]:
            small_dict[self.states[i].time] = self.states[i].values
        return small_dict

    def reduced_dict_log(self, size: int):
        """
        Method that creates a dictionary containing the condensed information of this object, using a log scale for
        timestamps, except the last one is always the highest timestamp available.
        :param size: Number of elements in dictionary.
        :return: Dictionary of object information, using each timestamp as key and the value array as dict value.
        """
        small_dict = {
            self.states[0].time: self.states[0].values,
            self.states[-1].time: self.states[-1].values,
        }
        for i in [int(len(self.states)*np.log(j)/(size-1)**2) - 1 for j in range(2, size)]:
            small_dict[self.states[i].time] = self.states[i].values
        return small_dict

    def value_at(self, time: float):
        # TODO: FINISH METHOD
        return

    def __getitem__(self, item):
        """
        Method used to help reproduce a list behavior.
        :param item: Item index.
        :return: Values array of the selected item.
        """
        return self.states[item].values

    def __len__(self):
        """
        Method that returns size of values array to allow list behavior.
        :return: Size of list of states. Number of frames present.
        """
        return len(self.states)

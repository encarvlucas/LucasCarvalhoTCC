import util


class ComplexPointList:
    """
    Class that defines a list of points with an index value and a property value each.
    """
    indexes = None
    values = None

    def __init__(self, _indexes, _values):
        """
        Class constructor that generates one list for each property
        :param _indexes: The point index in the ordered list of the mesh point coordinates
        :param _values: The value of the desired property in each point.
        """
        util.check_method_call(_indexes)

        self.indexes = _indexes

        try:
            if isinstance(list(_values), list):
                self.values = _values
            else:
                raise TypeError

        except TypeError:
            self.values = []
            for i in range(len(_indexes)):
                self.values.append(_values)

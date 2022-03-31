from asali.utils.input_parser import ReactorType, ResolutionMethod

import numpy as np


class SolutionParser:
    def __init__(self):
        """
        Class to parse solution of reactor model
        """
        self._reactor_type = None
        self._resolution_method = None
        self._gas = None
        self._surf = None
        self._y = None
        self._x = None
        self._length = None
        self._is_solved = False

    def get_reactor_type(self):
        """
        Reactor model getter
        :return: Reactor model
        """
        return self._reactor_type

    def set_reactor_type(self, value):
        """
        Reactor model setter
        :param value: Reactor model
        :return:
        """
        self._reactor_type = value

    # Creating a property object for Reactor model
    reactor_type = property(get_reactor_type, set_reactor_type)

    def get_resolution_method(self):
        """
        Resolution method getter
        :return: Resolution method
        """
        return self._resolution_method

    def set_resolution_method(self, value):
        """
        Resolution method setter
        :param value: Resolution method
        :return:
        """
        self._resolution_method = value

    # Creating a property object for Resolution method
    resolution_method = property(get_resolution_method, set_resolution_method)

    def get_gas(self):
        """
        Gas object getter
        :return: Gas object
        """
        return self._gas

    def set_gas(self, value):
        """
        Gas object setter
        :param value: Gas object
        :return:
        """
        self._gas = value

    # Creating a property object for Gas object
    gas = property(get_gas, set_gas)

    def get_surf(self):
        """
        Surface object getter
        :return: Surface object
        """
        return self._surf

    def set_surf(self, value):
        """
        Surface object setter
        :param value: Surface object
        :return:
        """
        self._surf = value

    # Creating a property object for Surface object
    surf = property(get_surf, set_surf)

    def get_x(self):
        """
        Independent variable getter
        :return: Independent variable
        """
        return self._x

    def set_x(self, value):
        """
        Independent variable setter
        :param value: Independent variable
        :return:
        """
        self._x = value

    # Creating a property object for Independent variable
    x = property(get_x, set_x)

    def get_y(self):
        """
        Dependent variables getter
        :return: Dependent variables
        """
        return self._y

    def set_y(self, value):
        """
        Dependent variables setter
        :param value: Dependent variables
        :return:
        """
        self._y = value

    # Creating a property object for Dependent variables
    y = property(get_y, set_y)

    def get_is_solved(self):
        """
        Solved bool getter
        :return: Solved bool
        """
        return self._is_solved

    def set_is_solved(self, value):
        """
        Solved bool variables setter
        :param value: Solved bool
        :return:
        """
        self._is_solved = value

    # Creating a property object for Dependent variables
    is_solved = property(get_is_solved, set_is_solved)

    def get_length(self):
        """
        Length getter
        :return: Length
        """
        return self._length

    def set_length(self, value):
        """
        Length setter
        :param value: Length
        :return:
        """
        self._length = value

    # Creating a property object for Length
    length = property(get_length, set_length)

    def get_mass_fraction(self):
        """
        Get mass fraction
        :return: Vector/Matrix representing the resulting mass fraction
        """
        if self._reactor_type == ReactorType.BATCH or self._reactor_type == ReactorType.CSTR:
            return self._y[:, :self._gas.n_species]

        if self._reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, :self.gas.n_species]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                mass_fraction = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.surf.n_species + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    mass_fraction[i] = sol_for_time[:, :self.gas.n_species]

                return mass_fraction

        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, :self.gas.n_species]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1
                mass_fraction = np.zeros([self._x.size], dtype=np.ndarray)
                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    mass_fraction[i] = sol_for_time[:, :self.gas.n_species]

                return mass_fraction

    def get_mass_fraction_wall(self):
        """
        Get mass fraction of wall phase
        :return: Vector/Matrix representing the resulting mass fraction
        """
        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, self.gas.n_species:self.gas.n_species + self.gas.n_species]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1
                mass_fraction = np.zeros([self._x.size], dtype=np.ndarray)
                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    mass_fraction[i] = sol_for_time[:, self.gas.n_species:self.gas.n_species + self.gas.n_species]

                return mass_fraction

    def get_mole_fraction(self):
        """
        Get mole fraction
        :return: Vector/Matrix representing the resulting mole fraction
        """
        if self._reactor_type == ReactorType.BATCH or self._reactor_type == ReactorType.CSTR:
            mole_fraction = np.zeros([len(self._x), self._gas.n_species], dtype=np.float64)

            for i in range(0, len(self._x)):
                self._gas.Y = self._y[i, :self._gas.n_species]
                mole_fraction[i, :self._gas.n_species] = self._gas.X

            return mole_fraction

        if self._reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                mole_fraction = np.zeros([len(self._x), self._gas.n_species], dtype=np.float64)

                for i in range(0, len(self._x)):
                    self._gas.Y = self._y[i, :self._gas.n_species]
                    mole_fraction[i, :self._gas.n_species] = self._gas.X

                return mole_fraction

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                mole_fraction = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.surf.n_species + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    mass_fraction_vector = sol_for_time[:, :self.gas.n_species]
                    mole_fraction[i] = np.zeros_like(mass_fraction_vector)
                    for j, mass_fraction in enumerate(mass_fraction_vector):
                        self.gas.Y = mass_fraction
                        mole_fraction[i][j, :] = self.gas.X

                return mole_fraction

        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                mole_fraction = np.zeros([len(self._x), self._gas.n_species], dtype=np.float64)

                for i in range(0, len(self._x)):
                    self._gas.Y = self._y[i, :self._gas.n_species]
                    mole_fraction[i, :self._gas.n_species] = self._gas.X

                return mole_fraction

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                mole_fraction = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1
                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    mass_fraction_vector = sol_for_time[:, :self.gas.n_species]
                    mole_fraction[i] = np.zeros_like(mass_fraction_vector)
                    for j, mass_fraction in enumerate(mass_fraction_vector):
                        self.gas.Y = mass_fraction
                        mole_fraction[i][j, :] = self.gas.X

                return mole_fraction

    def get_mole_fraction_wall(self):
        """
        Get mass fraction of wall phase
        :return: Vector/Matrix representing the resulting mass fraction
        """
        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                mole_fraction = np.zeros([len(self._x), self._gas.n_species], dtype=np.float64)

                for i in range(0, len(self._x)):
                    self._gas.Y = self._y[i, self.gas.n_species:self.gas.n_species + self.gas.n_species]
                    mole_fraction[i, :self._gas.n_species] = self._gas.X

                return mole_fraction

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                mole_fraction = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1
                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    mass_fraction_vector = sol_for_time[:, self.gas.n_species:self.gas.n_species + self.gas.n_species]
                    mole_fraction[i] = np.zeros_like(mass_fraction_vector)
                    for j, mass_fraction in enumerate(mass_fraction_vector):
                        self.gas.Y = mass_fraction
                        mole_fraction[i][j, :] = self.gas.X

                return mole_fraction

    def get_temperature(self):
        """
        Get temperature
        :return: Vector/Matrix representing the resulting temperature at a fixed axial/time position
        """
        if self._reactor_type == ReactorType.BATCH or self._reactor_type == ReactorType.CSTR:
            return self._y[:, -1]

        if self._reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, -1]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                temperature = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.surf.n_species + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    temperature[i] = sol_for_time[:, -1]

                return temperature

        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, -2]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                temperature = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    temperature[i] = sol_for_time[:, -2]

                return temperature

    def get_temperature_wall(self):
        """
        Get temperature of wall phase
        :return: Vector/Matrix representing the resulting temperature at a fixed axial/time position
        """
        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, -1]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                temperature = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    temperature[i] = sol_for_time[:, -1]

                return temperature

    def get_coverage(self):
        """
        Get coverage
        :return: Vector/Matrix representing the resulting coverage
        """
        if self._reactor_type == ReactorType.BATCH:
            return self._y[:, self._gas.n_species + 1:self._gas.n_species + 1 + self._surf.n_species]

        if self._reactor_type == ReactorType.CSTR:
            return self._y[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]

        if self._reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                coverage = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.surf.n_species + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    coverage[i] = sol_for_time[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]

                return coverage

        if self._reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self._resolution_method == ResolutionMethod.STEADYSTATE:
                return self._y[:,
                       self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species]

            if self._resolution_method == ResolutionMethod.TRANSIENT:
                coverage = np.zeros([self._x.size], dtype=np.ndarray)
                nv = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

                for i in range(0, self._y.shape[0]):
                    sol_for_time = self._y[i, :].reshape(-1, nv)
                    coverage[i] = sol_for_time[:,
                                  self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species]

                return coverage

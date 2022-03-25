import numpy as np
import os

from pint import UnitRegistry


class UnitConverter:
    def __init__(self):
        """
        Class to convert Unit dimensions
        """
        self.unit_register = UnitRegistry(system='mks')
        self.unit_register.load_definitions(['hour = 60 * minute = h = hr = H',
                                             'kmol = 1000 * mole',
                                             'mbar = 0.001 * bar',
                                             'dm = decimeter',
                                             'cm = centimeter',
                                             'mm = millimeter',
                                             'm3 = m ** 3',
                                             'dm3 = dm ** 3',
                                             'cm3 = cm ** 3',
                                             'mm3 = mm ** 3',
                                             ])

    @staticmethod
    def value_type_handler(value, converter):
        """
        Function that convert a value to another value by using a multiplier
        :param value: Value to be converted
        :param converter: Multiplier
        :return: Converted values
        """
        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def converter(self, value, start_ud, final_ud):
        """
        Convert any value from a starting unit dimension to a final unit dimension
        :param value: Value to be converted
        :param start_ud: Value acutal unit dimension
        :param final_ud: Value final unit dimension
        :return: Value in the new unit dimension
        """
        parsed_start_ud = 1. * self.unit_register.parse_expression(start_ud)
        converter = parsed_start_ud.to(final_ud).magnitude
        return UnitConverter.value_type_handler(value, converter)

    def convert_to_seconds(self, value, ud):
        """
        Convert to second
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in second
        """
        return self.converter(value, ud, 'seconds')

    def convert_to_meter(self, value, ud):
        """
        Convert to meter
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in meter
        """
        return self.converter(value, ud, 'meter')

    def convert_to_cubic_meter(self, value, ud):
        """
        Convert to cubic meter
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in cubic meter
        """
        return self.converter(value, ud, 'meter**3')

    def convert_to_pascal(self, value, ud):
        """
        Convert to Pascal
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in Pascal
        """
        return self.converter(value, ud, 'pascal')

    def convert_to_one_over_meter(self, value, ud):
        """
        Convert to 1/meter
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in 1/meter
        """
        return self.converter(value, ud, '1/meter')

    def convert_to_kg_per_seconds(self, value, ud):
        """
        Convert to kg/s
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in kg/s
        """
        return self.converter(value, ud, 'kg/s')

    def convert_to_cubic_meter_per_seconds(self, value, ud):
        """
        Convert to m3/s
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in m3/s
        """
        return self.converter(value, ud, 'meter**3/s')

    def convert_to_kelvin(self, value, ud):
        """
        Convert to K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in K
        """
        parsed_start_ud = self.unit_register.parse_expression(ud)

        if isinstance(value, float):
            return self.unit_register.Quantity(value, parsed_start_ud).to('kelvin').magnitude

        if isinstance(value, int):
            return self.unit_register.Quantity(value, parsed_start_ud).to('kelvin').magnitude

        return np.asarray(
            [self.unit_register(self.unit_register.Quantity(v, parsed_start_ud)).to('kelvin').magnitude for v in value])

    def convert_to_kg_per_cubic_meter(self, value, ud):
        """
        Convert to kg/m3
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in kg/m3
        """
        return self.converter(value, ud, 'kg/meter**3')

    def convert_to_joule_per_kg_per_kelvin(self, value, ud):
        """
        Convert to J/kg/K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in J/kg/K
        """
        return self.converter(value, ud, 'J/kg/degK')

    def convert_to_watt_per_meter_per_kelvin(self, value, ud):
        """
        Convert to W/m/K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in W/m/K
        """
        return self.converter(value, ud, 'W/meter/degK')

    def convert_from_seconds(self, value, ud):
        """
        Convert from seconds
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value, 'seconds', ud)

    def convert_from_cubic_meter(self, value, ud):
        """
        Convert from m3
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value, 'meter**3', ud)

    def convert_from_pascal(self, value, ud):
        """
        Convert from Pascal
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value, 'pascal', ud)

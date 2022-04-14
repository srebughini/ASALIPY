import numpy as np


class UnitConverter:
    def __init__(self):
        """
        Class to convert Unit dimensions
        """
        self.to_kelvin = {
            'degK': 0.,
            'K': 0.,
            'degC': 273.15}

        self.to_pascal = {
            'pascal': 1.,
            'Pa': 1.,
            'bar': 1e05,
            'MPa': 1e06,
            'kPa': 1e03,
            'torr': 133.322,
            'atm': 101325,
            'mmHg': 133.322,
            'mbar': 1e-03 * 1e05}

        self.from_pascal = {k: 1 / v for k, v in self.to_pascal.items()}

        self.to_seconds = {
            's': 1.,
            'seconds': 1.,
            'second': 1.,
            'min': 60.,
            'minutes': 60.,
            'minute': 60.,
            'h': 3600.,
            'hours': 3600.,
            'hour': 3600.,
            'd': 24. * 3600.,
            'days': 24 * 3600.,
            'day': 24 * 3600.
        }

        self.from_seconds = {k: 1 / v for k, v in self.to_seconds.items()}

        self.to_meter = {
            'meter': 1.,
            'm': 1.,
            'decimeter': 0.1,
            'dm': 0.1,
            'centimeter': 0.01,
            'cm': 0.01,
            'millimeter': 0.001,
            'mm': 0.001,
            'micron': 1e-06
        }

        self.from_meter = {k: 1 / v for k, v in self.to_meter.items()}

        self.to_cubic_meter = {
            'meter**3': 1.,
            'm3': 1.,
            'm**3': 1.,
            'decimeter**3': 1e-03,
            'dm3': 1e-03,
            'dm**3': 1e-03,
            'centimeter**3': 1e-06,
            'cm3': 1e-06,
            'cm**3': 1e-06,
            'millimeter**3': 1e-09,
            'mm3': 1e-09,
            'mm**3': 1e-09
        }

        self.from_cubic_meter = {k: 1 / v for k, v in self.to_cubic_meter.items()}

        self.to_one_over_meter = {'1/' + k: 1 / v for k, v in self.to_meter.items()}

        self.from_one_over_meter = {k: 1 / v for k, v in self.to_one_over_meter.items()}

        self.to_kilograms = {
            'kg': 1.,
            'g': 1e-03,
            'mg': 1e-06
        }

        self.from_kilograms = {k: 1 / v for k, v in self.to_kilograms.items()}

        self.to_joule = {
            'J': 1.,
            'kJ': 1e03,
            'MJ': 1e06,
            'cal': 4.184,
            'kcal': 4.184 * 1e03,
            'kWh': 3.6 * 1e06,
            'Wh': 3.6 * 1e03
        }

        self.from_joule = {k: 1 / v for k, v in self.to_joule.items()}

        self.to_watt = {
            'W': 1.,
            'kW': 1e03,
            'MW': 1e06,
            'cal/s': 4.184,
            'kcal/s': 4.184 * 1e03,
            'cal/min': 4.184 / 60.,
            'kcal/min': 4.184 * 1e03 / 60.,
            'cal/h': 4.184 / 3600.,
            'kcal/h': 4.184 * 1e03 / 3600.
        }

        self.from_watt = {k: 1 / v for k, v in self.to_watt.items()}

        self.to_kilograms_per_seconds = {}
        for kup, vup in self.to_kilograms.items():
            for kdown, vdown in self.to_seconds.items():
                self.to_kilograms_per_seconds[kup + "/" + kdown] = vup / vdown

        self.from_kilograms_per_seconds = {k: 1 / v for k, v in self.to_kilograms_per_seconds.items()}

        self.to_cubic_meter_per_seconds = {}
        for kup, vup in self.to_cubic_meter.items():
            for kdown, vdown in self.to_seconds.items():
                self.to_cubic_meter_per_seconds[kup + "/" + kdown] = vup / vdown

        self.from_cubic_meter_per_seconds = {k: 1 / v for k, v in self.to_cubic_meter_per_seconds.items()}

        self.to_kilograms_per_cubic_meter = {}
        for kup, vup in self.to_kilograms.items():
            for kdown, vdown in self.to_cubic_meter.items():
                self.to_kilograms_per_cubic_meter[kup + "/" + kdown] = vup / vdown

        self.from_kilograms_per_cubic_meter = {k: 1 / v for k, v in self.to_kilograms_per_cubic_meter.items()}

        self.to_joule_per_kilograms = {}
        for kup, vup in self.to_joule.items():
            for kdown, vdown in self.to_kilograms.items():
                self.to_joule_per_kilograms[kup + "/" + kdown] = vup / vdown

        self.from_joule_per_kilograms = {k: 1 / v for k, v in self.to_joule_per_kilograms.items()}

        self.to_watt_per_meter = {}
        for kup, vup in self.to_watt.items():
            for kdown, vdown in self.to_meter.items():
                self.to_watt_per_meter[kup + "/" + kdown] = vup / vdown

        self.from_watt_per_meter = {k: 1 / v for k, v in self.to_watt_per_meter.items()}

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

    @staticmethod
    def converter(value, start_ud, final_ud, to_parser, from_parser):
        """
        Convert any value from a starting unit dimension to a final unit dimension
        :param value: Value to be converted
        :param start_ud: Value actual unit dimension
        :param final_ud: Value final unit dimension
        :param to_parser: Dict that convert to the basic unit dimension (Pa, s, J, W, m, kg, m3, 1/m)
        :param from_parser: Dict that convert from the basic unit dimension (Pa, s, J, W, m, kg, m3, 1/m)
        :return: Value in the new unit dimension
        """
        converter = 1. * to_parser[start_ud] * from_parser[final_ud]
        return UnitConverter.value_type_handler(value, converter)

    def convert_to_seconds(self, value, ud):
        """
        Convert to second
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in second
        """
        return self.converter(value, ud, 'seconds', self.to_seconds, self.from_seconds)

    def convert_to_meter(self, value, ud):
        """
        Convert to meter
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in meter
        """
        return self.converter(value, ud, 'meter', self.to_meter, self.from_meter)

    def convert_to_cubic_meter(self, value, ud):
        """
        Convert to cubic meter
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in cubic meter
        """
        return self.converter(value, ud, 'meter**3', self.to_cubic_meter, self.from_cubic_meter)

    def convert_to_pascal(self, value, ud):
        """
        Convert to Pascal
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in Pascal
        """
        return self.converter(value, ud, 'pascal', self.to_pascal, self.from_pascal)

    def convert_to_one_over_meter(self, value, ud):
        """
        Convert to 1/meter
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in 1/meter
        """
        return self.converter(value, ud, '1/meter', self.to_one_over_meter, self.from_one_over_meter)

    def convert_to_kg_per_seconds(self, value, ud):
        """
        Convert to kg/s
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in kg/s
        """
        return self.converter(value,
                              ud,
                              'kg/s',
                              self.to_kilograms_per_seconds,
                              self.from_kilograms_per_seconds)

    def convert_to_cubic_meter_per_seconds(self, value, ud):
        """
        Convert to m3/s
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in m3/s
        """
        return self.converter(value,
                              ud,
                              'meter**3/s',
                              self.to_cubic_meter_per_seconds,
                              self.from_cubic_meter_per_seconds)

    def convert_to_kelvin(self, value, ud):
        """
        Convert to K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in K
        """
        if isinstance(value, float):
            return value + self.to_kelvin[ud]

        if isinstance(value, int):
            return value + self.to_kelvin[ud]

        return np.asarray([v + self.to_kelvin[ud] for v in value])

    def convert_to_kg_per_cubic_meter(self, value, ud):
        """
        Convert to kg/m3
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in kg/m3
        """
        return self.converter(value,
                              ud,
                              'kg/meter**3',
                              self.to_kilograms_per_cubic_meter,
                              self.from_kilograms_per_cubic_meter)

    def convert_to_joule_per_kg_per_kelvin(self, value, ud):
        """
        Convert to J/kg/K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in J/kg/K
        """

        ud_split = ud.split("/")
        ud_clean = ud_split[0] + "/" + ud_split[1]

        return self.converter(value,
                              ud_clean,
                              'J/kg',
                              self.to_joule_per_kilograms,
                              self.from_joule_per_kilograms)

    def convert_to_watt_per_meter_per_kelvin(self, value, ud):
        """
        Convert to W/m/K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in W/m/K
        """

        ud_split = ud.split("/")
        ud_clean = ud_split[0] + "/" + ud_split[1]

        return self.converter(value,
                              ud_clean,
                              'W/meter',
                              self.to_watt_per_meter,
                              self.from_watt_per_meter)

    def convert_from_seconds(self, value, ud):
        """
        Convert from seconds
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value,
                              'seconds',
                              ud,
                              self.to_seconds,
                              self.from_seconds)

    def convert_from_cubic_meter(self, value, ud):
        """
        Convert from m3
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value,
                              'meter**3',
                              ud,
                              self.to_cubic_meter,
                              self.from_cubic_meter)

    def convert_from_pascal(self, value, ud):
        """
        Convert from Pascal
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value,
                              'pascal',
                              ud,
                              self.to_pascal,
                              self.from_pascal)

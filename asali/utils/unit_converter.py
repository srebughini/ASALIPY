import numpy as np


class UnitConverter:
    def __init__(self):
        """
        Class to convert Unit dimensions
        """
        self.to_kelvin = UnitConverter.temperature_ud()
        self.to_pascal = UnitConverter.pressure_ud()
        self.to_seconds = UnitConverter.time_ud()
        self.to_meter = UnitConverter.length_ud()
        self.to_square_meter = UnitConverter.area_ud()
        self.to_cubic_meter = UnitConverter.volume_ud()
        self.to_kilograms = UnitConverter.mass_ud()
        self.to_kilomole = UnitConverter.mole_ud()
        self.to_joule = UnitConverter.energy_ud()
        self.to_watt = UnitConverter.power_ud()
        self.to_pascal_seconds = UnitConverter.viscosity_ud()

        self.from_pascal = UnitConverter.convert_from_to(self.to_pascal)
        self.from_seconds = UnitConverter.convert_from_to(self.to_seconds)
        self.from_meter = UnitConverter.convert_from_to(self.to_meter)
        self.from_square_meter = UnitConverter.convert_from_to(self.to_square_meter)
        self.from_cubic_meter = UnitConverter.convert_from_to(self.to_cubic_meter)
        self.from_kilograms = UnitConverter.convert_from_to(self.to_kilograms)
        self.from_kilomole = UnitConverter.convert_from_to(self.to_kilomole)
        self.from_joule = UnitConverter.convert_from_to(self.to_joule)
        self.from_watt = UnitConverter.convert_from_to(self.to_watt)
        self.from_pascal_seconds = UnitConverter.convert_from_to(self.to_pascal_seconds)

        self.to_one_over_meter = {'1/' + k: 1 / v for k, v in self.to_meter.items()}
        self.from_one_over_meter = UnitConverter.convert_from_to(self.to_one_over_meter)

        self.to_kilograms_per_seconds, self.from_kilograms_per_seconds = UnitConverter.create_fractional_ud(
            self.to_kilograms,
            self.to_seconds)

        self.to_cubic_meter_per_seconds, self.from_cubic_meter_per_seconds = UnitConverter.create_fractional_ud(
            self.to_cubic_meter,
            self.to_seconds)

        self.to_kilograms_per_cubic_meter, self.from_kilograms_per_cubic_meter = UnitConverter.create_fractional_ud(
            self.to_kilograms,
            self.to_cubic_meter)

        self.to_joule_per_kilograms, self.from_joule_per_kilograms = UnitConverter.create_fractional_ud(
            self.to_joule,
            self.to_kilograms)

        self.to_joule_per_kilomole, self.from_joule_per_kilomole = UnitConverter.create_fractional_ud(
            self.to_joule,
            self.to_kilomole)

        self.to_watt_per_meter, self.from_watt_per_meter = UnitConverter.create_fractional_ud(
            self.to_watt,
            self.to_meter)

        self.to_kilograms_per_kilomole, self.from_kilograms_per_kilomole = UnitConverter.create_fractional_ud(
            self.to_kilograms,
            self.to_kilomole)

        self.to_square_meter_per_seconds, self.from_square_meter_per_seconds = UnitConverter.create_fractional_ud(
            self.to_square_meter,
            self.to_seconds)

        self.to_meter_per_seconds, self.from_meter_per_seconds = UnitConverter.create_fractional_ud(
            self.to_meter,
            self.to_seconds)

    @staticmethod
    def temperature_ud():
        """
        Return temperature unit dimensions as dict
        :return: Temperature unit dimensions as dict
        """
        return {'K': 0.,
                'degC': 273.15}

    @staticmethod
    def pressure_ud():
        """
        Return pressure unit dimensions as dict
        :return: Pressure unit dimensions as dict
        """
        return {'pascal': 1.,
                'Pa': 1.,
                'bar': 1e05,
                'MPa': 1e06,
                'kPa': 1e03,
                'torr': 133.322,
                'atm': 101325,
                'mmHg': 133.322,
                'mbar': 1e-03 * 1e05}

    @staticmethod
    def time_ud():
        """
        Return time unit dimensions as dict
        :return: Time unit dimensions as dict
        """
        return {'s': 1.,
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
                'day': 24 * 3600.}

    @staticmethod
    def length_ud():
        """
        Return length unit dimensions as dict
        :return: Length unit dimensions as dict
        """
        return {'meter': 1.,
                'm': 1.,
                'decimeter': 0.1,
                'dm': 0.1,
                'centimeter': 0.01,
                'cm': 0.01,
                'millimeter': 0.001,
                'mm': 0.001,
                'micron': 1e-06}

    @staticmethod
    def area_ud():
        """
        Return area unit dimensions as dict
        :return: area unit dimensions as dict
        """
        return {'meter**2': 1.,
                'm2': 1.,
                'm**2': 1.,
                'decimeter**2': 1e-02,
                'dm2': 1e-02,
                'dm**2': 1e-02,
                'centimeter**2': 1e-04,
                'cm2': 1e-04,
                'cm**2': 1e-04,
                'millimeter**2': 1e-06,
                'mm2': 1e-06,
                'mm**2': 1e-06}

    @staticmethod
    def volume_ud():
        """
        Return volume unit dimensions as dict
        :return: Volume unit dimensions as dict
        """
        return {'meter**3': 1.,
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
                'mm**3': 1e-09}

    @staticmethod
    def mass_ud():
        """
        Return mass unit dimensions as dict
        :return: Mass unit dimensions as dict
        """
        return {'kg': 1.,
                'g': 1e-03,
                'mg': 1e-06}

    @staticmethod
    def mole_ud():
        """
        Return mole unit dimensions as dict
        :return: Mole unit dimensions as dict
        """
        return {'kmol': 1.,
                'mol': 1e-03}

    @staticmethod
    def energy_ud():
        """
        Return energy unit dimensions as dict
        :return: Energy unit dimensions as dict
        """
        return {'J': 1.,
                'kJ': 1e03,
                'MJ': 1e06,
                'cal': 4.184,
                'kcal': 4.184 * 1e03,
                'kWh': 3.6 * 1e06,
                'Wh': 3.6 * 1e03}

    @staticmethod
    def power_ud():
        """
        Return power unit dimensions as dict
        :return: Power unit dimensions as dict
        """
        return {'W': 1.,
                'kW': 1e03,
                'MW': 1e06,
                'cal/s': 4.184,
                'kcal/s': 4.184 * 1e03,
                'cal/min': 4.184 / 60.,
                'kcal/min': 4.184 * 1e03 / 60.,
                'cal/h': 4.184 / 3600.,
                'kcal/h': 4.184 * 1e03 / 3600.}

    @staticmethod
    def viscosity_ud():
        """
        Return viscosity unit dimensions as dict
        :return: Viscosity unit dimensions as dict
        """
        return {'Pas': 1.,
                'cP': 1000.,
                'g/cm/s': 10.}

    @staticmethod
    def convert_from_to(to_unit_dimension):
        """
        Convert from to unit
        :param to_unit_dimension: Convert TO unit dimension
        :return: Converter FROM unit dimension
        """
        return {k: 1 / v for k, v in to_unit_dimension.items()}

    @staticmethod
    def create_fractional_ud(numerator, denominator):
        """
        Create fraction unit dimensions
        :param numerator: Numerator unit dimension
        :param denominator: Denominator unit dimension
        :return: to unit dimension dict, from unit dimension dict
        """
        to_fraction_unit_dimension = {}
        for kup, vup in numerator.items():
            for kdown, vdown in denominator.items():
                to_fraction_unit_dimension[kup + "/" + kdown] = vup / vdown

        return to_fraction_unit_dimension, UnitConverter.convert_from_to(to_fraction_unit_dimension)

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

    def convert_to_pascal_seconds(self, value, ud):
        """
        Convert to Pascal*seconds
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in Pascal*seconds
        """
        return self.converter(value, ud, 'Pas', self.to_pascal_seconds, self.from_pascal_seconds)

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

    def convert_to_kg_per_kmol(self, value, ud):
        """
        Convert to kg/kmol
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in kg/kmol
        """
        return self.converter(value,
                              ud,
                              'kg/kmol',
                              self.to_kilograms_per_kilomole,
                              self.from_kilograms_per_kilomole)

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

    def convert_to_square_meter_per_seconds(self, value, ud):
        """
        Convert to m2/s
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in m2/s
        """
        return self.converter(value,
                              ud,
                              'meter**2/s',
                              self.to_square_meter_per_seconds,
                              self.from_square_meter_per_seconds)

    def convert_to_meter_per_seconds(self, value, ud):
        """
        Convert to m/s
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in m/s
        """
        return self.converter(value,
                              ud,
                              'meter/s',
                              self.to_meter_per_seconds,
                              self.from_meter_per_seconds)

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

    def convert_to_joule_per_kmol_per_kelvin(self, value, ud):
        """
        Convert to J/kmol/K
        :param value: Value to be converted
        :param ud: Value initial unit dimension
        :return: Value in J/kmol/K
        """

        ud_split = ud.split("/")
        ud_clean = ud_split[0] + "/" + ud_split[1]

        return self.converter(value,
                              ud_clean,
                              'J/kmol',
                              self.to_joule_per_kilomole,
                              self.from_joule_per_kilomole)

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

    def convert_from_kg_per_kmol(self, value, ud):
        """
        Convert from kg/kmol
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimension
        """
        return self.converter(value,
                              'kg/kmol',
                              ud,
                              self.to_kilograms_per_kilomole,
                              self.from_kilograms_per_kilomole)

    def convert_from_kg_per_cubic_meter(self, value, ud):
        """
        Convert from kg/m3
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value,
                              'kg/meter**3',
                              ud,
                              self.to_kilograms_per_cubic_meter,
                              self.from_kilograms_per_cubic_meter)

    def convert_from_joule_per_kg_per_kelvin(self, value, ud):
        """
        Convert from J/kg/K
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimension
        """

        ud_split = ud.split("/")
        ud_clean = ud_split[0] + "/" + ud_split[1]

        return self.converter(value,
                              'J/kg',
                              ud_clean,
                              self.to_joule_per_kilograms,
                              self.from_joule_per_kilograms)

    def convert_from_joule_per_kmol_per_kelvin(self, value, ud):
        """
        Convert from J/kmol/K
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimension
        """

        ud_split = ud.split("/")
        ud_clean = ud_split[0] + "/" + ud_split[1]

        return self.converter(value,
                              'J/kmol',
                              ud_clean,
                              self.to_joule_per_kilomole,
                              self.from_joule_per_kilomole)

    def convert_from_watt_per_meter_per_kelvin(self, value, ud):
        """
        Convert from W/m/K
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimension
        """

        ud_split = ud.split("/")
        ud_clean = ud_split[0] + "/" + ud_split[1]

        return self.converter(value,
                              'W/meter',
                              ud_clean,
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

    def convert_from_pascal_seconds(self, value, ud):
        """
        Convert from Pascal*seconds
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final unit dimensions
        """
        return self.converter(value,
                              ud,
                              'Pas',
                              self.to_pascal_seconds,
                              self.from_pascal_seconds)

    def convert_from_cubic_meter_per_seconds(self, value, ud):
        """
        Convert from m3/s
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final
        """
        return self.converter(value,
                              'meter**3/s',
                              ud,
                              self.to_cubic_meter_per_seconds,
                              self.from_cubic_meter_per_seconds)

    def convert_from_square_meter_per_seconds(self, value, ud):
        """
        Convert from m2/s
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final
        """
        return self.converter(value,
                              'meter**2/s',
                              ud,
                              self.to_square_meter_per_seconds,
                              self.from_square_meter_per_seconds)

    def convert_from_meter_per_seconds(self, value, ud):
        """
        Convert from m/s
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final
        """
        return self.converter(value,
                              'meter/s',
                              ud,
                              self.to_meter_per_seconds,
                              self.from_meter_per_seconds)

    def convert_from_meter(self, value, ud):
        """
        Convert from meter
        :param value: Value to be converted
        :param ud: Value final unit dimension
        :return: Value in final ud
        """
        return self.converter(value, 'meter', ud, self.to_meter, self.from_meter)

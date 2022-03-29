from asali.utils.unit_converter import UnitConverter


class SolidMaterial:
    def __init__(self):
        """
        Class to describe solid material properties
        """
        self.uc = UnitConverter()
        self._density = 2300.0
        self._thermal_conductivity = 2.5
        self._specific_heat = 750.0

    def get_density(self):
        """
        Density getter
        :return: Density tolerance
        """
        return self._density

    def set_density(self, value):
        """
        Density setter
        :param value: Density
        :return:
        """
        self._density = self.uc.convert_to_kg_per_cubic_meter(value[0], value[1])

    # Creating a property object for Density
    density = property(get_density, set_density)

    def get_thermal_conductivity(self):
        """
        Thermal conductivity getter
        :return: Thermal conductivity
        """
        return self._thermal_conductivity

    def set_thermal_conductivity(self, value):
        """
        Thermal conductivity setter
        :param value: Thermal conductivity
        :return:
        """
        self._thermal_conductivity = self.uc.convert_to_watt_per_meter_per_kelvin(value[0], value[1])

    # Creating a property object for Thermal conductivity
    thermal_conductivity = property(get_thermal_conductivity, set_thermal_conductivity)

    def get_specific_heat(self):
        """
        Specific heat getter
        :return: Specific heat
        """
        return self._specific_heat

    def set_specific_heat(self, value):
        """
        Specific heat setter
        :param value: Specific heat
        :return:
        """
        self._specific_heat = self.uc.convert_to_joule_per_kg_per_kelvin(value[0], value[1])

    # Creating a property object for Absolute tolerance
    specific_heat = property(get_specific_heat, set_specific_heat)

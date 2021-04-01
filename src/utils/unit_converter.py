import numpy as np
import os

from pint import UnitRegistry


class UnitConverter:
    def __init__(self):
        self.unit_register = UnitRegistry(system='mks')
        self.unit_register.load_definitions(os.path.join(os.path.dirname(__file__), 'ud.asali'))
        self.unit_quality = self.unit_register.Quantity

    def convert_to_seconds(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('second').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_meter(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('meter').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_cubic_meter(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('meter**3').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_pascal(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('pascal').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_one_over_meter(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('1/meter').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_kg_per_seconds(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('kg/s').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_cubic_meter_per_seconds(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('meter**3/s').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_kelvin(self, value, ud):
        converter = self.unit_quality(0., self.unit_register.parse_expression(ud)).to('kelvin').magnitude

        if isinstance(value, float):
            return value + converter

        if isinstance(value, int):
            return value + converter

        return np.asarray(value) + converter

    def convert_to_kg_per_cubic_meter(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('kg/meter**3').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_joule_per_kg_per_kelvin(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('J/kg/degK').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_to_watt_per_meter_per_kelvin(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression(ud)).to('W/meter/degK').magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_from_seconds(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression('second')).to(ud).magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_from_cubic_meter(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression('meter**3')).to(ud).magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

    def convert_from_pascal(self, value, ud):
        converter = self.unit_quality(1., self.unit_register.parse_expression('pascal')).to(ud).magnitude

        if isinstance(value, float):
            return value * converter

        if isinstance(value, int):
            return value * converter

        return np.asarray(value) * converter

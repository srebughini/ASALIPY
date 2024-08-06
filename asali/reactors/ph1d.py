from abc import ABC

from asali.reactors.basic1d import Basic1DReactor
from asali.utils.input_parser import ReactorType


class PseudoHomogeneous1DReactor(Basic1DReactor, ABC):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing PseudoHomogeneous 1D reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)

        self.solution_parser.reactor_type = ReactorType.PSEUDOHOMOGENEOUSPFR
        self.n_v = self.gas.n_species + self.surf.n_species + 1
        self._setup.n_v = self.n_v

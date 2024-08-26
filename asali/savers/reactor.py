from asali.savers.batch_and_cstr import BatchAndCstrSaver
from asali.savers.het1d import Heterogeneous1DReactorSaver
from asali.savers.ph1d import PseudoHomogeneous1DReactorSaver
from asali.utils.input_parser import ReactorType


class ReactorSaver:
    def __init__(self, cls):
        """
        Class to save reactor results
        :param cls: Reactor class object
        """
        if cls.solution_parser.reactor_type not in [r for r in ReactorType]:
            raise Exception("ASALI::ERROR::Unknown class type: ", str(type(cls)))

        if not cls.solution_parser.is_solved:
            raise Exception("ASALI::ERROR::Nothing to plot, no solution found")

        self.saver = ReactorSaver.get_saver(cls)

    @staticmethod
    def get_saver(cls):
        """
        Function to get the correct saver based on Reactor Type
        :param cls: Reactor class object
        :return: Saver object
        """
        if cls.solution_parser.reactor_type == ReactorType.BATCH:
            return BatchAndCstrSaver(cls)

        if cls.solution_parser.reactor_type == ReactorType.CSTR:
            return BatchAndCstrSaver(cls)

        if cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            return PseudoHomogeneous1DReactorSaver(cls)

        if cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            return Heterogeneous1DReactorSaver(cls)

    def save_using_mass_fraction(self, file_path, species_names=None, coverage_names=None):
        """
        Saving output to file
        :param file_path: File path where to save the results
        :param species_names: List of species to be saved
        :param coverage_names: List of coverage to be saved
        :return:
        """
        self.saver.save_using_mass_fraction(file_path, species_names, coverage_names)

    def save_using_mole_fraction(self, file_path, species_names=None, coverage_names=None):
        """
        Saving output to file
        :param file_path: File path where to save the results
        :param species_names: List of species to be saved
        :param coverage_names: List of coverage to be saved
        :return:
        """
        self.saver.save_using_mole_fraction(file_path, species_names, coverage_names)

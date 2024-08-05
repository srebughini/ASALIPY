from enum import IntEnum

import warnings


class ReactorType(IntEnum):
    """
    Reactor type
    """
    BATCH = 0
    CSTR = 1
    PSEUDOHOMOGENEOUSPFR = 2
    HETEROGENEOUSPRF = 3


class ResolutionMethod(IntEnum):
    """
    Resolution method
    """
    STEADYSTATE = 0
    TRANSIENT = 1


class ReactorSection(IntEnum):
    """
    Reactor section
    """
    CIRCLE = 0
    SQUARE = 1
    TRIANGLE = 2


class ReactorModel(IntEnum):
    """
    Reactor model
    """
    TUBULAR = 0
    PACKEDBED = 1
    HONEYCOMB = 2


class InputParser:
    """
    Class to parse user input during setting up reactor models
    """

    @staticmethod
    def print_warning(message, verbosity):
        """
        Print warning of screen
        ::param message: Message to be printed
        :return:
        """
        warning_format = 'ASALI::WARNING::{}'
        if verbosity > 20:
            warnings.warn(warning_format.format(message))

    @staticmethod
    def raise_error(message):
        """
        Raise error
        :param message: Message to be printed
        :return:
        """
        error_format = 'ASALI::ERROR::{}'
        raise Exception(error_format.format(message))

    @staticmethod
    def true_parser(v):
        """
        Parser of True
        :param v: Variable to be parsed
        :return: Bool rapresenting the variable
        """
        if isinstance(v, str):
            return v.lower().strip() in ['yes', 'y', '1', 'true', 'on']

        if v is True:
            return True

        if v == 1:
            return True

    @staticmethod
    def resolution_parser(method):
        """
        Parser of the resolution method
        :param method: Resolution method as string
        :return: ResolutionMethod object
        """
        if method in ["steadystate", "ss", "SteadyState", "STEADYSTATE", "steady state", "Steady State"]:
            return ResolutionMethod.STEADYSTATE

        if method in ["transient", "tt", "Transient", "TRANSIENT"]:
            return ResolutionMethod.TRANSIENT

        InputParser.raise_error("Unknown resolution method")

    @staticmethod
    def section_parser(method):
        """
        Parser of the reactor section
        :param method: Reactor section as string
        :return: ReactorSection
        """
        if method in ["circle", "Circle", "CIRCLE"]:
            return ReactorSection.CIRCLE

        if method in ["square", "Square", "SQUARE"]:
            return ReactorSection.SQUARE

        if method in ["triangle", "Triangle", "TRIANGLE"]:
            return ReactorSection.TRIANGLE

        InputParser.raise_error("Unknown section shape")

from asali.reactors.batch import BatchReactor
from asali.reactors.cstr import CstrReactor
from asali.reactors.het1d_steady_state import SteadyStateHeterogeneous1DReactor
from asali.reactors.het1d_transient import TransientHeterogeneous1DReactor
from asali.reactors.ph1d_steady_state import SteadyStatePseudoHomogeneous1DReactor
from asali.reactors.ph1d_transient import TransientPseudoHomogeneous1DReactor
from asali.utils.input_parser import InputParser
from tests.basic_unit_test import BasicUnitTest
from tests.reactor_unit_test import ReactorUnitTest
from asali.utils.unit_converter import UnitConverter
from asali.utils.cantera_file_converter import CanteraFileConverter
from termcolor import colored

import argparse


def unit_converter():
    """
    Run unit test for UnitConverter
    :return:
    """
    uc = UnitConverter()
    ut = BasicUnitTest(uc, "utils.json")
    ut.check_all()
    print(" ")


def cantera_file_converter():
    """
    Run unit test for CanteraFileConverter
    :return:
    """
    cfc = CanteraFileConverter()
    ut = BasicUnitTest(cfc, "utils.json")
    ut.check_all()
    print(" ")


def input_parser():
    """
    Run unit test for InputParser
    :return:
    """
    ip = InputParser()
    ut = BasicUnitTest(ip, "utils.json")
    ut.check_all()
    print(" ")


def batch_reactor(is_local):
    """
    Run unit test for BatchReactor
    :param is_local: Bool that shows if the test is performed locally or not
    :return:
    """
    br = BatchReactor('tests/H2-O2-Rh.yaml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(br, "reactors.json", is_local)
    ut.check_all()
    print(" ")


def cstr_reactor(is_local):
    """
    Run unit test for CstrReactor
    :param is_local: Bool that shows if the test is performed locally or not
    :return:
    """
    cr = CstrReactor('tests/H2-O2-Rh.yaml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(cr, "reactors.json", is_local)
    ut.check_all()
    print(" ")


def steady_state_pseudohomogeneous_reactor(is_local):
    """
    Run unit test for SteadyStatePseudoHomogeneous1DReactor
    :param is_local: Bool that shows if the test is performed locally or not
    :return:
    """
    pfr = SteadyStatePseudoHomogeneous1DReactor('tests/H2-O2-Rh.yaml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(pfr, "reactors.json", is_local)
    ut.check_all()
    print(" ")


def transient_pseudohomogeneous_reactor(is_local):
    """
    Run unit test for TransientPseudoHomogeneous1DReactor
    :param is_local: Bool that shows if the test is performed locally or not
    :return:
    """
    pfr = TransientPseudoHomogeneous1DReactor('tests/H2-O2-Rh.yaml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(pfr, "reactors.json", is_local)
    ut.check_all()
    print(" ")


def steady_state_heterogeneous_reactor(is_local):
    """
    Run unit test for SteadyStateHeterogeneous1DReactor
    :param is_local: Bool that shows if the test is performed locally or not
    :return:
    """
    pfr = SteadyStateHeterogeneous1DReactor('tests/H2-O2-Rh.yaml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(pfr, "reactors.json", is_local)
    ut.check_all()
    print(" ")

def transient_heterogeneous_reactor(is_local):
    """
    Run unit test for TransientHeterogeneous1DReactor
    :param is_local: Bool that shows if the test is performed locally or not
    :return:
    """
    pfr = TransientHeterogeneous1DReactor('tests/H2-O2-Rh.yaml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(pfr, "reactors.json", is_local)
    ut.check_all()
    print(" ")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--class', dest='class_to_check', help='Class to be check')
    parser.add_argument('--local', action='store_true', dest='is_local', help='Is the test performed locally')
    args = parser.parse_args()

    if args.class_to_check is None:
        class_to_check = "all"
    else:
        class_to_check = args.class_to_check.lower()

    if class_to_check == "all":
        input_parser()
        unit_converter()
        cantera_file_converter()
        batch_reactor(args.is_local)
        cstr_reactor(args.is_local)
        steady_state_pseudohomogeneous_reactor(args.is_local)
        transient_pseudohomogeneous_reactor(args.is_local)
        steady_state_heterogeneous_reactor(args.is_local)
        transient_heterogeneous_reactor(args.is_local)
    elif class_to_check == "unitconverter":
        unit_converter()
    elif class_to_check == "batchreactor":
        batch_reactor(args.is_local)
    elif class_to_check == "cstrreactor":
        cstr_reactor(args.is_local)
    elif class_to_check == "steadystatepseudohomogeneous1dreactor":
        steady_state_pseudohomogeneous_reactor(args.is_local)
    elif class_to_check == "transientpseudohomogeneous1dreactor":
        transient_pseudohomogeneous_reactor(args.is_local)
    elif class_to_check == "steadystateheterogeneous1dreactor":
        steady_state_heterogeneous_reactor(args.is_local)
    elif class_to_check == "transientheterogeneous1dreactor":
        transient_heterogeneous_reactor(args.is_local)
    elif class_to_check == "canterafileconverter":
        cantera_file_converter()
    elif class_to_check == "inputparser":
        input_parser()
    else:
        print(colored("ASALI::ERROR::UNKNOWN CLASS", "red"))

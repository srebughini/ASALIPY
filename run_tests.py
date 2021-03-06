from asali.reactors.batch import BatchReactor
from asali.reactors.cstr import CstrReactor
from asali.reactors.het1d import Heterogeneous1DReactor
from asali.reactors.ph1d import PseudoHomogeneous1DReactor
from tests.basic_unit_test import BasicUnitTest
from tests.reactor_unit_test import ReactorUnitTest
from asali.utils.unit_converter import UnitConverter
from termcolor import colored

import argparse


def unit_converter():
    uc = UnitConverter()
    ut = BasicUnitTest(uc, "utils.json")
    ut.check_all()
    print(" ")


def batch_reactor():
    br = BatchReactor('tests/H2-O2-Rh.xml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(br, "reactors.json")
    ut.check_all()
    print(" ")


def cstr_reactor():
    cr = CstrReactor('tests/H2-O2-Rh.xml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(cr, "reactors.json")
    ut.check_all()
    print(" ")


def pseudohomogeneous_reactor():
    pfr = PseudoHomogeneous1DReactor('tests/H2-O2-Rh.xml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(pfr, "reactors.json")
    ut.check_all()
    print(" ")


def heterogeneous_reactor():
    pfr = Heterogeneous1DReactor('tests/H2-O2-Rh.xml', 'gas', 'Rh_surface')
    ut = ReactorUnitTest(pfr, "reactors.json")
    ut.check_all()
    print(" ")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--class', dest='class_to_check', help='Class to be check')
    args = parser.parse_args()

    if args.class_to_check is None:
        class_to_check = "all"
    else:
        class_to_check = args.class_to_check.lower()

    if class_to_check == "all":
        unit_converter()
        batch_reactor()
        cstr_reactor()
        pseudohomogeneous_reactor()
        heterogeneous_reactor()
    elif class_to_check == "unitconverter":
        unit_converter()
    elif class_to_check == "batchreactor":
        batch_reactor()
    elif class_to_check == "cstrreactor":
        cstr_reactor()
    elif class_to_check == "pseudohomogeneous1dreactor":
        pseudohomogeneous_reactor()
    elif class_to_check == "heterogeneous1dreactor":
        heterogeneous_reactor()
    else:
        print(colored("ASALI::ERROR::Unknown class", "red"))

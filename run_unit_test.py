from asali.reactors.batch import BatchReactor
from asali.reactors.cstr import CstrReactor
from asali.reactors.het1d import Heterogeneous1DReactor
from asali.reactors.ph1d import PseudoHomogeneous1DReactor
from tests.basic_unit_test import BasicUnitTest
from tests.reactor_unit_test import ReactorUnitTest
from asali.utils.unit_converter import UnitConverter


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
    unit_converter()
    batch_reactor()
    cstr_reactor()
    pseudohomogeneous_reactor()
    heterogeneous_reactor()

from cantera import ck2cti
from cantera import ctml_writer
from cantera import ctml2yaml
from cantera import ck2yaml
from cantera import cti2yaml

import os


class CanteraFileConverter:

    @staticmethod
    def check_if_file_exist(file_path):
        if not os.path.isfile(file_path):
            raise Exception(file_path + " does not exist!!!")

    @staticmethod
    def replace_extension(file_path, new_extension):
        _, extension = os.path.splitext(file_path)
        if extension is None:
            return file_path + new_extension
        else:
            return file_path.replace(extension, new_extension)

    @staticmethod
    def check_file_extension(file_path, desired_extension):
        if desired_extension not in file_path:
            raise Exception(file_path, " wrong file format!!")

    @staticmethod
    def parse_chemkin_inputs(kinetic_file_path, thermodynamic_file_path, transport_file_path, surface_file_path,
                             output_file_path, output_extension):
        CanteraFileConverter.check_if_file_exist(kinetic_file_path)
        CanteraFileConverter.check_if_file_exist(transport_file_path)
        CanteraFileConverter.check_if_file_exist(thermodynamic_file_path)

        if output_file_path is None:
            output_file_path = CanteraFileConverter.replace_extension(kinetic_file_path, output_extension)
        else:
            output_file_path = CanteraFileConverter.replace_extension(output_file_path, output_extension)

        if surface_file_path is not None:
            CanteraFileConverter.check_if_file_exist(surface_file_path)
            return ["--input=" + kinetic_file_path,
                    "--transport=" + transport_file_path,
                    "--thermo=" + thermodynamic_file_path,
                    "--surface=" + surface_file_path,
                    "--output=" + output_file_path]
        else:
            return ["--input=" + kinetic_file_path,
                    "--transport=" + transport_file_path,
                    "--thermo=" + thermodynamic_file_path,
                    "--output=" + output_file_path]

    @staticmethod
    def from_chemkin_to_cti(kinetic_file_path, thermodynamic_file_path, transport_file_path, surface_file_path=None,
                            output_file_path=None):
        ck2cti.main(CanteraFileConverter.parse_chemkin_inputs(kinetic_file_path,
                                                              thermodynamic_file_path,
                                                              transport_file_path,
                                                              surface_file_path,
                                                              output_file_path,
                                                              ".cti"))

    @staticmethod
    def from_chemkin_to_yaml(kinetic_file_path, thermodynamic_file_path, transport_file_path, surface_file_path=None,
                             output_file_path=None):

        ck2cti.main(CanteraFileConverter.parse_chemkin_inputs(kinetic_file_path,
                                                              thermodynamic_file_path,
                                                              transport_file_path,
                                                              surface_file_path,
                                                              output_file_path,
                                                              ".cti"))

        CanteraFileConverter.from_cti_to_yaml(CanteraFileConverter.replace_extension(output_file_path, ".cti"),
                                              CanteraFileConverter.replace_extension(output_file_path, ".yaml"))
        """ TODO
        ck2yaml.main(CanteraFileConverter.parse_chemkin_inputs(kinetic_file_path,
                                                               thermodynamic_file_path,
                                                               transport_file_path,
                                                               surface_file_path,
                                                               output_file_path,
                                                               ".yaml"))
        """

    @staticmethod
    def from_cti_to_xml(input_file_path, output_file_path):
        CanteraFileConverter.check_if_file_exist(input_file_path)
        CanteraFileConverter.check_file_extension(input_file_path, ".cti")

        ctml_writer.convert(filename=input_file_path,
                            outName=CanteraFileConverter.replace_extension(output_file_path, ".xml"))

    @staticmethod
    def from_cti_to_yaml(input_file_path, output_file_path):
        CanteraFileConverter.check_if_file_exist(input_file_path)
        CanteraFileConverter.check_file_extension(input_file_path, ".cti")

        cti2yaml.convert(input_file_path,
                         CanteraFileConverter.replace_extension(output_file_path, ".yaml"))

    @staticmethod
    def from_xml_to_yaml(input_file_path, output_file_path):
        CanteraFileConverter.check_if_file_exist(input_file_path)
        CanteraFileConverter.check_file_extension(input_file_path, ".xml")

        ctml2yaml.convert(input_file_path,
                          CanteraFileConverter.replace_extension(output_file_path, ".yaml"))

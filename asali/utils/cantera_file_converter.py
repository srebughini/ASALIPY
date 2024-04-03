from cantera import ctml2yaml, ck2yaml
from cantera import cti2yaml

import os


class CanteraFileConverter:
    """
    Class collecting all the utils for convert Cantera input files
    """

    @staticmethod
    def check_if_file_exist(file_path):
        """
        Check if file exist, if not raise and Exception
        :param file_path: File path
        :return:
        """
        if not os.path.isfile(file_path):
            raise Exception(file_path + " does not exist!!!")

    @staticmethod
    def replace_extension(file_path, new_extension):
        """
        Replace file extension
        :param file_path: File path
        :param new_extension: New extension
        :return: File path with the new extension
        """
        _, extension = os.path.splitext(file_path)
        if extension is None:
            return file_path + new_extension

        return file_path.replace(extension, new_extension)

    @staticmethod
    def check_file_extension(file_path, desired_extension):
        """
        Check file extension, if not raise and Exception
        :param file_path: File path
        :param desired_extension: Extensions
        :return:
        """
        _, extension = os.path.splitext(file_path)
        if desired_extension not in extension:
            raise Exception(file_path, " wrong file format!!")

    @staticmethod
    def parse_chemkin_inputs(kinetic_file_path, thermodynamic_file_path, transport_file_path, surface_file_path,
                             output_file_path, output_extension):
        """
        Parse of CHEMKIN input file
        :param kinetic_file_path: CHEMKIN kinetic file path
        :param thermodynamic_file_path: CHEMKIN thermodynamic file path
        :param transport_file_path: CHEMKIN transport file path
        :param surface_file_path: CHEMKIN surface kinetic file path
        :param output_file_path: Output file path
        :param output_extension: Output file extension
        :return: List of commands to convert CHEMKIN file, Output file path with correct extension
        """
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
                    "--output=" + output_file_path,
                    "--quiet",
                    "--permissive"], output_file_path

        return ["--input=" + kinetic_file_path,
                "--transport=" + transport_file_path,
                "--thermo=" + thermodynamic_file_path,
                "--output=" + output_file_path,
                "--quiet",
                "--permissive"], output_file_path

    @staticmethod
    def from_chemkin_to_yaml(kinetic_file_path, thermodynamic_file_path, transport_file_path, surface_file_path=None,
                             output_file_path=None):
        """
        Convert from CHEMKIN input files to Cantera .yaml file
        :param kinetic_file_path: CHEMKIN kinetic file path
        :param thermodynamic_file_path: CHEMKIN thermodynamic file path
        :param transport_file_path: CHEMKIN transport file path
        :param surface_file_path: CHEMKIN surface kinetic file path
        :param output_file_path: Output file path
        :return: Cantera .yaml output file path
        """

        input_list, _ = CanteraFileConverter.parse_chemkin_inputs(kinetic_file_path,
                                                                  thermodynamic_file_path,
                                                                  transport_file_path,
                                                                  surface_file_path,
                                                                  output_file_path,
                                                                  ".yaml")
        ck2yaml.main(input_list)
        file_path = CanteraFileConverter.replace_extension(output_file_path, ".yaml")
        return file_path

    @staticmethod
    def from_cti_to_yaml(input_file_path, output_file_path):
        """
        Convert from Cantera .cti file to Cantera .yaml file
        :param input_file_path: Cantera .cti input file path
        :param output_file_path: Cantera .yaml output file path
        :return: Cantera .yaml output file path
        """
        CanteraFileConverter.check_if_file_exist(input_file_path)
        CanteraFileConverter.check_file_extension(input_file_path, ".cti")

        file_path = CanteraFileConverter.replace_extension(output_file_path, ".yaml")
        cti2yaml.convert(input_file_path, file_path)
        return file_path

    @staticmethod
    def from_xml_to_yaml(input_file_path, output_file_path):
        """
        Convert from Cantera .xml file to Cantera .yaml file
        :param input_file_path: Cantera .xml input file path
        :param output_file_path: Cantera .yaml output file path
        :return: Cantera .yaml output file path
        """
        CanteraFileConverter.check_if_file_exist(input_file_path)
        CanteraFileConverter.check_file_extension(input_file_path, ".xml")

        file_path = CanteraFileConverter.replace_extension(output_file_path, ".yaml")
        ctml2yaml.convert(input_file_path, file_path)
        return file_path

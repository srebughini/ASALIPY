import os

from asali.utils.cantera_file_converter import CanteraFileConverter

if __name__ == "__main__":
    # Convert from CHEMKIN format to YAML format
    CanteraFileConverter.from_chemkin_to_yaml(kinetic_file_path=os.path.join("files", "kinetic.kin"),
                                              thermodynamic_file_path=os.path.join("files", "thermo.tdc"),
                                              transport_file_path=os.path.join("files", "transport.tra"),
                                              surface_file_path=os.path.join("files", "surface.sur"),
                                              output_file_path=os.path.join("files", "output_v3.yaml"))

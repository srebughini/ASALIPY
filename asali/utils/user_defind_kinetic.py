import json
import re

import sympy as sp
import numpy as np

from enum import Enum

from asali.utils.input_parser import InputParser


class UserDefinedKineticKeys(Enum):
    """
    Enum class to handle the .json keys
    """
    REACTIONS = "reactions"
    SPECIES = "species"
    TYPE = "type"
    FORMULA = "formula"
    RATE = "rate"
    RATEUNITS = "rate_units"
    SPECIEUNITS = "specie_units"
    ID = "id"
    TEMPERATURE = "temperature"


class ReactionType(Enum):
    """
    Enum class to handle reaction types
    """
    HOMO = "homogeneous"
    HET = "heterogeneous"


class CompositionType(Enum):
    """
    Enum class to handle composition type
    """
    MOLE = "mole_fraction"
    MASS = "mass_fraction"


class UserDefinedKinetic:
    def __init__(self):
        """
        Class to handle user defined kinetic
        """
        self._species_list = None
        self._file_path = None
        self._coefficients_matrix = None
        self._homogeneous_reaction_array = None
        self._heterogeneous_reaction_array = None
        self._reaction_rate_array = None
        self._cantera_idx_to_udk_idx_dict = None
        self._is_set = False

    def get_file_path(self):
        """
        Get chemistry file path
        :return: Chemistry file path
        """
        return self._file_path

    def set_file_path(self, value):
        """
        Set chemistry file path
        :param value: Chemistry file path
        :return
        """
        self._is_set = True
        self._file_path = value

    # Creating a property object for chemistry file path
    file_path = property(get_file_path, set_file_path)

    def get_is_set(self):
        """
        Get bool to check if the user defined kinect is set
        :return: If True the udk is et, if False is not
        """
        return self._file_path

    def set_is_set(self, value):
        """
        Set bool to check if the user defined kinect is set
        :param value: Bool to set user defined kinetic
        :return
        """
        self._is_set = value

    # Creating a property object for 'is set' bool
    is_set = property(get_is_set, set_is_set)

    @staticmethod
    def load_file(file_path):
        """
        Load file
        :param file_path: File path for user defined kinetic .json
        :return: File as dict
        """
        with open(file_path, 'r') as file:
            data = json.load(file)

        return data

    @staticmethod
    def parse_species(term):
        """
        Parse species from reactant/product of a chemical formula
        :param term: Reactant or Product part of chemical formular
        :return: Stoichiometric coefficient as list, Species as list
        """
        match = re.match(r'(\d*)\s*([A-Za-z0-9_]+)', term)
        if match:
            coeff = int(match.group(1)) if match.group(1) else 1
            species = match.group(2)
            return coeff, species
        else:
            raise InputParser.raise_error(f"Invalid term: {term}")

    @staticmethod
    def parse_reaction_rate(reaction_rate_as_str, species_list, temperature_as_variable):
        """
        Parse the reaction rate into a callable function
        :param reaction_rate_as_str: Reaction rate as str
        :param species_list: List of species
        :param temperature_as_variable: Temperature variable name
        :return: Reaction rate as callable function
        """
        # Replace species concentrations with symbolic variables
        for species in species_list:
            reaction_rate_as_str = reaction_rate_as_str.replace(f'[{species}]', f'{species}')
        # Convert the expression into a sympy expression
        reaction_rate_expr = sp.sympify(reaction_rate_as_str)

        # Add temperature variable
        species_and_temperature_list = species_list.copy()
        species_and_temperature_list.append(temperature_as_variable)

        # Create a lambda function for the rate expression
        reaction_rate_function = sp.lambdify([sp.Symbol(v) for v in species_and_temperature_list],
                                             reaction_rate_expr,
                                             'numpy')
        return reaction_rate_function

    @staticmethod
    def extract_stoichiometric_coefficients(formula_as_str):
        """
        Extract stoichiometric coefficients from formula
        :param formula_as_str: Chemical formula
        :return: Stoichiometric coefficient as dict
        """
        coefficients = {}
        # Split the reaction into reactants and products
        reactants, products = formula_as_str.split('->')
        # Process the reactants
        for reactant in reactants.split('+'):
            coeff, species = UserDefinedKinetic.parse_species(reactant.strip())
            coefficients[species] = -coeff  # Negative for reactants
        # Process the products
        for product in products.split('+'):
            coeff, species = UserDefinedKinetic.parse_species(product.strip())
            coefficients[species] = coeff  # Positive for products
        return coefficients

    @staticmethod
    def validate_file(data):
        """
        Validate and parse .json file
        :param data: .json file in dict format
        :return: Species as list,
                 Stoichiometric coefficient as dict,
                 Reaction rates as callable functions,
                 Reaction types as Enum
        """
        coefficients = {}
        reaction_rates = {}
        reaction_type = {}

        # Check the species key
        if UserDefinedKineticKeys.SPECIES.value not in data or not isinstance(
                data[UserDefinedKineticKeys.SPECIES.value], list):
            raise InputParser.raise_error(f"Invalid JSON: '{UserDefinedKineticKeys.SPECIES.value}' must be a list.")

        # Check the reactions key
        if UserDefinedKineticKeys.REACTIONS.value not in data or not isinstance(
                data[UserDefinedKineticKeys.REACTIONS.value], list):
            raise InputParser.raise_error(f"Invalid JSON: '{UserDefinedKineticKeys.REACTIONS.value}' must be a list.")

        # Check the temperature key
        if UserDefinedKineticKeys.TEMPERATURE.value not in data or not isinstance(
                data[UserDefinedKineticKeys.TEMPERATURE.value], str):
            raise InputParser.raise_error(f"Invalid JSON: '{UserDefinedKineticKeys.TEMPERATURE.value}' must be a list.")

        # Check id key for each reaction
        for reaction in data[UserDefinedKineticKeys.REACTIONS.value]:
            if UserDefinedKineticKeys.ID.value not in reaction:
                raise InputParser.raise_error(
                    f"Invalid JSON: Each reaction must have a '{UserDefinedKineticKeys.ID.value}'.")

        # Check keys for each reaction
        for reaction in data[UserDefinedKineticKeys.REACTIONS.value]:
            id = reaction[UserDefinedKineticKeys.ID.value]

            for key in [UserDefinedKineticKeys.FORMULA,
                        UserDefinedKineticKeys.RATE,
                        UserDefinedKineticKeys.SPECIEUNITS,
                        UserDefinedKineticKeys.RATEUNITS,
                        UserDefinedKineticKeys.TYPE]:
                if key.value not in reaction:
                    raise InputParser.raise_error(
                        f"Invalid JSON: Reaction {id} must have an '{UserDefinedKineticKeys.FORMULA.value}'.")

            # Validate the reaction expression
            try:
                coefficients[id] = UserDefinedKinetic.extract_stoichiometric_coefficients(
                    reaction[UserDefinedKineticKeys.FORMULA.value])
            except Exception as e:
                raise InputParser.raise_error(
                    f"Invalid {UserDefinedKineticKeys.FORMULA.value} in reaction {id}: {e}")

            # Validate the reaction rate
            try:
                reaction_rates[id] = UserDefinedKinetic.parse_reaction_rate(reaction[UserDefinedKineticKeys.RATE.value],
                                                                            data[UserDefinedKineticKeys.SPECIES.value],
                                                                            data[UserDefinedKineticKeys.TEMPERATURE.value])
            except Exception as e:
                raise InputParser.raise_error(f"Invalid reaction rate in reaction {id}: {e}")

            # Extract reaction type
            try:
                reaction_type[reaction[UserDefinedKineticKeys.ID.value]] = ReactionType(
                    reaction[UserDefinedKineticKeys.TYPE.value])
            except Exception as e:
                raise InputParser.raise_error(
                    f"Invalid {UserDefinedKineticKeys.TYPE.value} in reaction {id}: it can be only homogeneous/heterogeneous")

            # Extract composition type
            try:
                composition_type = CompositionType(reaction[UserDefinedKineticKeys.SPECIEUNITS.value])
                if composition_type != CompositionType.MASS:
                    InputParser.raise_error(
                        f"Invalid {UserDefinedKineticKeys.SPECIEUNITS.value} in reaction {id}: the only accepted value is mass_fraction")
            except:
                raise InputParser.raise_error(
                    f"Invalid {UserDefinedKineticKeys.SPECIEUNITS.value} in reaction {id}: the only accepted value is mass_fraction")

            # Validate unit dimensions
            try:
                reaction_rate_ud = reaction[UserDefinedKineticKeys.RATEUNITS.value]
                if reaction_type[id] == ReactionType.HOMO:
                    target_ud = "kmol/m3/s"
                else:
                    target_ud = "kmol/m2/s"

                if reaction_rate_ud != target_ud:
                    InputParser.raise_error(
                        f"Invalid {UserDefinedKineticKeys.RATEUNITS.value} in reaction {id}: the only accepted unit dimension is {target_ud}")
            except:
                raise InputParser.raise_error(
                    f"Invalid {UserDefinedKineticKeys.RATEUNITS.value} in reaction {id}: the only accepted unit dimension is {target_ud}")

        return data[UserDefinedKineticKeys.SPECIES.value], coefficients, reaction_rates, reaction_type

    def _convert_from_dict_to_matrix(self,
                                     gas,
                                     species_list,
                                     coefficients_dict,
                                     reaction_rates_dict,
                                     reaction_type_dicts):
        """
        Convert dic extracted from the .json in to matrix and vectors to better improve code performance
        :param gas: Cantera Solution object
        :param species_list: User defined kinetic species list
        :param coefficients_dict: Stoichiometric coefficient as dict
        :param reaction_rates_dict: Reaction rates as callable functions
        :param reaction_type_dicts: Reaction types as Enum
        :return:
        """
        self._coefficients_matrix = np.zeros([gas.n_species, len(coefficients_dict.keys())], dtype=np.float64)
        self._homogeneous_reaction_array = np.zeros([len(coefficients_dict.keys())], dtype=np.float64)
        self._heterogeneous_reaction_array = np.zeros([len(coefficients_dict.keys())], dtype=np.float64)
        self._reaction_rate_array = [None] * len(coefficients_dict.keys())

        for i, r in enumerate(coefficients_dict.keys()):
            self._reaction_rate_array[i] = reaction_rates_dict[r]

            if reaction_type_dicts[r] == ReactionType.HOMO:
                self._homogeneous_reaction_array[i] = 1.0

            if reaction_type_dicts[r] == ReactionType.HET:
                self._heterogeneous_reaction_array[i] = 1.0

            for j, n in enumerate(gas.species_names):
                if n in species_list:
                    self._coefficients_matrix[j, i] = coefficients_dict[r][n]

        self._cantera_idx_to_udk_idx_dict = {}
        for i, n in enumerate(self._species_list):
            idx = gas.species_names.index(n)
            self._cantera_idx_to_udk_idx_dict[idx] = i

    def _convert_mass_fraction_to_udk(self, gas):
        """
        Convert mass fraction from Cantera to User Defined Kinetic
        :param gas: Cantera Solution object
        :return: Mass fraction in User Defined Kinetic format
        """
        udk_mass_fraction = np.zeros([len(self._species_list)], dtype=np.float64)
        for cantera_idx, udk_idx in self._cantera_idx_to_udk_idx_dict.items():
            udk_mass_fraction[udk_idx] = gas.Y[cantera_idx]

        return udk_mass_fraction

    def _calculate_reaction_rates(self, gas):
        """
        Estimate reaction rates for each reaction
        :param gas: Cantera Solution object
        :return: Array of reaction rate for each reaction in kmol/m2/s or kmol/m3/s
        """
        udk_mass_fraction = self._convert_mass_fraction_to_udk(gas)
        udk_mass_fraction_and_temperature = np.append(udk_mass_fraction, gas.T)
        return np.asarray([r(*udk_mass_fraction_and_temperature) for r in self._reaction_rate_array])

    def _calculate_delta_enthalpy(self, gas):
        """
        Estimate delta enthalpy for each reaction
        :param gas: Cantera Solution object
        :return: Array of delta enthalpy for each reaction in J/kmol
        """
        return np.dot(self._coefficients_matrix.T, gas.partial_molar_enthalpies)

    def load_and_validate(self, gas):
        """
        Load and validate the .json file
        :param gas: Cantera Solution object
        :return:
        """
        data = UserDefinedKinetic.load_file(self._file_path)
        self._species_list, coefficients_dict, reaction_rates_dict, reaction_type_dicts = UserDefinedKinetic.validate_file(
            data)

        self._convert_from_dict_to_matrix(gas,
                                          self._species_list,
                                          coefficients_dict,
                                          reaction_rates_dict,
                                          reaction_type_dicts)

    def get_homogeneous_reaction_rates(self, gas):
        """
        Get homogeneous reaction rates for each species [kmol/m3/s]
        :param gas: Cantera Solution object
        :return: Homogeneous reaction rates for each species
        """
        reaction_rate_array = np.dot(self._coefficients_matrix,
                                     np.multiply(self._calculate_reaction_rates(gas),
                                                 self._homogeneous_reaction_array))

        return reaction_rate_array

    def get_heterogeneous_reaction_rates(self, gas):
        """
        Get heterogeneous reaction rates for each species [kmol/m2/s]
        :param gas: Cantera Solution object
        :return: Homogeneous reaction rates for each species
        """
        reaction_rate_array = np.dot(self._coefficients_matrix,
                                     np.multiply(self._calculate_reaction_rates(gas),
                                                 self._heterogeneous_reaction_array))

        return reaction_rate_array

    def get_homogeneous_heat_of_reaction(self, gas):
        """
        Get homogeneous heat of reaction in J/m3/s
        :param gas: Cantera Solution object
        :return: Float representing the heat of reaction
        """
        return -np.dot(np.multiply(self._calculate_reaction_rates(gas),
                                   self._homogeneous_reaction_array), self._calculate_delta_enthalpy(gas))

    def get_heterogeneous_heat_of_reaction(self, gas):
        """
        Get heterogeneous heat of reaction in J/m3/s
        :param gas: Cantera Solution object
        :return: Float representing the heat of reaction
        """
        return -np.dot(np.multiply(self._calculate_reaction_rates(gas),
                                   self._heterogeneous_reaction_array), self._calculate_delta_enthalpy(gas))

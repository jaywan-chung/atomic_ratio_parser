"""
    This module provides :py:func:`parse_atomic_ratio` function to parse atomic ratios from a chemical formula.

    :copyright: © 2025 by Jaywan Chung
    :license: MIT
"""
from decimal import Decimal
import re
import unittest

CHEMICAL_SYMBOLS = (
     'H', 'He', 'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',  'K', 'Ca',
    'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',  'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te',  'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
)
CHEMICAL_SYMBOL_PATTERN = f"({'|'.join(sorted(CHEMICAL_SYMBOLS, reverse=True))})"
POSITIVE_FLOAT_PATTERN = r'[0-9]+(\.[0-9]+)?'

_prog_chemical_symbol = re.compile(CHEMICAL_SYMBOL_PATTERN)
_prog_positive_float = re.compile(POSITIVE_FLOAT_PATTERN)


def parse_atomic_ratio(string):
    """Parse atomic ratios from the string.
    Only the first chemical formula is parsed if there are many.
    If there is no valid chemical formula in the string, return an empty dictionary.

    :param string: A string containing a chemical formula
    :type string: str
    :return: A dictionary whose key is chemical symbol and the value is the corresponding atomic ratio
    :rtype: dict

    Examples:
        >>> from atomic_ratio_parser import parse_atomic_ratio
        >>> ratio_dict = parse_atomic_ratio('Chemical formula for aluminium sulfate is Al2(SO4)3')
        >>> assert(ratio_dict == {'Al': Decimal('2'), 'S': Decimal('3'), 'O': Decimal('12')})
        >>> ratio_dict = parse_atomic_ratio('Pb0.95Na0.04Te, Bi2Te3')  # process only the first formula
        >>> assert(ratio_dict == {'Pb': Decimal('0.95'), 'Na': Decimal('0.04'), 'Te': Decimal('1')})
    """
    # remove non-necessary characters
    string = string.replace('_', '')
    string = string.replace('{', '')
    string = string.replace('}', '')
    string = string.replace(',', ' ')

    result = {}
    for token in string.split():
        try:
            # expand the formula and parse it
            expanded_formula = _get_expanded_chemical_formula(token)
            result = _parse_atomic_ratio_from_expanded_chemical_formula(expanded_formula)
        except ValueError:  # not a chemical formula
            continue
        break  # stop if a chemical formula is found
    return result


def _get_expanded_chemical_formula(chemical_formula):
    """Return an expanded chemical formula that has no parentheses.

    *Warning* The result can have the same chemical symbol several times.

    :param chemical_formula: A chemical formula that has no non-necessary symbols like `_`, `{`, `}`.
    :type chemical_formula: str
    :return: An expanded chemical formula
    :rtype: str

    :raises ValueError: When the `chemical_formula` is not a valid chemical formula.
    """
    prev_formula = None
    formula = chemical_formula

    while formula != prev_formula:  # make the chemical formula as simple as possible
        prev_formula = formula
        formula = _expand_innermost_atomic_ratio(formula)

    return formula


def _expand_innermost_atomic_ratio(chemical_formula):
    """Remove the first, innermost parentheses of the chemical formula.

    :param chemical_formula: A chemical formula that has no non-necessary symbols like `_`, `{`, `}`.
    :type chemical_formula: str
    :return: The chemical formula with the innermost formula expanded.
    :rtype: str

    :raises ValueError: When the `chemical_formula` is not a valid chemical formula.
    """
    # find the innermost parentheses
    right_parenthesis_idx = chemical_formula.find(')')
    left_parenthesis_idx = chemical_formula.rfind('(', 0, right_parenthesis_idx)
    if right_parenthesis_idx < 0:
        if left_parenthesis_idx < 0:  # if there is no parentheses, nothing to do
            return chemical_formula
        else:
            raise ValueError('Invalid chemical formula: missing right parenthesis')
    if left_parenthesis_idx < 0:
        raise ValueError('Invalid chemical formula: missing left parenthesis')

    ratio_match = _prog_positive_float.match(chemical_formula[right_parenthesis_idx+1:])
    parentheses_ratio = Decimal('1')
    if ratio_match:
        parentheses_ratio = Decimal(ratio_match.group())

    ratio_dict = _parse_atomic_ratio_from_expanded_chemical_formula(
        chemical_formula[left_parenthesis_idx+1:right_parenthesis_idx])
    for chemical_symbol in ratio_dict.keys():
        ratio_dict[chemical_symbol] *= parentheses_ratio

    expanded_str = _convert_ratio_dict_to_str(ratio_dict)

    return chemical_formula[:left_parenthesis_idx] + expanded_str \
        + chemical_formula[right_parenthesis_idx+len(ratio_match.group())+1:]


def _convert_ratio_dict_to_str(ratio_dict):
    """Convert a dictionary of atomic ratios into a string.
    Chemical symbols are sorted by atomic number, and invalid chemical symbols are ignored.

    :param ratio_dict: A dictionary of atomic ratios
    :type ratio_dict: dict
    :return: A chemical formula
    :rtype: str
    """
    result = ''
    for chemical_symbol in CHEMICAL_SYMBOLS:
        atomic_ratio = ratio_dict.get(chemical_symbol)
        if atomic_ratio:
            result += chemical_symbol + str(atomic_ratio)

    return result


def _parse_atomic_ratio_from_expanded_chemical_formula(expanded_chemical_formula):
    """Return a dictionary whose key is an atomic symbol and the value is the corresponding atomic ratio.
    The input must be a simple, expanded chemical formula that has no parentheses and no non-necessary symbols like `_`.
    For example, `Bi2Te3` is a simple, expanded chemical formula but `(BiTe)2Te` and `Bi_2Te_3` are not.

    :param expanded_chemical_formula: A chemical formula that has no parentheses and no non-necessary symbols.
    :type expanded_chemical_formula: str
    :return: A dictionary of atomic symbol-atomic ratio pairs.
    :rtype: dict

    :raises ValueError: When the `simple_chemical_formula` is not a simple, expanded chemical formula.
    """
    formula = expanded_chemical_formula
    result = {}

    while len(formula) > 0:
        # find chemical symbol
        match = _prog_chemical_symbol.match(formula)
        if not match:  # chemical symbol is not found
            break
        chemical_symbol = match.group()
        formula = formula[len(chemical_symbol):]

        # find atomic ratio
        match = _prog_positive_float.match(formula)
        atomic_ratio = Decimal('1')
        if match:  # atomic ratio is stated
            atomic_ratio_str = match.group()
            atomic_ratio = Decimal(atomic_ratio_str)
            formula = formula[len(atomic_ratio_str):]

        result[chemical_symbol] = result.get(chemical_symbol, Decimal('0')) + atomic_ratio

    if len(formula) > 0:
        raise ValueError(f'{repr(expanded_chemical_formula)} is not a simple, expanded chemical formula')

    return result


class AtomicRatioParserTest(unittest.TestCase):
    def test_parse_atomic_ratio(self):
        input_and_output_pairs = (
            ('(PbTe)0.7(PbS)0.3', {'Pb': Decimal('1.0'), 'Te': Decimal('0.7'), 'S': Decimal('0.3')}),
            ('(PbTe)_0.7(PbS)_0.3', {'Pb': Decimal('1.0'), 'Te': Decimal('0.7'), 'S': Decimal('0.3')}),
            ('(PbTe)_{0.7}(PbS)_{0.3}', {'Pb': Decimal('1.0'), 'Te': Decimal('0.7'), 'S': Decimal('0.3')}),
            ('Bi2Te3(Pb)0.3', {'Bi': Decimal('2'), 'Te': Decimal('3'), 'Pb': Decimal('0.3')}),
            ('(TePb(PbS)0.3)0.1(BiTe)2',
             {'Te': Decimal('2.1'), 'Pb': Decimal('0.13'), 'S': Decimal('0.03'), 'Bi': Decimal('2')}),
            ('(C_{10}H_{12})_{2}(Bi_2Te_3)_2.5',
             {'C': Decimal('20'), 'H': Decimal('24'), 'Bi': Decimal('5.0'), 'Te': Decimal('7.5')}),
            ('BiTe)2Te', {}),
            ('Al2(SO4_3', {}),
            ('Chemical formula for aluminium sulfate is Al2(SO4)3',
             {'Al': Decimal('2'), 'S': Decimal('3'), 'O': Decimal('12')}),
            ('Pb0.95Na0.04Te, Bi2Te3', {'Pb': Decimal('0.95'), 'Na': Decimal('0.04'), 'Te': Decimal('1')}),
        )
        for (input, output) in input_and_output_pairs:
            self.assertEqual(parse_atomic_ratio(input), output)

    def test_get_expanded_chemical_formula(self):
        input_and_output_pairs = (
            ('(PbTe)0.7(PbS)0.3', 'Te0.7Pb0.7S0.3Pb0.3'),
            ('(TePb(PbS)0.3)0.1(BiTe)2', 'S0.03Te0.1Pb0.13Te2Bi2'),
        )
        for (input, output) in input_and_output_pairs:
            self.assertEqual(_get_expanded_chemical_formula(input), output)

        inputs_raising_error = ('BiTe)2Te', 'Al2(SO4_3', '(Ale)2')
        for input in inputs_raising_error:
            with self.assertRaises(ValueError):
                _expand_innermost_atomic_ratio(input)

    def test_expand_innermost_atomic_ratio(self):
        input_and_output_pairs = (
            ('(PbTe)0.7(PbS)0.3', 'Te0.7Pb0.7(PbS)0.3'),
            ('Bi2Te3(Pb)0.3', 'Bi2Te3Pb0.3'),
            ('(TePb(PbS)0.3)0.1(BiTe)2', '(TePbS0.3Pb0.3)0.1(BiTe)2'),
        )
        for (input, output) in input_and_output_pairs:
            self.assertEqual(_expand_innermost_atomic_ratio(input), output)

        inputs_raising_error = ('BiTe)2Te', 'Al2(SO4_3', '(Ale)2')
        for input in inputs_raising_error:
            with self.assertRaises(ValueError):
                _expand_innermost_atomic_ratio(input)

    def test_convert_ratio_dict_to_str(self):
        input_and_output_pairs = (
            ({'Bi': Decimal('2'), 'Te': Decimal('3')}, 'Te3Bi2'),
            ({'Pb': Decimal('1'), 'Te': Decimal('0.7'), 'S': Decimal('0.3')}, 'S0.3Te0.7Pb1'),
        )
        for (input, output) in input_and_output_pairs:
            self.assertEqual(_convert_ratio_dict_to_str(input), output)

    def test_parse_atomic_ratio_from_expanded_chemical_formula(self):
        formula_ratios_pairs = (
            ('Bi2Te3', {'Bi': Decimal('2'), 'Te': Decimal('3')}),
            ('PbTe0.7S0.3', {'Pb': Decimal('1'), 'Te': Decimal('0.7'), 'S': Decimal('0.3')}),
        )
        for (formula, ratios) in formula_ratios_pairs:
            self.assertEqual(_parse_atomic_ratio_from_expanded_chemical_formula(formula), ratios)

        formulas_raising_error = ('(BiTe)2Te', 'Bi_2Te_3', 'Al_2(SO_4)_3')
        for formula in formulas_raising_error:
            with self.assertRaises(ValueError):
                _parse_atomic_ratio_from_expanded_chemical_formula(formula)

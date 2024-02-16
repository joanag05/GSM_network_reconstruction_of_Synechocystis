
import re

def get_lipid_formula( lipid_class, lipid):
    """Return the molecular formula of a lipid.

    Parameters
    ----------
    lipid : str
        The lipid name.

    Returns
    -------
    str
        The molecular formula of the lipid.

    """
    if lipid_class == "sqdg":
        return sqdg_formula(lipid)
    if lipid_class == "mgdg":
        return mgdg_formula(lipid)

def sqdg_formula(lipid):
    """Return the molecular formula of SQDG.

    Returns
    -------
    str
        The molecular formula of SQDG.

    """
    head_formula = 'C9H15O8S'
    acyl_1 = lipid[4:7]
    acyl_2 =lipid[7:]

    acyl_1_formula = calculate_fatty_acid_formula(int(acyl_1[0:2]), int(acyl_1[2:]))
    acyl_2_formula = calculate_fatty_acid_formula(int(acyl_2[0:2]), int(acyl_2[2:]))


    final_formula = sum_formula([head_formula, acyl_1_formula, acyl_2_formula])
    print(final_formula)

def add_formula(formula_dict, formula):
    """Add a formula to a formula dictionary.

    Parameters
    ----------
    formula_dict : dict
        The formula dictionary.
    formula : str
        The formula to add.

    Returns
    -------
    dict
        The formula dictionary with the added formula.

    """
    if formula == '':
        return formula_dict
    
    pattern = r'([A-Z][a-z]?)(\d*)'
    elements = re.findall(pattern, formula)
    
    for element, count in elements:
        element = element.capitalize()
        count = int(count) if count else 1
        formula_dict[element] = formula_dict.get(element, 0) + count
    
    return formula_dict


def formula_dict_to_formula(formula_dict):
    """Return the formula of a formula dictionary.

    Parameters
    ----------
    formula_dict : dict
        The formula dictionary.

    Returns
    -------
    str
        The formula of the formula dictionary.

    """
    formula = ''
    for element, count in sorted(formula_dict.items()):
        formula += element
        if count > 1:
            formula += str(count)
    return formula          


def sum_formula(formulas):
    """Return the sum formula of a list of formulas.

    Parameters
    ----------
    formulas : list of str
        The list of formulas.

    Returns
    -------
    str
        The sum formula of the list of formulas.

    """
    formula_dict = {}
    for formula in formulas:
        formula_dict = add_formula(formula_dict, formula)
    return formula_dict_to_formula(formula_dict)

def calculate_fatty_acid_formula(carbon_count, double_bond_count):

    hydrogen_count = 2 * carbon_count - double_bond_count * 2 - 1

    formula = f"C{carbon_count}H{hydrogen_count}O2"

    return formula


def mgdg_formula(lipid):
    """Return the molecular formula of MGDG.

    Returns
    -------
    str
        The molecular formula of MGDG.

    """
    head_formula = 'C9H16O6'
    acyl_1 = lipid[4:7]
    acyl_2 =lipid[7:]

    acyl_1_formula = calculate_fatty_acid_formula(int(acyl_1[0:2]), int(acyl_1[2:]))
    acyl_2_formula = calculate_fatty_acid_formula(int(acyl_2[0:2]), int(acyl_2[2:]))


    final_formula = sum_formula([head_formula, acyl_1_formula, acyl_2_formula])
    print(final_formula)


if __name__ == '__main__':
    get_lipid_formula('sqdg', 'sqdg160183')
    #get_lipid_formula('mgdg', 'mgdg160183')


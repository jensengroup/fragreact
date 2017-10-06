#!/home/charnley/opt/anaconda/envs/my-rdkit-env/bin/python

import numpy as np
import re
from rdkit import Chem


def canonical(smiles):
    """
    SMILES provided is canonical, so the output should be the same no matter
    how a particular molecule is input
    """
    m = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(m)
    return smiles


def count_smiles(smiles_list):

    smiles_dict = {}

    components, components_count = np.unique(smiles_list, return_counts=True)

    for comp, count in zip(components, components_count):
        smiles_dict[comp] = count

    return smiles_dict


def substract_smiles(A, B):
    """
    A - B = Cp + Cn

    where Cp has positive results
    and Cn has negative results

    """
    Cp = []
    Cn = []
    A = count_smiles(A)
    B = count_smiles(B)
    for key in np.unique(A.keys() + B.keys()):

        if key not in A:
            Cn += [key] * B[key]
            continue

        if key not in B:
            Cp += [key] * A[key]
            continue

        diff = A[key] - B[key]

        if diff == 0:
            continue
        elif diff > 0:
            Cp += [key]*diff
        elif diff < 0:
            Cn += [key]*abs(diff)

    return Cp, Cn


def tuning(left_side, right_side):

    corrected_left = []
    corrected_right = []

    left_side = count_smiles(left_side)
    right_side = count_smiles(right_side)

    for key in np.unique(left_side.keys() + right_side.keys()):

        if key not in left_side:
            print "hello"

        if key not in right_side:
            print "hello2"

        diff = right_side[key] - left_side[key]

        if diff == 0:
            continue

        elif diff > 0:
            corrected_left += [key] * diff

        elif diff < 0:
            corrected_right += [key] * diff

    return corrected_left, corrected_right




def get_bond_type(m, a, b):

    try:
        bond_type = str(m.GetBondBetweenAtoms(a, b).GetBondType())

    except AttributeError:
        return False

    if bond_type == "SINGLE":
        bond = ""

    elif bond_type == "DOUBLE":
        bond = "="

    elif bond_type == "TRIPLE":
        bond = "#"

    else:
        bond = False

    return bond


def get_components_scheme1(smiles):

    c1 = Chem.MolFromSmarts("[*]~[*]")

    m = Chem.MolFromSmiles(smiles)
    substructures = m.GetSubstructMatches(c1)

    components = []

    for sub in substructures:

        a, b = sub

        ab = get_bond_type(m, a, b)

        a = m.GetAtomWithIdx(a).GetSymbol()
        b = m.GetAtomWithIdx(b).GetSymbol()

        component = a + ab + b
        components.append(component)

    components = [canonical(component) for component in components]

    return components


def get_components_scheme2(smiles):

    c2 = "[*]~[D2]~[*]"
    c3 = "[*]~[D3](~[*])~[*]"
    c4 = "[*]~[*](~[*])(~[*])~[*]"

    c2 = Chem.MolFromSmarts(c2)
    c3 = Chem.MolFromSmarts(c3)
    c4 = Chem.MolFromSmarts(c4) # TODO

    m = Chem.MolFromSmiles(smiles)
    substructures = m.GetSubstructMatches(c2)

    components = []

    for sub in substructures:

        a, b, c = sub

        ab = get_bond_type(m, a, b)
        ac = get_bond_type(m, a, c)
        bc = get_bond_type(m, b, c)

        a = m.GetAtomWithIdx(a).GetSymbol()
        b = m.GetAtomWithIdx(b).GetSymbol()
        c = m.GetAtomWithIdx(c).GetSymbol()

        component = a + ab + b + bc + c
        components.append(component)


    substructures = m.GetSubstructMatches(c3)

    for sub in substructures:

        a, b, c, d = sub

        ab = get_bond_type(m, a, b)
        ac = get_bond_type(m, a, c)
        ad = get_bond_type(m, a, d)
        bc = get_bond_type(m, b, c)
        bd = get_bond_type(m, b, d)
        cd = get_bond_type(m, c, d)

        a = m.GetAtomWithIdx(a).GetSymbol()
        b = m.GetAtomWithIdx(b).GetSymbol()
        c = m.GetAtomWithIdx(c).GetSymbol()
        d = m.GetAtomWithIdx(d).GetSymbol()

        component = a + ab + b + "(" + bc + c + ")" + bd + d
        components.append(component)

    components = [canonical(component) for component in components]

    return components


def decompontent_scheme2(smiles):
    """
    Tune the equation
    A (bb) => aa

    where
    A (target) is big smiles
    aa (scheme2 components) is scheme2 components
    bb (single bonds) is additional bonds required, to have equald bonds on each side

    this is done for each A which consists of len(aa) > 0

    """

    components = get_components_scheme2(smiles)

    if len(components) == 0:
        return [], []

    bonds_leftside = get_components_scheme1(smiles)
    bonds_rightside = []

    for component in components:
        bonds_rightside += get_components_scheme1(component)

    left, right = tuning(bonds_leftside, bonds_rightside)

    right += components

    return left, right


def resultant(reactants, products, scheme=1):

    # Clean format

    for i, reactant in enumerate(reactants):
        reactant = reactant.split(".")
        if len(reactant) > 1:
            reactants[i] = reactant[0]
            reactants += reactant[1:]

    for i, product in enumerate(products):
        product = product.split(".")
        if len(product) > 1:
            products[i] = product[0]
            products += product[1:]

    reactants = [canonical(reactant) for reactant in reactants]
    products = [canonical(product) for product in products]

    # TODO Add different schemes

    reactants_leftside = []
    reactants_rightside = []
    products_leftside = []
    products_rightside = []

    reactants_missing = []
    products_missing = []

    for reactant in reactants:
        left, right = decompontent_scheme2(reactant)

        if len(left) == 0 and len(right) == 0:
            reactants_missing += [reactant]

        reactants_leftside += left
        reactants_rightside += right

    for product in products:
        left, right = decompontent_scheme2(product)

        if len(left) == 0 and len(right) == 0:
            products_missing += [product]

        products_leftside += left
        products_rightside += right


    left_positive, left_negative = substract_smiles(products_leftside, reactants_leftside)
    right_positive, right_negative = substract_smiles(products_rightside, reactants_rightside)

    left = left_positive + right_negative + reactants_missing
    right = right_positive + left_negative + products_missing

    left, right = substract_smiles(left, right)

    return left, right



if __name__ == "__main__":

    import argparse

    description = """Example:
fragmentreaction -r "C1=CC=CC1" "C=C" -p "C1=C[C@H]2C[C@@H]1CC2"
fragmentreaction -f filename.csv"""

    epilog = """ """

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-s', '--scheme', type=int, help='Level of fragmentation', metavar='int')

    parser.add_argument('-r', '--reactants', nargs='+', type=str, help='Reactants of the reaction', metavar='SMILES')
    parser.add_argument('-p', '--products', nargs='+', type=str, help='Products of the reaction', metavar='SMILES')

    parser.add_argument('-f', '--filename', type=str, help='File with reaction smiles', metavar='file')

    args = parser.parse_args()

    if args.reactants and args.products:

        # TODO test
        # get C6H8 for both sides to test

        print resultant(args.reactants, args.products, scheme=2)



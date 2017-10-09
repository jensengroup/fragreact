#!/home/charnley/opt/anaconda/envs/my-rdkit-env/bin/python

from __future__ import print_function

import numpy as np
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def save_reaction(reactants, products, filename):
    rxn = ".".join(reactants) + ">>"+".".join(products)
    rxn = AllChem.ReactionFromSmarts(rxn)
    img = Draw.ReactionToImage(rxn, subImgSize=(200, 200))
    img.save(filename)
    return


def print_smiles(smiles_list):

    smiles_dict = count_smiles(smiles_list)
    keys = smiles_dict.keys()
    keys.sort()

    out = []

    for key in keys:
        out += [str(smiles_dict[key]) + " " + key]

    return " ".join(out)


def print_reaction(reactants, products, human=False):

    if not human:
        print(">>".join([".".join(left), ".".join(right)]))

    else:
        reactants = print_smiles(reactants)
        products = print_smiles(products)

        print(reactants, ">>", products)

    return


def canonical(smiles):
    """
    SMILES provided is canonical, so the output should be the same no matter
    how a particular molecule is input
    """
    m = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(m)
    return smiles


def kekulize(smiles):
    m = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(m)
    smiles = Chem.MolToSmiles(m, kekuleSmiles=True)
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
            print("hello")
            quit()

        if key not in right_side:
            print("hello2")
            quit()

        diff = right_side[key] - left_side[key]

        if diff == 0:
            continue

        elif diff > 0:
            corrected_left += [key] * diff

        elif diff < 0:
            corrected_right += [key] * diff

    return corrected_left, corrected_right


def get_bond_type(m, a, b):

    # NOTE
    # If m is not kekulized then bonds can be AROMATIC
    # which is a problem for the component schemes

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


def get_atoms(smiles, ignore_hydrogen=True):

    smiles = kekulize(smiles)

    p = re.compile(r"[A-Z][a-z]?")
    atoms = p.findall(smiles)

    if ignore_hydrogen:
        atoms = [atom for atom in atoms if atom != "H"]

    return atoms


def get_components_scheme1(smiles):

    c1 = Chem.MolFromSmarts("[*]~[*]")

    m = Chem.MolFromSmiles(smiles)

    if kekulize:
        Chem.Kekulize(m)

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


def get_components_scheme2(smiles, kekulize=True):

    c2 = "[*]~[D2]~[*]"
    c3 = "[*]~[D3](~[*])~[*]"
    c4 = "[*]~[*](~[*])(~[*])~[*]"

    c2 = Chem.MolFromSmarts(c2)
    c3 = Chem.MolFromSmarts(c3)
    c4 = Chem.MolFromSmarts(c4)

    m = Chem.MolFromSmiles(smiles)

    if kekulize:
        Chem.Kekulize(m)

    substructures = m.GetSubstructMatches(c2)

    components = []

    for sub in substructures:

        a, b, c = sub

        ab = get_bond_type(m, a, b)
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
        bc = get_bond_type(m, b, c)
        bd = get_bond_type(m, b, d)

        a = m.GetAtomWithIdx(a).GetSymbol()
        b = m.GetAtomWithIdx(b).GetSymbol()
        c = m.GetAtomWithIdx(c).GetSymbol()
        d = m.GetAtomWithIdx(d).GetSymbol()

        component = a + ab + b + "(" + bc + c + ")" + bd + d
        components.append(component)

    substructures = m.GetSubstructMatches(c4)

    for sub in substructures:

        a, b, c, d, e = sub

        ab = get_bond_type(m, a, b)
        bc = get_bond_type(m, b, c)
        bd = get_bond_type(m, b, d)
        be = get_bond_type(m, b, e)

        a = m.GetAtomWithIdx(a).GetSymbol()
        b = m.GetAtomWithIdx(b).GetSymbol()
        c = m.GetAtomWithIdx(c).GetSymbol()
        d = m.GetAtomWithIdx(d).GetSymbol()
        e = m.GetAtomWithIdx(e).GetSymbol()

        component = a + ab + b
        component += "(" + bc + c + ")"
        component += "(" + bd + d + ")"
        component += be + e

        components.append(component)

    components = [canonical(component) for component in components]

    return components


def decompontent_scheme1(smiles):
    """
    Tune the equation
    A (bb) => aa

    where
    A (target) is big smiles
    aa (scheme1 components) is scheme2 components
    bb (atoms) is additional bonds required, to have equald bonds on each side

    this is done for each A which consists of len(aa) > 0

    """

    components = get_components_scheme1(smiles)

    if len(components) == 0:
        return [], []

    bonds_leftside = get_atoms(smiles)
    bonds_rightside = []

    for component in components:
        bonds_rightside += get_atoms(component)

    left, right = tuning(bonds_leftside, bonds_rightside)

    right += components

    return left, right


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
    """
    assummed that smiles lists are both split(".") and canonical at this point

    """

    reactants_leftside = []
    reactants_rightside = []
    products_leftside = []
    products_rightside = []

    reactants_missing = []
    products_missing = []

    if scheme == 1:
        decompontent_scheme = decompontent_scheme1
    elif scheme == 2:
        decompontent_scheme = decompontent_scheme2

    for reactant in reactants:
        left, right = decompontent_scheme(reactant)

        if len(left) == 0 and len(right) == 0:
            reactants_missing += [reactant]

        reactants_leftside += left
        reactants_rightside += right

    for product in products:
        left, right = decompontent_scheme(product)

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


def split_smiles(smiles_list):

    for i, smiles in enumerate(smiles_list):
        smiles = smiles.split(".")
        if len(smiles) > 1:
            smiles_list[i] = smiles[0]
            smiles_list += smiles[1:]

    return smiles_list


def fragreact(reactants, products, scheme, save_image=False, name="out"):

    reactants = [canonical(smiles) for smiles in reactants]
    products = [canonical(smiles) for smiles in products]

    ratoms = []
    for smiles in reactants:
        ratoms += get_atoms(smiles)

    patoms = []
    for smiles in products:
        patoms += get_atoms(smiles)

    ratoms.sort()
    patoms.sort()

    if ratoms != patoms:
        print(name, "atom mismatch:", count_smiles(ratoms), ">>", count_smiles(patoms))
        quit()

    if save_image:
        save_reaction(reactants, products, "reaction-"+name+".png")

    left, right = resultant(reactants, products, scheme=scheme)

    if save_image:
        save_reaction(left, right, "reaction-out-s"+str(scheme)+".png")

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

    parser.add_argument('-s', '--scheme', type=int, help='Level of fragmentation', metavar='int', default=1)
    parser.add_argument('-i', '--image', action='store_true', help='Save image of reaction')
    parser.add_argument('-u', '--human', action='store_true', help='Human readable output')

    parser.add_argument('-r', '--reactants', nargs='+', type=str, help='Reactants of the reaction', metavar='SMILES')
    parser.add_argument('-p', '--products', nargs='+', type=str, help='Products of the reaction', metavar='SMILES')

    parser.add_argument('-f', '--filename', type=str, help='File with reaction smiles', metavar='file')

    parser.add_argument('-n', '--reaction', type=str, help='Reactants and products in SMILES reaction format', metavar='SMILES')

    args = parser.parse_args()

    if args.reaction:
        reaction = args.reaction.split(">>")
        reactants = reaction[0].split(".")
        products = reaction[1].split(".")

        left, right = fragreact(reactants, products, args.scheme, save_image=args.image)
        print_reaction(left, right, human=args.human)

    if args.reactants and args.products:
        reactants = split_smiles(args.reactants)
        products = split_smiles(args.products)

        left, right = fragreact(reactants, products, args.scheme, save_image=args.image)
        print_reaction(left, right, human=args.human)

    if args.filename:
        with open(args.filename, 'r') as f:
            for line in f:
                line = line.split()
                if len(line) == 0: continue
                name = line[0]
                if ">>" in line[1]:
                    reaction = line[1].split(">>")
                    reactants = reaction[0].split(".")
                    products = reaction[1].split(".")
                else:
                    reactants = line[1].split(".")
                    products = line[2].split(".")

                left, right = fragreact(reactants, products, args.scheme, save_image=args.image, name=name)
                print(name, end=" ")
                print_reaction(left, right, human=args.human)


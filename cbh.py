#!/usr/bin/env python

import numpy as np
import re
from rdkit import Chem

def print_smiles(smiles_list, human=False):

    smiles_dict = count_smiles(smiles_list)
    keys = smiles_dict.keys()
    keys.sort()

    out = []

    for key in keys:
        out += [str(smiles_dict[key]) + " " + key]

    return " ".join(out)


def print_reaction(reactants, products, human=False):

    if not human:
        reaction = ">>".join([".".join(reactants), ".".join(products)])

    else:
        reactants = print_smiles(reactants)
        products = print_smiles(products)
        reaction =  reactants+ ">>"+ products

    return reaction


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
    """
    Count SMILES by creating a dictionary with SMILES as keys, point to the
    number of that particular SMILES.

    e.i. dict[smiles] = # of smiles

    """

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


def get_components(smiles, smart, kekulize=True):

    m = Chem.MolFromSmiles(smiles)
    smart = Chem.MolFromSmarts(smart)

    if kekulize:
        Chem.Kekulize(m)

    substructures = m.GetSubstructMatches(smart)

    components = []

    for sub in substructures:

        component = Chem.MolFragmentToSmiles(m,
            atomsToUse=sub,
            isomericSmiles=True,
            kekuleSmiles=True,
            canonical=True)

        mc = Chem.MolFromSmiles(component)
        n_atoms = mc.GetNumAtoms()
        n_bonds = len(mc.GetBonds())

        component = Chem.MolToSmiles(mc)

        if "+" in component or "-" in component:

            # Very awful hack to fix the charged molecules and their explicit
            # hydrogens

            charges = np.zeros(n_atoms, dtype=int)

            for idx in range(n_atoms):
                atom = mc.GetAtomWithIdx(idx)
                atom.SetNumExplicitHs(0)
                charge = atom.GetFormalCharge()
                charges[idx] = charge
                atom.SetFormalCharge(0)

            component = Chem.MolToSmiles(mc, canonical=False)
            component = component.replace("[", "").replace("]","")

            mc = Chem.MolFromSmiles(component)

            for idx, charge in zip(range(n_atoms), charges):
                atom = mc.GetAtomWithIdx(idx)
                atom.SetFormalCharge(charge)

            component = Chem.MolToSmiles(mc)


        if n_atoms <= n_bonds:

            mw = Chem.RWMol(m)

            if len(sub) == 3:
                mw.RemoveBond(sub[0], sub[-1])

            elif len(sub) == 4 or len(sub) == 5:
                for i in range(0, n_atoms):
                    for j in range(i+1, n_atoms):
                        if i == 1 or j == 1: continue
                        mw.RemoveBond(sub[i], sub[j])

            component = Chem.MolFragmentToSmiles(mw,
                    atomsToUse=sub,
                    isomericSmiles=True,
                    kekuleSmiles=True,
                    canonical=True)

            # print(sub, Chem.MolToSmiles(mc), component)

            if "1" in component:
                quit("Error connectivity")

        else:
            component = Chem.MolToSmiles(mc)

        # charge = Chem.GetFormalCharge(mc)
        #
        # if not charge == 0:
        #     # NOTE
        #     # Lots of lots of if case down this road
        #
        #     n_atoms = mc.GetNumAtoms()
        #
        #     for i in range(n_atoms):
        #
        #         atom = mc.GetAtomWithIdx(i)
        #         charge = atom.GetFormalCharge()
        #
        #         if not charge == 0:
        #             atom.SetFormalCharge(0)

        component = canonical(component)
        components += [component]

    return components


def get_components_scheme1(smiles, kekulize=True):

    c1 = "[*]~[*]"

    if "+" in smiles or "-" in smiles:
        pass
    else:
        return get_components(smiles, c1)

    # The code below doesn't get charges
    return get_components(smiles, c1)

    c1 = Chem.MolFromSmarts(c1)
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

    # if "+" in smiles or "-" in smiles:
    #     pass
    # else:
    components = []
    components += get_components(smiles, c2)
    components += get_components(smiles, c3)
    components += get_components(smiles, c4)
    return components

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


def decompontent(smiles, scheme=1):

    if scheme == 1: decompontent_scheme = decompontent_scheme1
    elif scheme == 2: decompontent_scheme = decompontent_scheme2

    left, right = decompontent_scheme(smiles)

    return left, right


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


def split_smiles(smiles, num_sep=None):
    """
    number seperator num_sep (e.g. 3xCC, num_spe="x")
    """

    if type(smiles) == type(""):
        smiles_list = smiles.split(".")

    else:
        smiles_list = smiles
        for i, smiles in enumerate(smiles_list):
            smiles = smiles.split(".")
            if len(smiles) > 1:
                smiles_list[i] = smiles[0]
                smiles_list += smiles[1:]

    if num_sep:
        for i, smiles in enumerate(smiles_list):
            if num_sep in smiles:
                num, smiles = smiles.split(num_sep)
                num = int(num)

                smiles_list[i] = smiles
                smiles_list += [smiles]*(num-1)

    return smiles_list


def cbh_n(reactants, products, scheme, do_canonical=True):
    """
    Use connectivity-based hieracy for reaction (reactants -> products)

    in:
        reactants -- list of SMILES
        products -- list of SMILES
        scheme -- int level of connecitivty

    out:
        left -- list of smiles for the reactant part of the CBHn reaction
        right -- list of smiles for the product part of the CBHn reaction
    """

    if do_canonical:
        reactants = [canonical(smiles) for smiles in reactants]
        products = [canonical(smiles) for smiles in products]

    left, right = resultant(reactants, products, scheme=scheme)

    return left, right


def check_reaction(reactants, products):
    """
    Check the validity of the reaction.

    Reaction should have eq. no. of atoms for both reactants and products.

    """

    ratoms = [get_atoms(smiles) for smiles in reactants]
    patoms = [get_atoms(smiles) for smiles in products]

    ratoms.sort()
    patoms.sort()

    return ratoms == patoms


if __name__ == "__main__":

    import argparse


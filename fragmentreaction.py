#!/home/charnley/opt/anaconda/envs/my-rdkit-env/bin/python

import numpy as np
import re

# rdkit
from rdkit import Chem


def fragment(smiles, scheme=1):

    m = Chem.MolFromSmiles(smiles)
    single_bond = Chem.MolFromSmarts("[*]~[*]")
    substructures = m.GetSubstructMatches(single_bond)

    components = []

    for sub in substructures:
        bond_type = m.GetBondBetweenAtoms(sub[0], sub[1]).GetBondType()
        symbols = [m.GetAtomWithIdx(sub[1]).GetSymbol(), m.GetAtomWithIdx(sub[1]).GetSymbol()]
        bond_type = str(bond_type)

        if bond_type == "SINGLE":
            connector = ""

        elif bond_type == "DOUBLE":
            connector = "="

        elif bond_type == "TRIPLE":
            connector = "#"

        else:
            print "ERROR BOND TYPE"
            quit()

        component = connector.join(symbols)
        components.append(component)

    return components


def deconstruct(smiless, ignore_hydrogen=True):

    if type(smiless) == type("x"):
        smiless = [smiless]

    # split multi smiles strings
    local_smiles = []
    for smiles in smiless:
        local_smiles += smiles.split("\.")

    p = re.compile(r"[A-Z][a-z]?")

    atoms = []
    components = []

    for smiles in local_smiles:
        atoms += p.findall(smiles)
        components += fragment(smiles)

    if ignore_hydrogen:
        atoms = [atom for atom in atoms if atom != "H"]

    atoms, atoms_count = np.unique(atoms, return_counts=True)
    components, components_count = np.unique(components, return_counts=True)

    compcount = {}
    for comp, count in zip(components, components_count):
        compcount[comp] = count

    return components, compcount, atoms, atoms_count


def get_atoms(smiles, ignore_hydrogen=True, group=False):

    p = re.compile(r"[A-Z][a-z]?")
    atoms = p.findall(smiles)

    if ignore_hydrogen:
        atoms = [atom for atom in atoms if atom != "H"]

    return atoms


def reaction(reactants, products):

    re_components, re_components_dict, re_atoms, re_atoms_count = deconstruct(reactants)
    pr_components, pr_components_dict, pr_atoms, pr_atoms_count = deconstruct(products)

    reactants = {}
    products = {}

    reactants_atoms = []
    products_atoms = []

    for key in np.unique(re_components_dict.keys() + pr_components_dict.keys()):

        atoms = get_atoms(key)

        if key not in pr_components_dict:
            reactants[key] = re_components_dict[key]
            reactants_atoms += re_components_dict[key] * atoms
            continue

        if key not in re_components_dict:
            products[key] = pr_components_dict[key]
            products_atoms += pr_components_dict[key] * atoms
            continue

        sign = re_components_dict[key] - pr_components_dict[key]

        if sign == 0: continue

        if sign > 0:
            reactants[key] = sign
            reactants_atoms += sign * atoms
        else:
            products[key] = abs(sign)
            products_atoms += abs(sign) * atoms

    for atom in np.unique(reactants_atoms + products_atoms):

        n_reactants = reactants_atoms.count(atom)
        n_products = products_atoms.count(atom)

        diff = n_reactants - n_products

        if diff == 0:
            continue

        elif diff < 0:
            reactants[atom] = abs(diff)
        else:
            products[atom] = abs(diff)

    return reactants, products


if __name__ == "__main__":

    import argparse

    description = """

For example
fragmentreaction -r "C1=CC=CC1" "C=C" -p "C1=C[C@H]2C[C@@H]1CC2"

    """

    epilog = """
    """

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-s', '--scheme', type=int, help='Level of fragmentation', metavar='int')

    parser.add_argument('-r', '--reactants', nargs='+', type=str, help='Reactants of the reaction', metavar='rect')
    parser.add_argument('-p', '--products', nargs='+', type=str, help='Products of the reaction', metavar='prod')

    parser.add_argument('-f', '--filename', type=str, help='File with reaction smiles', metavar='file')

    args = parser.parse_args()

    # if args.reactants == None or args.products == None:
    #     print "No enough arguments"
    #     print
    #     print parser.print_help()
    #     quit()


    if args.reactants and args.products:
        print reaction(args.reactants, args.products)

    if args.filename:
        with open(args.filename, 'r') as f:
            for line in f:
                line = line.split()
                if len(line) == 0: continue

                name = line[0]
                reactants = line[1]
                products = line[2]

                cr, cp = reaction(reactants, products)
                out = ""
                out += name
                for r in cr:
                    out += " " + str(cr[r]) + "*" + r

                out += " ->"
                for p in cp:
                    out += " " + str(cp[p]) + "*" + p

                print out






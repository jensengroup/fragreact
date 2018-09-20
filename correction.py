
# from rdkit import Chem
# from rdkit.Chem import AllChem

import numpy as np


def make_it_float_or_not(val):

    try:
        val = float(val)
    except ValueError:
        pass

    return val

def get_database(filename):

    f = open(filename, 'r')

    header = next(f)
    header = header.split(", ")
    header = [head.strip() for head in header]

    idx_smiles = header.index("smiles")

    data = {}
    for line in f:
        line = line.split(",")
        line = [col.strip() for col in line]
        line = [make_it_float_or_not(val) for val in line]

        smiles = line[0]
        data[smiles] = line

    return header, data


def get_energy_smiles(smi, idx, DB, debug=False):

    R = 1.9872036 * 10**-3 # kcal K−1 mol−1
    T = 298.0

    if smi == "[H+]":
        energy = 2.5 * R * T
    else:
        energy = DB[smi][idx]

        if debug:
            if np.isnan(energy):
                print("  ", smi, energy)

    return energy


def get_energy_rxn(rxnsmi, idx_reference, idx_methods, DB, debug=False):

    if rxnsmi == ">>": return 0.0

    # reaction
    rxnsmi = rxnsmi.split(">>")
    rxnsmi = [rxn.split(".") for rxn in rxnsmi]
    sign = [-1, +1]

    # start correcting
    energy_ref = 0.0
    for side, one in zip(rxnsmi, sign):
        for smi in side:
            energy_ref += one*get_energy_smiles(smi, idx_reference, DB, debug=debug)

    energy_met = 0.0
    for idx in idx_methods:
        for side, one in zip(rxnsmi, sign):
            for smi in side:
                energy_met += one*get_energy_smiles(smi, idx, DB, debug=debug)

    return energy_ref - energy_met


if __name__ == "__main__":

    """
    enthalpi(H+) = 2.5RT

    """

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--database', type=str, help='SMILES and energy database', metavar='file')

    parser.add_argument('-s', '--scheme', type=str, help='correction', metavar='str')
    parser.add_argument('-f', '--filename', type=str, help='cbh file', metavar='str')

    parser.add_argument('-r', '--reference', type=str, help='Reference to use', default='G4')
    parser.add_argument('-m', '--methods', nargs='+', type=str, help='Columns to use from database')

    args = parser.parse_args()

    DEBUG = True

    header, DB = get_database(args.database)

    idx_reference = header.index(args.reference)
    idx_methods = []
    for x in args.methods:
        idx_methods.append(header.index(x))

    if args.filename:
        with open(args.filename, 'r') as f:
            for line in f:
                name, rxnsmi = line.split()
                print(name, get_energy_rxn(rxnsmi, idx_reference, idx_methods, DB, debug=DEBUG))

    if args.scheme:
        print(get_energy(args.scheme, idx_reference, idx_methods, DB))



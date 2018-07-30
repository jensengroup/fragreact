
import numpy as np
import hashlib
import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

import shell as sh
import re

def read_components(filename):
    smiles_list = []
    f = open(filename, 'r')
    for line in f:
        smiles_list.append(line.strip())
    return smiles_list


def get_mopac_energy(filename):

    regex = r'[\-]*\d+\.\d+[eE\-]*\d*'

    try:
        energy = sh.shell('grep "HEAT OF FORMATION" '+filename, shell=True)
        energy = re.findall(regex, energy)
        energy = energy[0] # kcal/mol
        energy = float(energy)
    except:
        energy = float("nan")

    return energy


def smiles2hash(smiles):
    molstr = hashlib.md5(smiles.encode('utf-8')).hexdigest()
    return molstr


def get_conformations(smiles, max_conf=20):

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(m)

    confs = min(1 + 3*rot_bond, max_conf)

    AllChem.EmbedMultipleConfs(m, numConfs=confs,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True)

    conf_list = []

    for i, conf in enumerate(m.GetConformers()):
        tm = Chem.Mol(m, False, conf.GetId())
        conf_list.append(tm)

    return conf_list

def generate_conformations(args):

    smiles_list = read_components(args.filename)

    for j, smiles in enumerate(smiles_list):
        m_list = get_conformations(smiles)

        molname = args.prefix

        if args.md5:
            molname += smiles2hash(smiles)
        else:
            molname += str(j)

        if os.path.isfile(args.folder + molname + "_0.sdf"):
            print molname, smiles
            continue

        for k, m in enumerate(m_list):

            name = molname + "_" + str(k)
            print name, smiles

            writer = Chem.SDWriter(args.folder + name + ".sdf")
            writer.write(m)

    return


def find_all_in_folder(folder, srch):

    print srch

    files = [f for f in os.listdir(folder) if re.match(srch, f)]
    files = [f for f in files if "out" in f]
    files.sort()

    return files


def already_found_low(folder):

    files = os.listdir(folder)
    files = [f.replace(".sdf", "") for f in files]

    return files


def find_lowest(args):

    f = open(args.filename, 'r')

    # already in low folder
    lowfldr = already_found_low(args.folder + "low")


    data = {}

    for line in f:
        line = line.split()
        key = line[0]
        smiles = line[1]

        if key in lowfldr: continue

        # find all out files
        files = find_all_in_folder(args.folder, key)
        files = [f.replace(".out", "") for f in files]

        for molname in files:

            energy = get_mopac_energy(args.folder + molname + ".out")
            key, idx = molname.split("_")

            if key not in data:
                data[key] = {}
                data[key]['energy'] = energy
                data[key]['idx'] = idx
                data[key]['smiles'] = smiles

            else:
                if energy > data[key]['energy']:
                    data[key]['energy'] = energy
                    data[key]['idx'] = idx
                else:
                    continue

    keys = data.keys()
    keys.sort()

    for key in keys:
        # cp sdf to new name
        idx = data[key]['idx']
        energy = data[key]['energy']
        smiles = data[key]['smiles']

        cmd = "cp {:} {:}".format(args.folder + key + "_" + idx + ".sdf", args.folder + "low/" + key + ".sdf")
        sh.shell(cmd)

        if np.isnan(energy):
            print key, smiles, "nan"
        else:
            print key, smiles

    return


def main():

    import argparse
    import sys

    description = """ """

    epilog = """ """

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-f', '--filename', type=str, help='CSV file of either components or conformation list', metavar='file')
    parser.add_argument('-o', '--folder', type=str, default="", help='Output folder', metavar='str')
    parser.add_argument('-p', '--prefix', type=str, default="", help='String prefix on the conformations', metavar='str')
    parser.add_argument('-m', '--md5', action="store_true", help='Use MD5 of the SMILES string, instead of id')

    parser.add_argument('-s', '--find_lowest', action="store_true", help='find lowest conformation from list (assumes MOPAC output)')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if args.folder != "":
        if args.folder[-1] != "/":
            args.folder += "/"

    #

    if args.find_lowest:
        find_lowest(args)

    else:
        generate_conformations(args)

    return


if __name__ == "__main__":
    main()



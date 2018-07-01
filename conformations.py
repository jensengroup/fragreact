
import numpy as np
import hashlib
import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem


def read_components(filename):
    smiles_list = []
    f = open(filename, 'r')
    for line in f:
        smiles_list.append(line.strip())
    return smiles_list


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

    parser.add_argument('-f', '--filename', type=str, help='Components list', metavar='file')
    parser.add_argument('-o', '--folder', type=str, default="", help='Output folder', metavar='str')
    parser.add_argument('-p', '--prefix', type=str, default="", help='String prefix on the conformations', metavar='str')
    parser.add_argument('-m', '--md5', action="store_true", help='Use MD5 of the SMILES string, instead of id')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    foldername = args.folder

    if foldername != "":
        if foldername[-1] != "/":
            foldername += "/"


    smiles_list = read_components(args.filename)

    for j, smiles in enumerate(smiles_list):
        m_list = get_conformations(smiles)

        molname = args.prefix

        if args.md5:
            molname += smiles2hash(smiles)
        else:
            molname += str(j)

        if os.path.isfile(foldername + molname + "_0.sdf"):
            print "already exists:", molname, smiles
            quit()

        for k, m in enumerate(m_list):

            name = molname + "_" + str(k)
            print name, smiles

            writer = Chem.SDWriter(foldername + name + ".sdf")
            writer.write(m)

    return


if __name__ == "__main__":
    main()



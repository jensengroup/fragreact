
import numpy as np
import re
from rdkit import Chem

if __name__ == "__main__":

    import sys
    args = sys.argv[1:]

    smiless = args

    for smiles in smiless:
        m = Chem.MolFromSmiles(smiles)



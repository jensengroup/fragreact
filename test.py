import numpy as np
import cbh

def test_fragmentation_cbh1():
    test_smiles = "C1=C[C@H]2C[C@@H]1CC2"
    assert sorted(cbh.get_components_scheme1(test_smiles)) == sorted("CC C=C CC CC CC CC CC CC".split())
    return

def test_fragmentation_cbh2():
    test_smiles = "C1=C[C@H]2C[C@@H]1CC2"
    assert sorted(cbh.get_components_scheme2(test_smiles)) == sorted("C=CC C=CC CCC CCC CCC CC(C)C CC(C)C".split())
    return

def test_fragmentation_reaction_cbh1():

    reactants, products = ["C1=CC=CC1", "C=C"], ["C1=C[C@H]2C[C@@H]1CC2"]
    left, right = cbh.cbh_n(reactants, products, 1)

    assert sorted(left) == sorted(['C', 'C', 'C', 'C', 'C=C', 'C=C'])
    assert sorted(right) == sorted(['CC', 'CC', 'CC', 'CC'])

    return

def test_fragmentation_reaction_cbh2():

    reactants, products = ["C1=CC=CC1", "C=C"], ["C1=C[C@H]2C[C@@H]1CC2"]
    left, right = cbh.cbh_n(reactants, products, 2)

    assert sorted(left) == sorted(['C=CC', 'C=CC', 'CC', 'CC', 'CC', 'CC'])
    assert sorted(right) == sorted(['CC(C)C', 'CC(C)C', 'CCC', 'CCC'])

    return

def test_split_smiles():

    assert cbh.split_smiles("CC.CC") == ["CC", "CC"]
    assert cbh.split_smiles("2;CC", num_sep=";") == ["CC", "CC"]
    assert cbh.split_smiles(["CC", "CC.CC"]) == ["CC", "CC", "CC"]

    return


def test_get_components_scheme1():

    smiles = "C=[NH+]C"
    components = ["C=[NH2+]", "C[NH3+]"]
    output = cbh.get_components_scheme1("C=[NH+]C")
    assert sorted(components) == sorted(output)

    smiles = "C#[N+]C"
    components = ["C#[NH+]", "C[NH3+]"]
    output = cbh.get_components_scheme1(smiles)
    assert sorted(components) == sorted(output)

    smiles = "CC(=O)[O-]"
    components = ["CC", "C=O", "C[O-]"]
    output = cbh.get_components_scheme1(smiles)
    assert sorted(components) == sorted(output)

    smiles = "C[S+](C)C"
    components = ['C[SH2+]', 'C[SH2+]', 'C[SH2+]']
    output = cbh.get_components_scheme1(smiles)
    assert sorted(components) == sorted(output)

    return

def test_get_components_scheme2():

    fun = cbh.get_components_scheme2

    # getting the right number of H on N
    smiles = "CCc1c[nH]c2ccccc12"
    components = ['CCC', 'C=CN', 'CNC', 'C=CC', 'C=CC', 'C=CC', 'C=CC', 'C=C(C)C', 'C=C(C)C', 'C=C(C)N']
    output = fun(smiles)
    assert sorted(components) == sorted(output)

    # connected smiles
    smiles = "C1CO1"
    components = ['CCO', 'CCO', 'COC']
    output = fun(smiles)
    assert sorted(components) == sorted(output)

    return


if __name__ == "__main__":
    print("use python3 -m pytest test.py")


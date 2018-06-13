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

if __name__ == "__main__":
    test_split_smiles()
    test_fragmentation_cbh1()
    test_fragmentation_cbh2()
    test_fragmentation_reaction_cbh1()
    test_fragmentation_reaction_cbh2()


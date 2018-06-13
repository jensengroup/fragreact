# fragreact

Molecular reaction fragmentation scheme towards improving the accuracy of enthalpy calculation.
Following the procedure of 10.1021/acs.orglett.7b00891 and dx.doi.org/10.1021/ct200279q.
For example the reaction

    C1=CC=CC1 C=C  >> C1=C[C@H]2C[C@@H]1CC2

should be called as

    ./fragmentreaction.py -r "C1=CC=CC1" "C=C" -p "C1=C[C@H]2C[C@@H]1CC2"

and will give, for scheme 1

    4 C + 2 C=C >> 4 CC

and for scheme 2

    2 C=C + 4 CC >> 2 CC(C)C + 2 CCC


## Example flow using a list of reactions

using a file with list of reactions in SMILES format, for example `example/parent_reactions.csv`, we can create the corresponding correction reaction.

    anaconda fragreact.py --filename ./example/parent_reactions.csv --scheme 1 > ./example/corr_1.csv

If we want all the SMILES needed to create a database for these reactions we
can first run with scheme 1 and then scheme 2 (piping it down in the same file)

    anaconda fragreact.py --filename ./example/parent_reactions.csv --scheme 1 --save_database > ./example/components.csv
    anaconda fragreact.py --filename ./example/parent_reactions.csv --scheme 2 --save_database >> ./example/components.csv

we can then remove dubplicates using `sort` 

    sort -u ./examples/components.csv -o ./examples/components.csv

Using this list of SMILES you generate conformations and convert it from XYZ to any QM software.

    mkdir xyz
    anaconda conformations.py --filename ./examples/components.csv --prefix xyz/g_

Then use the SDF files to do the needed calculations at the chosen QM level.
Remmeber to find the conormation for each SMILES with lowest energy and use it for all levels.
From the SMILES/Calculations you create a CSV file in the format of

    SMILES, Reference, Method 1, Method 2, Method 3

In the following example the method G4 is used as reference to correct a list of other QM methods.

    SMILES,            G4,      B3LYP/6-311G(dp),        HF/6-311G(dp),            HF/STO-3G,          AM1,                  PM3,                  PM6
         C,  -25389.96705,          -25405.49289,         -25199.73379,         -24892.69780,     -8.78979,            -13.02534,            -12.28825
    CC(C)C,  -99328.65850,          -99379.16674,         -98641.86883,         -97453.81644,    -29.42010,            -29.58106,            -27.50766

You can add as many method/coloumns as needed.
Using the correction CSV file created early, we can get correction energies for each reaction using any method.
This is done by using the Coloumn names of the database CSV file for reference and method.

    anaconda correct.py --database ./example/database.csv --filename ./example/corr_1.csv --reference "G4" --method "HF/6-311G(dp)"

This will output the reaction name and the correction, which are added to the method energy.

    energy_ref ~= energy_met + energy_corr


## I don't have a SMILES string for my structure

There is a app for that. To get the proper SMILES you need to do a QM calculation.




# fragreact

Molecular reaction fragmentation scheme towards improving the accuracy of enthalpy calculation.

Following the procedure of 10.1021/acs.orglett.7b00891 and dx.doi.org/10.1021/ct200279q.


for example the reaction

    C1=CC=CC1 C=C  >> C1=C[C@H]2C[C@@H]1CC2

should be called as 

    ./fragmentreaction.py -r "C1=CC=CC1" "C=C" -p "C1=C[C@H]2C[C@@H]1CC2"

and will give, for scheme 1

    4 C + 2 C=C >> 4 CC

and for scheme 2

    2 C=C + 4 CC >> 2 CC(C)C + 2 CCC





import cbh
import numpy as np



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

    parser.add_argument('-s', '--scheme', type=int, help='Level of fragmentation', metavar='int', default=1)
    parser.add_argument('-i', '--image', action='store_true', help='Save image of reaction')
    parser.add_argument('-u', '--human', action='store_true', help='Human readable output')
    parser.add_argument('-d', '--save_database', action='store_true', help='Save all components and print out list')

    parser.add_argument('-x', '--decomponent', nargs='+', type=str, help='Decompenent single one SMILES', metavar='SMILES')

    parser.add_argument('-r', '--reactants', nargs='+', type=str, help='Reactants of the reaction', metavar='SMILES')
    parser.add_argument('-p', '--products', nargs='+', type=str, help='Products of the reaction', metavar='SMILES')
    parser.add_argument('-n', '--reaction', type=str, help='Reactants and products in SMILES reaction format', metavar='SMILES')

    parser.add_argument('-f', '--filename', type=str, help='File with reaction smiles', metavar='file')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if args.save_database: components = []


    if args.reaction:
        reaction = args.reaction.split(">>")
        reactants = reaction[0].split(".")
        products = reaction[1].split(".")
        left, right = cbh.cbh_n(reactants, products, args.scheme)


        if args.save_database:
            components += left + right
        else:
            print(cbh.print_reaction(left, right, human=args.human))


    elif args.reactants:
        pass


    elif args.filename:

        with open(args.filename, 'r') as f:

            for line in f:
                line = line.split()

                if len(line) == 0: continue

                name = line[0]

                if ">>" in line[1]:
                    reaction = line[1].split(">>")
                    reactants = reaction[0].split(".")
                    products = reaction[1].split(".")
                else:
                    reactants = line[1].split(".")
                    products = line[2].split(".")

                left, right = cbh.cbh_n(reactants, products, args.scheme)

                if args.save_database:
                    components += left + right
                else:
                    print(name, cbh.print_reaction(left, right, human=args.human))

    elif args.decomponent:
        for smiles_list in args.decomponent:
            smiles_list = smiles_list.split(".")
            for smiles in smiles_list:
                print(smiles)
                left, right = cbh.decompontent(smiles, scheme=args.scheme)
                print(" ", " ".join(left))
                print(" ", " ".join(right))


    if args.save_database:
        for comp in np.unique(components):
            print(comp)



if __name__ == "__main__":
    main()


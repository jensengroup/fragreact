
example:
	anaconda fragreact.py --filename ./example/parent_reactions.csv --scheme 1 > ./example/corr_1.csv
	anaconda fragreact.py --filename ./example/parent_reactions.csv --scheme 1 --save_database > ./example/components.csv
	anaconda fragreact.py --filename ./example/parent_reactions.csv --scheme 2 --save_database >> ./example/components.csv
	sort -u ./example/components.csv -o ./example/components.csv
	mkdir -p xyz
	anaconda conformations.py --filename ./example/components.csv --prefix xyz/g_ > ./example/components_conformations.csv


test:
	anaconda -m pytest test.py

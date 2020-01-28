DNA_IGEM    = ./data/features/igem/dna.db
PROTEIN_IGEM = ./data/features/igem/protein.db
DNA_SG      = ./data/features/snapgene/dna.db
PROTEIN_SG  = ./data/features/snapgene/protein.db

.PHONY: test docs build

test:
	mkdir -p ./tests/output
	python3 -m unittest discover tests -p '*_test.py'

docs:
	cd docs && make html
	git add .
	git commit -m "update docs"

reset:
	rm -f data/features/**/dna*
	rm -f data/features/**/protein*
	rm -f data/features/igem/igem.pickle

features: reset
	python3 -c "import synbio.features.cluster as c; c.cluster()"
	python2 ./synbio/features/snapgene.py
	cat ${DNA_SG} > ./data/features/dna.db
	cat ${DNA_IGEM} >> ./data/features/dna.db
	cat ${PROTEIN_SG} > ./data/features/protein.db
	cat ${PROTEIN_IGEM} >> ./data/features/protein.db
	python3 -c "import synbio.features.seed as s; s.seed()"

minor: test
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

patch: test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

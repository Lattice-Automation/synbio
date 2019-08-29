.PHONY: test docs build

test:
	mkdir -p ./tests/output
	python3 -m unittest discover tests -p '*_test.py'

docs:
	cd docs && make html
	git add .
	git commit -m "update docs"

build:
	conda build . -c jtimmons
	conda build purge

export:
	conda env export --name synbio > environment.yml

minor: test
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing
	$(MAKE) docs
	$(MAKE) build

patch: test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing
	$(MAKE) docs
	$(MAKE) build

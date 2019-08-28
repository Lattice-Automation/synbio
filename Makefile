.PHONY: docs test build

export:
	conda env export --name synbio > environment.yml

build:
	conda build . -c jtimmons

test:
	mkdir -p ./tests/output
	python3 -m unittest discover tests -p '*_test.py'

install:
	pip install -e .

minor: test
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing
	$(MAKE) docs
	$(MAKE) install
	$(MAKE) build

patch: test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing
	$(MAKE) docs
	$(MAKE) install
	$(MAKE) build

docs:
	cd docs && make html
	git add .
	git commit -m "update docs"

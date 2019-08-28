.PHONY: docs test

env:
	conda env create

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

patch: test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing
	$(MAKE) docs
	$(MAKE) install

docs:
	cd docs && make html
	git add .
	git commit -m "update docs"

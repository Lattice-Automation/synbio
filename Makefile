.PHONY: docs test

test:
	mkdir -p ./tests/output
	python3 -m unittest discover tests -p '*_test.py'

install:
	pip install -e .

minor: test docs
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

patch: test docs
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

docs:
	cd docs && make html

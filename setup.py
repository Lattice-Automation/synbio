from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="synbio",
    version="0.4.6",
    author="JJTimmons",
    author_email="jtimmons@latticeautomation.com",
    description="Synbio design and build library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 2 - Pre-Alpha",
    ],
    url="https://github.com/Lattice-Automation/synbio",
    test_suite="tests.suite",
)

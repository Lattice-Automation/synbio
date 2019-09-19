from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = ["biopython>=1.74.0", "networkx>=2.3.0", "primer3-py"]

setup(
    name="synbio",
    version="0.6.0",
    python_requires=">=3",
    author="JJTimmons",
    author_email="jtimmons@latticeautomation.com",
    url="https://github.com/Lattice-Automation/synbio",
    description="Synbio design and build library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="synbio, lab automation, dna assembly",
    license="GNU General Public License v2 (GPLv2)",
    packages=find_packages(),
    install_requires=requirements,
    test_suite="tests.suite",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

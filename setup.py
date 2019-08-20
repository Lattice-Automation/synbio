from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="synbio",
    version="0.4.5",
    author="JJTimmons",
    author_email="jtimmons@latticeautomation.com",
    description="Synbio design and build library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    url="https://latticeautomation.com/",
    test_suite="tests.suite",
)

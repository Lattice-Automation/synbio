# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import glob
import os
import shutil
import sys

sys.path.append(os.path.abspath(os.path.join("..")))

for filename in glob.glob(os.path.join(".", "*.py")):
    if "conf" not in filename:  # don't delete this...
        os.remove(filename)

for filename in glob.glob(os.path.join("..", "examples", "*.py")):
    shutil.copy(filename, ".")  # copy example file to this doc dir


# -- Project information -----------------------------------------------------

project = "synbio"
author = "jjtimmons"

# The full version, including alpha/beta/rc tags
release = "0.6.2"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "m2r",
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# see: https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html
html_theme_options = {"display_version": True, "collapse_navigation": False}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

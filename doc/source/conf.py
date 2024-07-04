# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the
# documentation:
# 
# https://www.sphinx-doc.org/en/master/usage/configuration.html.


#######################################################################


# Import from the standard library.
import os
import sys


#######################################################################


# Get the absolute path to the package.
path = os.path.abspath("../../pdbcraft")

# Insert the path to the package into the PATH.
sys.path.insert(0, path)


#######################################################################


# https://www.sphinx-doc.org/en/master/usage/configuration.html
# #project-information


# Set the name of the project.
project = "pdbcraft"

# Set the copyright of the project.
copyright = "2024, Valentina Sora"

# Set the name(s) of the project's author(s).
author = "Valentina Sora"


#######################################################################


# https://www.sphinx-doc.org/en/master/usage/configuration.html
# #general-configuration

# Set a list of strings that are module names of extensions:
#
# - 'myst_parser' is needed to parse Markwodn files, and it needs to
#   be installed separately (it can be done with 'pip install
#   myst-parser').
#
# - 'sphinx.ext.autodoc' is needed to automatically generate
#   documentation from source code files, and it is installed together
#   with Sphinx.
#
# - 'numpydoc' is needed to parse NumPy-style docstrings, and it needs
#   to be installed separately (it can be done with
#   'pip install numpydoc').
extensions = \
   ["myst_parser",
    "sphinx.ext.autodoc",
    "numpydoc",
    "sphinx_design"]

# Set how to order the functions/classes in the modules'
# documentations.
autodoc_member_order = "bysource"

# Set whether to show all class members by default.
numpydoc_show_class_members = False

# Set whether to overwrite '.rst' files when creating autosummaries.
autosummary_generate_overwrite = False

# Set the file extensions of source files. Sphinx considers the files
# with this suffix as sources. The value can be a dictionary mapping
# file extensions to file types.
source_suffix = \
	{".rst" : "restructuredtext",
	 ".txt" : "restructuredtext",
	 ".md" : "markdown"}

# Set a list of paths that contain extra templates (or templates that
# overwrite builtin/theme-specific templates). Relative paths are
# taken as relative to the configuration directory.
#
# As these files are not meant to be built, they are automatically
# added to 'exclude_patterns'.
templates_path = ["_templates"]


#######################################################################


# Set which 'myst-parser' extensions should be enabled.
myst_enable_extensions = ["amsmath", "dollarmath"]

# Set whether to create automatic heading anchors for heading up to
# level 3.
myst_heading_anchors = 3

# Set whether to raise the warning regarding having the document
# ending with footnotes preceded by a heading.
myst_footnote_transition = False


#######################################################################


# https://www.sphinx-doc.org/en/master/usage/configuration.html
# #options-for-html-output

# Set a list of paths that contain custom static files (such as style
# sheets or script files).
#
# Relative paths are taken as relative to the configuration directory.
#
# They are copied to the output’s '_static directory' after the
# theme’s static files, so a file named 'default.css' will overwrite
# the theme’s 'default.css'.
html_static_path = ["_static"]

# Set the HTML theme to be used.
# - If 'sphinx_rtd_theme', it needs to be installed with
#   'pip install sphinx-rtd-theme'.
# - If 'pydata_sphinx_theme', it needs to be installed with
#   'pip install pydata_sphinx_theme'.
html_theme = "pydata_sphinx_theme"

# Set a dictionary of values to pass into the template engine’s context
# for all pages.
html_context = { 
   # ...
   "default_mode": "light"
}

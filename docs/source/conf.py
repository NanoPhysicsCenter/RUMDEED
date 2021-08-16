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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
from datetime import date


# -- Project information -----------------------------------------------------

# Function borrowed from Lammps conf.py.in
def get_git_commit():
    import subprocess, time
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if commit.returncode == 0:
            return commit.stdout.decode()
    except:
        pass
    return ''

project = 'Vacuum-MD'
copyright = '{}, Kristinn Torfason'.format(date.today().year)
author = 'Kristinn Torfason'
release = get_git_commit()


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
        'sphinxfortran.fortran_domain',
        'sphinxfortran.fortran_autodoc',
        'sphinxcontrib.bibtex',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Master document
master_doc = 'manual'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
#html_theme = 'alabaster'

html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'prev_next_buttons_location': 'both',
    'style_external_links': True,
}

html_show_sourcelink = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Math - KT
numfig = True
math_numfig = True
numfig_secnum_depth = 2
math_eqref_format = "Eq. {number}"

def setup(app):
    app.add_css_file('css/custom.css')

# -- Fortran code
fortran_src = ['../src/']
fortran_ext = ['F90', 'f90']

# -- Bibtext refs
bibtex_bibfiles = ['refs.bib']
bibtex_default_style = 'unsrt'
bibtex_reference_style = 'label'

# -- Latex

latex_engine = 'xelatex'
latex_elements = {
    'papersize': 'a4paper',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
}
latex_show_urls = 'footnote'
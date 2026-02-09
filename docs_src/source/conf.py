from datetime import datetime
from importlib.metadata import version as get_version

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'Community Toolkit'
copyright = f'{datetime.now().year}, Sea-Bird Scientific'
author = 'Sea-Bird Scientific'
release = get_version('seabirdscientific')
version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx_mdinclude',
]

templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_title = f'{project} {version}'
html_permalinks_icon = '<span>#</span>'
html_theme = 'sbsfuro'
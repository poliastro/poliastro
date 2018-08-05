import os
from os.path import join, dirname, abspath
from datetime import datetime

import alabaster


# Alabaster theme + mini-extension
html_theme_path = [alabaster.get_path()]
extensions = [
    'alabaster',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
]
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'astropy': ('http://docs.astropy.org/en/stable/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
    'matplotlib': ('http://matplotlib.org', None)
}
#Nbsphinx configuration
if os.environ.get('READTHEDOCS') == 'True':
    nbsphinx_execute = 'never'
else:
    nbsphinx_execute = 'always'

    # Controls when a cell will time out (defaults to 30; use -1 for no timeout):
    nbsphinx_timeout = 60


# Paths relative to invoking conf.py - not this shared file
html_static_path = [join("..", "_static")]
html_theme = "alabaster"
html_favicon = 'favicon.ico'
html_theme_options = {
    'logo': 'logo_trans.png',
    'logo_name': True,
    'logo_text_align': 'center',
    'travis_button' : True,
    'codecov_button': True,
    'description':'Astrodynamics in Python',
    'body_text_align': 'left',
    'github_user': 'poliastro',
    'github_repo': 'poliastro',
    'show_relbars': True,
    'show_powered_by': False,
    'page_width': '80%',
    'github_banner': True,
    'extra_nav_links' : { 'Benchmarks': 'https://blog.poliastro.space/poliastro-benchmarks/',
                          'Blog': 'https://blog.poliastro.space/',
                        },
}
html_sidebars = {
    "**": ["about.html", "navigation.html", "searchbox.html", "donate.html"]
}
# Regular settings
pygments_style = 'sphinx'
project = "poliastro"
year = datetime.now().year
copyright = "2013 - %d,  Juan Luis Cano Rodríguez and the poliastro development team" % year
master_doc = "index"
templates_path = ["_templates"]
exclude_trees = ["_build"]
exclude_patterns = ['_build', '**.ipynb_checkpoints']
source_suffix = ".rst"
default_role = "obj"
version = '0.11'
release = '0.11.dev0'
autodoc_member_order = 'bysource'
suppress_warnings = ['image.nonlocal_uri']


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  ('index', 'poliastro.tex', 'poliastro Documentation',
   'Juan Luis Cano Rodríguez', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'poliastro', 'poliastro Documentation',
     ['Juan Luis Cano Rodríguez'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'poliastro', 'poliastro Documentation',
   'Juan Luis Cano Rodríguez', 'poliastro', 'One line description of project.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#texinfo_no_detailmenu = False

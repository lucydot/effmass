import os
from effmass import __version__

try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

try:
	import pypandoc
	long_description = pypandoc.convert('README.md', 'rst')
except ImportError:
	long_description = open('README.md').read()

config = {
    'description': 'An effective mass package',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'author': 'Lucy Whalley',
    'author_email': 'l.whalley@northumbria.ac.uk',
    'url': 'https://github.com/lucydot/effmass',
    'download_url': "https://github.com/lucydot/effmass/archive/%s.tar.gz" % (__version__),
    'version': __version__,
    'install_requires': [ 'vasppy>=0.5.0.0', 
			  'scipy',
			  'numpy', 
                          'matplotlib',
                          'adjustText',
                          'ase>=3.21.1',
                          'questionary>=1.9.0',
                          'prettytable>=2.1.0'],
    'extras_require': {
	    "docs": [
		    "sphinx >=3.2.1",
		    "sphinx_rtd_theme>=0.5.0",
	    ],
	    "tests": [
		    "pytest",
		    "pytest-lazy-fixture",
		    "code-climate-test-reporter",
		    "coverage==4.3.4"
	    ],
	    "dev": ["black"],
    },
    'python_requires': '>=3.6',
    'license': 'MIT',
    'packages': [ 'effmass' ],
    'scripts': [],
    'name': 'effmass',
    'entry_points': {"console_scripts": ["effmass = effmass.cli:cli"]}
}

setup(**config)

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
    'author_email': 'lucywhalley@gmail.com',
    'url': 'https://github.com/lucydot/effmass',
    'download_url': "https://github.com/lucydot/effmass/archive/%s.tar.gz" % (__version__),
    'version': __version__,
    'install_requires': [ 'vasppy>=0.5.0.0', 
			  'scipy',
			  'numpy', 
                          'matplotlib',
                          'adjustText',
                          'codeclimate-test-reporter',
                          'pytest-lazy-fixture',
                          'coverage==4.3.4' ],
    'python_requires': '>=3.6',
    'license': 'MIT',
    'packages': [ 'effmass' ],
    'scripts': [],
    'name': 'effmass'
}

setup(**config)

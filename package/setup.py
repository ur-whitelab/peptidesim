import os
from glob import glob
from setuptools import setup

exec(open('peptidesim/version.py').read())

setup(name = 'peptidesim',
      version = __version__,
      scripts = glob(os.path.join('scripts', '*')),
      description = 'An automated peptide simulator',
      author = 'Dilnoza Amirkulova, Andrew White',
      author_email = 'andrew.white@rochester.edu',
      url = 'http://thewhitelab.org/Software',
      license = 'GPL3',
      packages = ['peptidesim'],
      install_requires=['PeptideBuilder', 'GromacsWrapper', 'dill', 'traitlets', 'requests', 'future', 'pytest', 'biopython==1.72','pandas'],
      dependency_links=['https://github.com/mtien/PeptideBuilder/archive/master.zip#egg=PeptideBuilder-1.0.4'],
      test_suite='tests',
      zip_safe=True,
      package_data = {'peptidesim': ['templates/*.mdp']}
)

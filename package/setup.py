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
      install_requires=[
          'PeptideBuilder@https://api.github.com/repos/mtien/PeptideBuilder/tarball/master', 
          'GromacsWrapper@https://api.github.com/repos/whitead/GromacsWrapper/tarball/master', 
          'dill', 
          'traitlets', 
          'requests', 
          'future', 
          'pytest', 
          'biopython', 
          'pandas', 
          'tqdm'],
      test_suite='tests',
      zip_safe=True,
      package_data={'peptidesim': ['templates/*.mdp', 'templates:/*db']}
)

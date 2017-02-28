import os
from glob import glob
from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

exec(open('peptidesim/version.py').read())

setup(name = 'peptidesim', 
      version = __version__,
      scripts = glob(os.path.join('scripts', '*')),
      description = 'An automated peptide simulator',
      long_description=readme(),
      author = 'Dilnoza Amirkulova, Andrew White', 
      author_email = 'andrew.white@rochester.edu', 
      url = 'http://thewhitelab.org/Software',
      license = 'GPL3',
      packages = ['peptidesim'],
      install_requires=['PeptideBuilder==1.0.1', 'GromacsWrapper>=0.5.1-dev', 'dill', 'traitlets', 'requests', 'future'],
      dependency_links=['https://github.com/mtien/PeptideBuilder/archive/master.zip#egg=PeptideBuilder-1.0.1', 'https://github.com/Becksteinlab/GromacsWrapper/archive/develop.zip#egg=GromacsWrapper-0.5.1-dev'],
      test_suite='tests',
      zip_safe=True,
      package_data = {'peptidesim': ['templates/*.mdp']}
)

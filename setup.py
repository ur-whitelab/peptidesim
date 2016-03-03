from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name = 'peptidesim', 
      version = '0.1', 
      description = 'An automated peptide simulator',
      long_description=readme(),
      author = 'Dinloza Amirkulova, Andrew White', 
      author_email = 'andrew.white@rochester.edu', 
      url = 'http://thewhitelab.org/Software',
      license = 'GPL3',
      packages = ['peptidesim'],
      install_requires=['https://github.com/mtien/PeptideBuilder/master#egg=package-1.0.1', 'GromacsWrapper'],
      zip_safe=True
)

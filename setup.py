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
      install_requires=['PeptideBuilder==1.0.1', 'GromacsWrapper== 0.4.0'],
      dependency_links=['https://github.com/mtien/PeptideBuilder/archive/master.zip#egg=PeptideBuilder-1.0.1', 'https://github.com/Becksteinlab/GromacsWrapper/archive/develop.zip#egg=GromacsWrapper-0.4.1'],
      test_suite='tests',
      zip_safe=True
)

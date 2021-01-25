from distutils.core import setup

from setuptools import find_packages

from storage import __version__

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: Apache Software License
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3.8
Operating System :: POSIX :: Linux
""".strip().split('\n')

setup(name='staramr',
      version=__version__,
      description='Indexes genomes.',
      author='Aaron Petkau',
      author_email='aaron.petkau@gmail.com',
      url='https://github.com/apetkau/thesis-index',
      license='Apache v2.0',
      classifiers=classifiers,
      install_requires=[
          'biopython>=1.70',
          'pandas>=1.0.0',
          'bitarray',
          'pyvcf',
          'pytest',
      ],
      packages=find_packages(),
      include_package_data=True,
      )

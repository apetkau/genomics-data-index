from distutils.core import setup

from setuptools import find_packages

from genomics_data_index import __version__

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

setup(name='thesis-index',
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
          'pyvcf',
          'sqlalchemy',
          'pymysql',
          'requests',

          # ete3 requires PyQt5 for some functionality <https://github.com/etetoolkit/ete/issues/354>
          'ete3',
          'PyQt5',

          'ga4gh.vrs[extras]==0.6.2',
          'biocommons.seqrepo',
          'click',
          'click-config-file',
          'coloredlogs',
          'pyyaml',
          'pyroaring',
          'sourmash',
          'pybedtools',
          'pytest',
      ],
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['variants=genomics_data_index.cli.variants:main'],
      },
      )

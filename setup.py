from distutils.core import setup

from setuptools import find_packages

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

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='genomics-data-index',
      version='0.9.0',
      description='Indexes genomics data (mutations, kmers, MLST) for fast querying of features.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Aaron Petkau',
      author_email='aaron.petkau@gmail.com',
      url='https://github.com/apetkau/genomics-data-index',
      license='Apache v2.0',
      classifiers=classifiers,
      install_requires=[
          'biopython>=1.70',
          'pandas>=1.0.0',
          'numpy',
          'scikit-bio',

          # Need to restrict scipy due to this issue in scikit-bio (https://github.com/biocore/scikit-bio/issues/1818)
          'scipy<1.9',

          'vcfpy',
          'sqlalchemy',
          'pymysql',
          'requests',

          # ete3 requires PyQt5 for some functionality <https://github.com/etetoolkit/ete/issues/354>
          'ete3',
          'PyQt5',

          'biocommons.seqrepo',
          'click',
          'click-config-file',
          'coloredlogs',
          'pyyaml',
          'pyroaring',
          'sourmash',
          'pybedtools',
          'snakemake',
          'Jinja2',
          'pathvalidate',
          'pytest',
          'zipp',
          'packaging',
      ],
      # I need to restrict to less then 3.10 now due to this issue in ete3
      # https://github.com/etetoolkit/ete/issues/635
      # This stems from issues in PyQt5
      python_requires=">=3.8,<3.10",
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['gdi=genomics_data_index.cli.gdi:main'],
      },
      )

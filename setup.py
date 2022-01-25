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
      version='0.5.0.dev7',
      description='Indexes genomics data (mutations, kmers, MLST) for fast querying of features.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Aaron Petkau',
      author_email='aaron.petkau@gmail.com',
      url='https://github.com/apetkau/genomics-data-index',
      license='Apache v2.0',
      classifiers=classifiers,
      install_requires=[
          # pyvcf uses the option use_2to3 which is not compatible with setuptools>=58
          'setuptools<58',

          'biopython>=1.70',
          'pandas>=1.0.0',
          'numpy',
          'scikit-bio',
          'scipy',
          'pyvcf',
          'sqlalchemy',
          'pymysql',
          'requests',

          # ete3 requires PyQt5 for some functionality <https://github.com/etetoolkit/ete/issues/354>
          'ete3',
          'PyQt5',

          # ga4gh cannot use recent versions of jsonschema
          'jsonschema==3.2.0',
          'ga4gh.vrs[extras]==0.6.2',
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

          # I do not know exactly why these two statements are required and aren't 
          # picked up by the dependency resolver, but they seem to be necessary now.
          # 'tomli' required by 'black' and the wrong version is being installed. This fixes it.
          'tomli<2.0.0,>=0.2.6',
          'zipp',
      ],
      python_requires="<3.9",
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['gdi=genomics_data_index.cli.gdi:main'],
      },
      )

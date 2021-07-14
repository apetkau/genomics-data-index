# Genomics data index
[![Build Status](https://github.com/apetkau/thesis-index/workflows/Integration%20Tests/badge.svg?branch=development)](https://github.com/apetkau/thesis-index/actions?query=branch/development)
[![pypi](https://badge.fury.io/py/genomics-data-index.svg)](https://pypi.python.org/pypi/genomics-data-index/)
[![Binder](https://mybinder.org/badge_logo.svg)][tutorial1-binder]

This project is to design a system which can index large amounts of genomics data and enable rapid querying of this
data.

**Indexing** breaks genomes up into individual features (*nucleotide mutations*, *k-mers*, or *genes/MLST*)
and stores the index in a directory which can easily be shared with other people. Indexes can be generated 
direct from sequence data or loaded from existing intermediate files (e.g., VCF files, MLST results).

```bash
# Analyze sequence data (reads/assemblies, compressed/uncompressed)
gdi analysis --reference-file genome.gbk.gz *.fasta.gz *.fastq.gz

# (Alternatively) Index features in previously computed VCF files listed in vcf-files.txt
gdi load vcf vcf-files.txt
```

**Querying** provides both a *Python API* and *Command-line interface* to select sets of samples using this index
or attached external data (e.g., phylogenetic trees or DataFrames of metadata).

```python
# Select samples with a 26568 C > A mutation
r = s.hasa('MN996528.1:26568:C:A')

# Select samples with a D614G mutation on gene S
r = s.hasa('hgvs:MN996528.1:S:D614G')
```

**Summaries of the features** (mutations, kmers, MLST) can be exported from a set of samples alongside
*nucleotide alignments*, *distance matrices* or *trees* constructed from subsets of features.

```python
r.summary_features()
```

| Mutation | Count |
|----------|-------|
| 10 G>T   | 1     |
| 20 C>T   | 3     |
| 30 A>G   | 5     |

**Visualization** of trees and sets of selected samples can be constructed using the provided Python API and the
visualization tools provided by the [ETE Toolkit][].

```python
r.tree_styler() \
 .highlight(set1) \
 .highlight(set2) \
 #...
 .render()
```

![tree-visualization.png][]

You can see more examples of this software in action in the provided [Tutorials](#5-tutorial).

# Table of contents

- [1. Overview](#1-overview)
  * [1.1. Indexing](#11-indexing)
    + [1.1.1. Naming features](#111-naming-features)
  * [1.2. Querying](#12-querying)
    + [1.2.1. Python API](#121-python-api)
- [2. Background](#2-background)
- [3. Installation](#3-installation)
  * [3.1. Conda](#31-conda)
  * [3.2. PyPI/pip](#32-pypipip)
  * [3.3. From GitHub](#33-from-github)
- [4. Usage](#4-usage)
  * [4.1.1](#411-indexing)
  * [4.1.2](#412-querying)
  * [4.1.3](#413-main-usage-statement)
- [5. Tutorial](#5-tutorial)
- [6. Acknowledgements](#6-acknowledgements)

# 1. Overview

The software is divided into two main components: *(1) Indexing* and *(2) Querying*.

## 1.1. Indexing

![figure-index.png][]

The *indexing* component provides a mechanism to break genomes up into individual features and store these features
in a database. The types of features supported include: **Nucleotide mutations**, **K-mers**, and **Genes/MLST**.

### 1.1.1. Naming features

Indexing assigns names to the individual features, represented as strings inspired by the
[Sequence Position Deletion Insertion (SPDI)][] model.

1. **Nucleotide mutations**: `sequence:position:deletion:insertion` (e.g., `ref:100:A:T`)
2. **Genes/MLST**: `scheme:locus:allele` (e.g., `ecoli:adk:100`)
3. **Kmers**: *Not implemented yet*

Alternatively, for **Nucleotide mutations** names can be given using [hgvs](https://varnomen.hgvs.org/) (as output by [SnpEff](http://pcingola.github.io/SnpEff/)).

1. **Nucleotide mutations**: `hgvs:sequence:gene:p.protein_change` (e.g., `hgvs:ref:geneX:p.P20H`).

## 1.2. Querying

![figure-queries.png][]

The *querying* component provides a Python API or command-line interface for executing queries on the genomics index. The primary type of query is a 
**Samples query** which returns sets of samples based on different criteria. These criteria
are grouped into different **Methods**. Each method operates on a particular type of **Data**
which could include features stored in the *genomics index* as well as trees or external metadata.

### 1.2.1. Python API

An example query on an existing set of samples `s` would be:

```python
r = s.isa('B.1.1.7', isa_column='lineage') \
     .isin(['SampleA'], distance=1, units='substitutions') \
     .hasa('MN996528.1:26568:C:A')
```

This would be read as:

>Select all samples in `s` which are a **B.1.1.7 lineage** as defined in some attached DataFrame (`isa()`) *AND*
> which are within **1 substitution** of **SampleA** as defined on a phylogenetic tree (`isin()`) *AND* 
> which have a **MN996528.1:26568:C:A** mutation (`hasa()`).

*Note: I have left out some details in this query. Full examples for querying are available at [Tutorial 1: Salmonella dataset][tutorial1].*

# 2. Background

This is still an ongoing project. A lot of background material is found in my [Thesis proposal][thesis-proposal].

# 3. Installation

## 3.1. Conda

[Conda][conda] is a package and environment management software which makes it very easy to install and maintain dependencies of software without requiring administrator/root access. Packages from conda are provided using different *channels* and the [bioconda][] channel contains a very large collection of bioinformatics software which can be automatically installed. To make use of conda you will have to first download and [install conda][]. Once installed you can use the command `conda` to install software/manage conda environments.

To install this software, we will first, create a conda environment with the necessary dependencies as follows (a full conda package is not available yet https://github.com/apetkau/genomics-data-index/issues/51 ).

```bash
conda create -c bioconda -c conda-forge -c defaults --name gdi python=3.8 pyqt bedtools iqtree bcftools htslib

# Activate environment. Needed to install additional Python dependencies below.
conda activate gdi
```

Now, you can install with:

```bash
pip install genomics-data-index
```

If everything is working you should be able to run:

```bash
gdi --version
```

You should see `gdi, version 0.1.0` printed out.

### Additional dependencies

For snpeff to work you will need to install the package `mkisofs` on Ubuntu (e.g., `sudo apt install mkisofs`).
I do not know the exact package name on other systems.

## 3.2. PyPI/pip

To install just the Python component of this project from [PyPI][gdi-pypi] you can run the following:

```bash
pip install genomics-data-index
```

Note that you will have to install some additional dependencies separately. Please see the [conda-env.yaml][] environment file for details.

## 3.3. From GitHub

To install the project from the source on GitHub please first clone the git repository:

```bash
git clone https://github.com/apetkau/genomics-data-index.git
cd genomics-data-index
```

Now install all the dependencies using [conda][] and [bioconda][] with:

```bash
conda env create -f conda-env.yaml
conda activate gdi
```

Once these are installed you can setup the Python package with:

```bash
pip install .
```

# 4. Usage

The main command is called `gdi`. A quick overview of the usage is as follows:

## 4.1. Indexing

```bash
# Create new index in `index/`
# cd to `index/` to make next commands easier to run
gdi init index
cd index

# Creates an index of mutations (VCF files) and kmer sketches (sourmash)
gdi analysis --use-conda --include-kmer --kmer-size 31 --reference-file genome.gbk.gz *.fastq.gz

# (Optional) build tree from mutations (against reference genome `genome`) for phylogenetic querying
gdi rebuild tree --align-type full genome
```

The produced index will be in the directory `index/`.

## 4.2. Querying

```bash
# List indexed samples
gdi list samples

# Query for genomes with mutation
gdi query mutation 'genome:10:A:T'
```

## 4.3. Main usage statement

```
Usage: gdi [OPTIONS] COMMAND [ARGS]...

Options:
  --project-dir TEXT              A project directory containing the data and
                                  connection information.

  --ncores INTEGER RANGE          Number of cores for any parallel processing
                                  [default: 8]

  --log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  Sets the log level  [default: INFO]
  --version                       Show the version and exit.
  --config FILE                   Read configuration from FILE.
  --help                          Show this message and exit.

Commands:
  analysis
  build
  db
  export
  init
  input
  list
  load
  query
  rebuild
```

# 5. Tutorial

Tutorials and a demonstration of the software are available below (code in [separate repository][gdi-examples]). You can select the **[launch | binder]** badge to launch each of these tutorials in an interactive [Jupyter][] environment within the cloud environment using [Binder][].

1. [Tutorial 1: Querying (Salmonella)][tutorial1] - [![Binder](https://mybinder.org/badge_logo.svg)][tutorial1-binder] 
    * _In case GitHub link is not rendering try [here][tutorial1-nbviewer]_
2. [Tutorial 2: Indexing assemblies (SARS-CoV-2)][tutorial2] - [![Binder](https://mybinder.org/badge_logo.svg)][tutorial2-binder]
    * _In case GitHub link is not rendering try [here][tutorial2-nbviewer]_
3. [Tutorial 3: Querying overview][tutorial3] - [![Binder](https://mybinder.org/badge_logo.svg)][tutorial3-binder]
    * _In case GitHub link is not rendering try [here][tutorial3-nbviewer]_

Alternatively, you can run these tutorials on your local machine. In order to run these tutorials you will first have to install the `genomics-data-index` software (see the [Installation](#3-installation) section for details). In addition, you will have to install [Jupyter Lab][]. If you have already installed the `genomics-data-index` software with conda you can install Jupyter Lab as follows:

```bash
conda activate gdi
conda install jupyterlab
```

To run Jupyter you can run the following:

```bash
# QT_QPA_PLATFORM The below is useful to avoid having to set the DISPLAY env variable for Qt
# You can ignore setting this environment variable if you are running on a machine with an X server installed and configured

QT_QPA_PLATFORM="offscreen" jupyter lab
```

Please see the instructions for [Jupyter Lab][jupyter-docs] for details.

# 6. Acknowledgements

I would like to acknowledge the [Public Health Agency of Canada][], the [University of Manitoba][], 
and the [VADA Program][] for providing me with the opportunity, resources and training for working on this project.

Some icons used in this documentation are provided by [Font Awesome][] and licensed under a 
[Creative Commons Attribution 4.0][] license.

[thesis-proposal]: https://drive.google.com/file/d/1sd0WjmwO_KU5wacfpUiPGT20xVOwBc8i/view?usp=sharing
[gdi-examples]: https://github.com/apetkau/genomics-data-index-examples
[BIGSI]: https://bigsi.readme.io/
[COBS]: https://github.com/bingmann/cobs
[ci-dependencies]: .github/workflows/ci-test.yml#L37
[tutorial1]: https://github.com/apetkau/genomics-data-index-examples/blob/main/examples/tutorial-1-salmonella.ipynb
[tutorial2]: https://github.com/apetkau/genomics-data-index-examples/blob/main/examples/tutorial-2-sars-cov-2.ipynb
[tutorial3]: https://github.com/apetkau/genomics-data-index-examples/blob/main/examples/tutorial-3-querying-overview.ipynb
[snippy]: https://github.com/tseemann/snippy
[ETE Toolkit]: http://etetoolkit.org/
[tree-visualization.png]: docs/images/tree-visualization.png
[figure-index.png]: docs/images/figure-index.png
[figure-queries.png]: docs/images/figure-queries.png
[Public Health Agency of Canada]: https://www.canada.ca/en/public-health.html
[University of Manitoba]: https://umanitoba.ca/
[VADA Program]: http://vada.cs.umanitoba.ca/
[conda]: https://docs.conda.io/en/latest/
[install conda]: https://docs.conda.io/en/latest/miniconda.html
[bioconda]: https://bioconda.github.io/
[Font Awesome]: https://fontawesome.com/
[Creative Commons Attribution 4.0]: https://fontawesome.com/license/free
[Sequence Position Deletion Insertion (SPDI)]: https://doi.org/10.1093/bioinformatics/btz856
[pypi-gdi]: https://pypi.org/project/genomics-data-index/
[conda-env.yaml]: conda-env.yaml
[Jupyter Lab]: https://jupyter.org/
[Jupyter]: https://jupyter.org/
[jupyter-docs]: https://jupyterlab.readthedocs.io/en/latest/
[Binder]: https://mybinder.org/
[tutorial1-nbviewer]: https://nbviewer.jupyter.org/github/apetkau/genomics-data-index-examples/blob/main/examples/tutorial-1-salmonella.ipynb
[tutorial2-nbviewer]: https://nbviewer.jupyter.org/github/apetkau/genomics-data-index-examples/blob/main/examples/tutorial-2-sars-cov-2.ipynb
[tutorial3-nbviewer]: https://nbviewer.jupyter.org/github/apetkau/genomics-data-index-examples/blob/main/examples/tutorial-3-querying-overview.ipynb
[tutorial1-binder]: https://mybinder.org/v2/gh/apetkau/genomics-data-index-examples/main?urlpath=lab%2Ftree%2Fexamples%2Ftutorial-1-salmonella.ipynb
[tutorial2-binder]: https://mybinder.org/v2/gh/apetkau/genomics-data-index-examples/main?urlpath=lab%2Ftree%2Fexamples%2Ftutorial-2-sars-cov-2.ipynb
[tutorial3-binder]: https://mybinder.org/v2/gh/apetkau/genomics-data-index-examples/main?urlpath=lab%2Ftree%2Fexamples%2Ftutorial-3-querying-overview.ipynb
[gdi-pypi]: https://pypi.python.org/pypi/genomics-data-index/


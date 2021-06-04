# Genomics data index
[![Build Status](https://github.com/apetkau/thesis-index/workflows/Integration%20Tests/badge.svg?branch=development)](https://github.com/apetkau/thesis-index/actions?query=branch/development)
[![pypi](https://badge.fury.io/py/genomics-data-index.svg)](https://pypi.python.org/pypi/genomics-data-index/)

This project is to design a system which can index large amounts of genomics data and enable rapid querying of this
data.

**Indexing** breaks genomes up into individual features (*nucleotide mutations*, *k-mers*, or *genes/MLST*)
and stores the index in a directory which can easily be shared with other people. Indexes can be generated 
direct from sequence data or loaded from existing intermediate files (e.g., VCF files).

```bash
# Index features in VCF files listed in vcf-files.txt
gdi load vcf vcf-files.txt
```

**Querying** provides both a *Python API* and *Command-line interface* to select sets of samples using this index
or attached external data (e.g., phylogenetic trees or DataFrames of metadata).

```python
# Select samples with a 26568 C > A mutation
r = s.hasa('MN996528.1:26568:C:A')
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
  * [3.1. PyPI/pip](#31-pypipip)
  * [3.2. From source](#32-from-source)
- [4. Usage](#4-usage)
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

1. **Nucleotide mutations**: `reference:position:deletion:insertion` (e.g., `ref:100:A:T`)
2. **Genes/MLST**: `scheme:locus:allele` (e.g., `ecoli:adk:100`)
3. **Kmers**: *Not implemented yet*

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

## 3.1. PyPI/pip

To install this project from [PyPI][gdi-pypi] you can run the following:

```bash
pip install genomics-data-index
```

Note that you will have to install some additional dependencies separately. Please see the [conda-env.yaml][] environment file for details.

## 3.2. From source

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

The main command is called `gdi`:

```
Usage: gdi [OPTIONS] COMMAND [ARGS]...

Options:
  --project-dir TEXT              A project directory containing the data and
                                  connection information.

  --ncores INTEGER RANGE          Number of cores for any parallel processing
                                  [default: 8]

  --log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  Sets the log level  [default: INFO]
  --config FILE                   Read configuration from FILE.
  --help                          Show this message and exit.

Commands:
  build
  db
  export
  init
  list
  load
  query
  rebuild
```

If you've previously ran [snippy][] you can load up this data (i.e., the SNPs in VCF format) as follows:

```bash
# Initialize and cd to project
gdi init proj1
cd proj1

# Load data
gdi load snippy --reference-file reference.fasta snippy-analysis/
```

Where `snippy-analysis/` contains directories like `SampleA`, `SampleB`, etc.

# 5. Tutorial

Tutorials and a demonstration of the software is available at:

1. [Tutorial 1: Querying (Salmonella)][tutorial1]
2. [Tutorial 2: Indexing assemblies (SARS-CoV-2)][tutorial2]

# 6. Acknowledgements

I would like to acknowledge the [Public Health Agency of Canada][], the [University of Manitoba][], 
and the [VADA Program][] for providing me with the opportunity, resources and training for working on this project.

Some icons used in this documentation are provided by [Font Awesome][] and licensed under a 
[Creative Commons Attribution 4.0][] license.

[thesis-proposal]: https://drive.google.com/file/d/1sd0WjmwO_KU5wacfpUiPGT20xVOwBc8i/view?usp=sharing
[BIGSI]: https://bigsi.readme.io/
[COBS]: https://github.com/bingmann/cobs
[ci-dependencies]: .github/workflows/ci-test.yml#L37
[tutorial1]: docs/tutorial/tutorial-1-salmonella.ipynb
[tutorial2]: docs/tutorial/tutorial-2-sars-cov-2.ipynb
[snippy]: https://github.com/tseemann/snippy
[ETE Toolkit]: http://etetoolkit.org/
[tree-visualization.png]: docs/images/tree-visualization.png
[figure-index.png]: docs/images/figure-index.png
[figure-queries.png]: docs/images/figure-queries.png
[Public Health Agency of Canada]: https://www.canada.ca/en/public-health.html
[University of Manitoba]: https://umanitoba.ca/
[VADA Program]: http://vada.cs.umanitoba.ca/
[conda]: https://docs.conda.io/en/latest/
[bioconda]: https://bioconda.github.io/
[Font Awesome]: https://fontawesome.com/
[Creative Commons Attribution 4.0]: https://fontawesome.com/license/free
[Sequence Position Deletion Insertion (SPDI)]: https://doi.org/10.1093/bioinformatics/btz856
[pypi-gdi]: https://pypi.org/project/genomics-data-index/


[![Build Status](https://github.com/apetkau/thesis-index/workflows/Integration%20Tests/badge.svg?branch=development)](https://github.com/apetkau/thesis-index/actions?query=branch/development)

# Genomics data index

This project is to design a system which can index large amounts of genomics data and enable rapid querying of this
data. This is an ongoing (in-development) project and so not all (or even most) of the ideas below are implemented yet.

A short tutorial and example of this software is available at [Tutorial 1: Salmonella dataset][tutorial1].

## Background

The type of genomic features that I will index are:

1. Single nucleotide variants/mutations (and small indels)
2. Kmers
3. Genes (e.g., MLST)

The goal is to abstract out the common methods to index these different types of features into a single interface.

The types of queries I wish to perform are as follows:

1. Sample-based queries (or sample-relatedness queries)

   These could consist of questions such as "list samples related to Sample X"

2. Feature-based queries

   These could consist of questions such as "list all samples with an A -> T mutation at position 500" or "list all
   samples containing gene xyz (query based on k-mers similar to [BIGSI][] or [COBS][]).

3. Cluster-based queries

   These could consist of questions such as "list all samples in cluster X" (where cluster X is defined based on shared
   features or placement in a phylogenetic tree).

You can see more details in my [Thesis proposal][thesis-proposal].

# Dependencies

This project requires Python, MariaDB (or some other relational database) and the following other software:

1. bcftools
2. bedtools
3. sourmash
4. iqtree

# Installation

It's best to install everything in a conda environment:

```bash
conda env create -f conda-env.yaml
conda activate gdi
```

This is an on-going project which will undergo a lot of changes and so the exact conda dependnecies are probably best
found in the [CI test script][ci-dependencies] (for now).

Once these are installed you can setup the Python package with:

```bash
pip install .
```

# Usage

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

# Tutorial

A tutorial and demonstration of the software is available at:

1. [Tutorial 1: Salmonella dataset][tutorial1]


[thesis-proposal]: https://drive.google.com/file/d/1sd0WjmwO_KU5wacfpUiPGT20xVOwBc8i/view?usp=sharing
[BIGSI]: https://bigsi.readme.io/
[COBS]: https://github.com/bingmann/cobs
[ci-dependencies]: .github/workflows/ci-test.yml#L37
[tutorial1]: docs/tutorial/tutorial-1-salmonella.ipynb
[snippy]: https://github.com/tseemann/snippy



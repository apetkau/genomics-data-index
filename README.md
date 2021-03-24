# A framework for the indexing and clustering of pathogen genomes

This project is to design a system which can index large amounts of genomics data and enable rapid querying of this data.
This is an ongoing (in-development) project and so not all (or even most) of the ideas below are implemented yet.

The type of genomic features that I will index are:

1. Single nucleotide variants/mutations (and small indels)
2. Kmers
3. Genes (e.g., MLST)

The goal is to abstract out the common methods to index these different types of features into a single interface.

The types of queries I wish to perform are as follows:

1. Sample-based queries (or sample-relatedness queries)
   
   These could consist of questions such as "list samples related to Sample X"

2. Feature-based queries

   These could consist of questions such as "list all samples with an A -> T mutation at position 500" or "list all samples containing gene xyz (query based on k-mers similar to [BIGSI][] or [COBS][]".

3. Cluster-based queries

   These could consist of questions such as "list all samples in cluster X" (where cluster X is defined based on shared features or placement in a phylogenetic tree).

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
conda create --name index python=3.8 htslib bcftools==1.12 bedtools fasttree iqtree
conda activate index
```

This is an on-going project which will undergo a lot of changes and so the exact conda dependnecies are probably best found in the [CI test script][ci-dependencies] (for now).

Once these are installed you can setup the Python package with:

```bash
pip install .
```

In addition, you will also have to setup a relational database (I've been using MariaDB).

Data is stored in two separate locations: in the database and in a directory on the filesystem.
To configure both and not have to pass command-line arguments all the time you can define a file `config.yaml`.

```yaml
database_connection: mysql+pymysql://test:test@localhost/thesis?charset=utf8mb4
database_dir: /home/apetkau/workspace/thesis-index/tmp/data
```

# Usage

The main command is called `variants`:

```bash
variants --config config.yaml --help
Usage: variants [OPTIONS] COMMAND [ARGS]...

Options:
  --database-connection TEXT  A connection string for the database.
  --database-dir PATH         The root directory for the database files.
  --verbose / --no-verbose    Turn up verbosity of command
  --config FILE               Read configuration from FILE.
  --help                      Show this message and exit.

Commands:
  alignment
  export
  list
  load-kmer
  load-snippy
  load-vcf
  query
  tree
```

If you've previously ran [snippy][] you can load up this data (i.e., the SNPs in VCF format) as follows:

```bash
variants --config config.yaml load-snippy --build-tree --align-type full --reference-file reference.fasta snippy-dir/
```

You can then query like:

```bash
 variants query --type mutation reference:528:C:CAG | column -s $'\t' -t
Type      Feature              Sample Name  Sample ID  Status
mutation  reference:528:C:CAG  SampleB      2          Present
```

[thesis-proposal]: https://drive.google.com/file/d/1sd0WjmwO_KU5wacfpUiPGT20xVOwBc8i/view?usp=sharing
[BIGSI]: https://bigsi.readme.io/
[COBS]: https://github.com/bingmann/cobs
[ci-dependencies]: https://github.com/apetkau/thesis-index/blob/development/.github/workflows/ci-test.yml#L37
[snippy]: https://github.com/tseemann/snippy
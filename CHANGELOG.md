# 0.4.0

* [analysis]: Switched all steps to use conda in Snakemake pipeline (0.3.1).
* [api]: Added `rotation`, `allow_face_overlap`, `show_branch_length`, and `show_branch_support` to tree viewing API (0.4.0.dev1).
* [api]: Fixed up `features_comparison` for joined dataframe to query so it uses proper subset of dataframe (0.4.0.dev1).
* [api]: Implemented a `prune()` method for pruning a tree down to only the selected samples (0.4.0.dev1).
* [cli]: Fixed bug where `--features-summary mlst` would not work in command-line interface (0.4.0.dev1).

# 0.3.0

* [api]: Added ability to generate a DataFrame of percents of features present in selected samples. To be used for comparing different categories of samples (0.3.0.dev1).
* [api]: Added ability to handle joining larger dataframes to sets of samples by batching SQL queries (0.3.0.dev2).
* [api]: Added mutation "Type" to the output of summaries/comparison table (0.3.0.dev3).
* [analysis]: Added ability to select to include MNP as well as SNPs when building an alignment/tree (0.3.0.dev2).
* [analysis]: Santizing sample names for analysis and restoring afterwards. This way a greater variety of sample names is possible (0.3.0.dev3, 0.3.0.dev12).
* [analysis]: Automatically split a single multi-FASTA file into separate files per sequence (used primarily for SARS-CoV-2 data) (0.3.0.dev4).
* [analysis]: Added ability to handle lzma and bzip2 compressed sequence files (0.3.0.dev4).
* [api]: Added a method to select a random subsample of a query (0.3.0.dev4).
* [api]: Set outgroup of tree query (0.3.0.dev5).
* [api]: Added experimental class to cluster samples within different categories (e.g., lineages) by distances between the proportion of samples having particular features (0.3.0.dev5).
* [analysis]: All stages of the Snakemake workflow now have assigned conda environments (0.3.0.dev5).
* [cli]: Added support for loading mask files as BED files in addition to sequence files (0.3.0.dev6).
* [analysis]: Fixed issue where there was overlap between unknown/missing positions and VCF files (0.3.0.dev7).
* [cli/api]: Added `--include-variants DELETION` and `--include-variants DELETION_OTHER` to include deletions in the alignment that is generated (represented by gaps `-`) (0.3.0.dev8).
* [api]: Added method to define custom `q.isa()` methods and implemented a method to type SARS-CoV-2 genomes using [constellations](https://github.com/cov-lineages/constellations) of mutations (the SARS-CoV-2 typer should be considered experimental at this stage and is not guaranteed to work in all cases) (0.3.0.dev9).
* [api]: Loading both VCFs and Kmer sketches listed in the same file (0.3.0.dev10).
* [analysis]: Changed default option to `--use-conda` for analysis.
* [cli]: Added general query command for multiple types of queries (`gdi query hasa:feature` or `gdi query isa:sample`) as well as summarizing features (`--features-summary`) (0.3.0.dev11).
* [cli]: Renamed `--summarize` to `--summary` for CLI (0.3.0.dev11).
* [install]: Restricted `setuptools<58` since `pyvcf` uses the option `use_2to3` which is no longer compatible with setuptools >= 58 <https://setuptools.readthedocs.io/en/latest/history.html#v58-0-0> (0.3.0.dev13).

# 0.2.0

* [doc]: Updates to readme and other documentation.
* [doc]: Moving tutorials to separate repository <https://github.com/apetkau/genomics-data-index-examples>.
    * Added method to launch tutorials using Binder.
* [api]: Adding ability to summarize only features unique to a selected set (for `q.features_summary()`).
* [api]: Adding `Total` and `Percent` columns to `q.features_summary()` dataframe.
* [analysis]: Added ability to load VCF files with [SnpEff](http://pcingola.github.io/SnpEff/) annotations and associte these with the nucleotide identifiers.
    * `features_summary()` and `hasa()` works with the SnpEff variant identifiers (HGVS).
    * SnpEff results will be included when loading assemblies.
* [analysis]: Added the ability to analyze and index sequence reads (by mapping to a reference genome using [snippy](https://github.com/tseemann/snippy)).
* [analysis]: Added ability to insert new genomes into an existing index.
* [analysis]: Added ability to index and query missing data.
    * The resulting samples for each query can be divided up into 3 categories (**True/Present**, **False/Absent**, **Unknown/Unknown**). Different sample sets track these different categories (`q.sample_set`, `q.unknown_set`, `q.absent_set`).
* [developer]: Added custom TRACE log level.
* [analysis]: Added parallel processing for construction of features table.
* [analysis]: Switched command-line interface to use my Python API for queries.
* [api]: `features_summary()` now works with MLST results.
* [api]: Improved performance/redesigned `features_summary()` for mutation results.
* [analysis]: Added ability to batch up loading of samples into database.

# 0.1.0

* Initial release of entire project.
    * Able to index assembled genomes.
    * Index by mutations (VCF files) and kmers.
    * Partial indexing of MLST.
* Adding method to assign a score to clusters in a tree based on how well sets of samples are grouped together.
* Changed up applying visual styles to sets of samples in a tree to only apply when rendering.
    * This fixes an issue where copying trees was crashing when rendering large trees (I copy only once using 'newick').
* Wrote a large amount of documentation (not complete yet).
* Wrote tutorials for usage of this software.

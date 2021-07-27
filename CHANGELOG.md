# 0.2.0

* [doc]: Updates to readme and other documentation.
* [doc]: Moving tutorials to separate repository <https://github.com/apetkau/genomics-data-index-examples>.
    * Added method to launch tutorials using Binder.
* [api]: Adding ability to summarize only features unique to a selected set (for `q.summary_features()`).
* [api]: Adding `Total` and `Percent` columns to `q.summary_features()` dataframe.
* [analysis]: Added ability to load VCF files with [SnpEff](http://pcingola.github.io/SnpEff/) annotations and associte these with the nucleotide identifiers.
    * `summary_features()` and `hasa()` works with the SnpEff variant identifiers (HGVS).
    * SnpEff results will be included when loading assemblies.
* [analysis]: Added the ability to analyze and index sequence reads (by mapping to a reference genome using [snippy](https://github.com/tseemann/snippy)).
* [analysis]: Added ability to insert new genomes into an existing index.
* [analysis]: Added ability to index and query missing data.
    * The resulting samples for each query can be divided up into 3 categories (**True/Present**, **False/Absent**, **Unknown/Unknown**). Different sample sets track these different categories (`q.sample_set`, `q.unknown_set`, `q.absent_set`).
* [developer]: Added custom TRACE log level.
* [analysis]: Added parallel processing for construction of features table.

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

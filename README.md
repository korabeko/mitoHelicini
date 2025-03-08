# mitoHelicini

## An R package containing a collection of geo-referrenced sequences of complete mitogenomes and partial mitochondrial genes of the land snail family Helicidae, tribe Helicini, and fuctions for curating and exporting the data.
maintained by Ondřej Korábek (Department of Zoology, Faculty of Science, Charles University, Praha, Czechia, ondrej.korabek@natur.cuni.cz)

### For the impatient
To export an alignment ready for analysis: run `help(representative.alignment)` and then try the function with different settings. To see the whole data, go to the /inst/extdata folder or use `export.Helicini()` to copy its contents to your working directory.

### The data
#### What are Helicini?
The tribe [Helicini](https://en.wikipedia.org/wiki/Helicini) is a clade of land snails from the Western Palearctic (Gastropoda: Heterobranchia: Pneumopulmonata: Stylommatophora: Helicidae: Helicini). It comprises 10 genera of large snails (shell diameter ~2-6 cm). The most diverse genus is [*Helix*](https://en.wikipedia.org/wiki/Helix_(gastropod)), which includes [*Helix pomatia*](https://en.wikipedia.org/wiki/Helix_pomatia) Linnaeus, 1758, one of the best studied land snail species.
#### What data are included?
The dataset is an updated and expanded version of that published by Korábek et al. 2022, *Diversity* 14: 24 [https://doi.org/10.3390/d14010024](https://doi.org/10.3390/d14010024). The original data from that paper are available from [Dryad](https://doi.org/10.5061/dryad.pnvx0k6p5). The data consist of nucleotide sequences of partial mitochondrial genes obtained by PCR and assembled sequences of complete and nearly complete mitogenomes. These are data either obtained by O. Korábek or downloaded from public databases (mostly [NCBI's Nucleotide database/GenBank](https://www.ncbi.nlm.nih.gov/nuccore), but also [BOLD](https://v4.boldsystems.org/) and [GBOL](https://gbol1.bolgermany.de/en/)). Only sequences accompanied by reasonably precise geographic location are included. Most of the data come from PCR targeting selected parts of the mitochondrial genome and followed by Sanger sequencing. For many individuals, one or more fragments of the mitogenome are thus available with gaps of unknown length in between. Complete mitogenome sequences, mostly assembled from Illumina reads (WGS or RNAseq) are available for a small (but increasing) proportion of the sequenced individuals.
The sequences are provided with complete metadata (locality, date of collection, shell vouchers, the person responsible for identification, etc.). Furthermore, an estimate of the phylogenetic tree of the included samples is provided, although not for all the individuals.
#### How do the date relate to records in the public databases?
As detailed in [Korábek et al. (2022)](https://doi.org/10.3390/d14010024), some of the sequences dowloaded from GenBank were edited (usually remains of primers or sequence ends suspected to be of a flow quality were trimmed, some positions with suspicious substitutions have been changed to ambiguity codes). In case of sequences submitted to GenBank by O. Korábek, the version deposited in GenBank may in some cases be shorter than the currently available sequence or they may be minor differences in the sequences. Updating GenBank records is a tedious and inflexible process. Furthermore, they sometimes do weird things there, like renaming sequences of Helix straminea as Biomphalaria straminea... It is therefore strongly recommended to use the latest version of sequence as provided here; the user may check which seqeuences differ from the GenBank record using the R function `GB.compare`. Finally, there are also sequence that were not yet submitted to GenBank that may be identified with the function `which.submit`. This pacakage thus provides collated, curated and more up-to-date data than you can get from GenBank.
Not all Helicini sequences from GenBank were incorporated. Some were ommited as their geographic origin was not known with enough precision or the sequence quality was doubted. These filters were not applied consistently. The goal is to mirror the current state of knowledge as far as possible in terms of geographic distribution and mitochondrial lineages, so a handful of somewhat dubiously localized records were included as these were unique sequences. Checking the comments collumn of the localities table is always recommended.
#### How are the sequence data structured?
Aligning the fragments takes some effort, especially for the rRNA genes (see the <version>_methods.odt file for alignment methods), so the data are structured to provide ready-to-use alignments of these genes in a simple fasta format.
Both complete and incomplete mitogenome sequences are split into the individual alignment files: full rRNA and protein-coding genes are stored together with incomplete sequences (the latter represent the vast majority), tRNA genes are kept only if there is a complete sequence spanning the region between two rRNA or protein-coding genes, which ensures some anchor for the sequence ends here, since there are often some non-conserved nucleotides between genes or at ends of the tRNAs. This solution allows to recover contiguous sequences by concatenation. There is an exception to this definition of the individual tRNA files. Due to the position of primer used to amplify a fragment of 12S rRNA within trnM, the 12S gene is kept in one file with trnM and trnM is not required to be complete if part of the 12S fragment amplified by PCR.
Alignment files for the protein-coding genes should start and end at positions homologous between species, but this is not always possible: see comments for the individual files below! Annotation file is available, indicating the start and end positions of the individual alignment files with respect to a complete mitogenome sequence of Helix pomatia (GenBank MK347426).
All alignments are stored in the direction of the heavy/plus strand of the reference mitogenome of Helix pomatia, i.e. some genes have to be reverse-complemented before translation (and 12S is also coded on the light/minus strand!). Similarly, 5' and 3' ends in the list of files below refer to the plus strand, not the direction of the gene if coded on the light strand.
Contents of each alignment file are described below. File names begin with the date of the last major revision.





### The functions
The package provides a set of utilities used to curate and export the data provided with the package and some functions used for analysis of the data in published studies.
The most important function is `representative.alignment`, which helps to export phylogenetically representative (at the level of intraspecific clades, species, or genera) or complete (all individuals) datasets, exporting a concantenate of selected loci or complete mitogenome sequences. As the data provided in the package are aligned, the exported dataset can be radily used for analyses.

### Do you want to contribute?

### Funding



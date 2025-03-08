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
As detailed in [Korábek et al. (2022)](https://doi.org/10.3390/d14010024), some of the sequences dowloaded from GenBank were edited (usually remains of primers or sequence ends suspected to be of a flow quality were trimmed, some positions with suspicious substitutions have been changed to ambiguity codes). In case of sequences submitted to GenBank by O. Korábek, the version deposited in GenBank may in some cases be shorter than the currently available sequence or they may be minor differences in the sequences. Updating GenBank records is a tedious and inflexible process. Furthermore, they sometimes do weird things there, like renaming sequences of *Helix straminea* as *Biomphalaria straminea*... It is therefore strongly recommended to use the latest version of sequence as provided here; the user may check which seqeuences differ from the GenBank record using the R function `GB.compare`. Finally, there are also sequence that were not yet submitted to GenBank that may be identified with the function `which.submit`. This pacakage thus provides collated, curated, and more up-to-date data with better taxonomy than you can get from GenBank.

Not all Helicini sequences from GenBank were incorporated. Some were ommited as their geographic origin was not known with enough precision or the sequence quality was doubted. These filters were not applied consistently. The goal is to mirror the current state of knowledge as far as possible in terms of geographic distribution and mitochondrial lineages, so a handful of somewhat dubiously localized records were included as these were unique sequences. Checking the comments collumn of the localities table is always recommended.
#### How are the sequence data structured?
Aligning the sequences takes some effort, especially for the rRNA genes, and is much easier to be done peacemeal. Moreover, from most individuals only short gene fragments were available. Therefore the data are stored in FASTA alignment files separately for individual protein-coding and rRNA genes and for each interval between these (these intervals largely consist of tRNA genes). That means that complete mitogenome sequences are also split into the individual alignment files. The intervals between PCGs and rRNAs are only kept for a given individual if there is a complete sequence spanning the region between two rRNA or protein-coding genes. That ensures some anchor for the sequence ends here, since there are often some non-conserved nucleotides between genes or at ends of the tRNAs. There is an exception to this definition of the individual alignment files. Due to the position of primer used to amplify a fragment of *rrnS* (12S rRNA) within *trnM*, the 12S gene is kept in one file with *trnM* and *trnM* is not required to be complete if part of the 12S fragment amplified by PCR.

Alignment files for the protein-coding genes should start and end at positions homologous between species, but this is not always possible: see comments for the individual files below! Annotation file  in gff3 format is available, indicating the start and end positions of the individual alignment files with respect to a complete mitogenome sequence of *Helix pomatia* (GenBank MK347426).

All alignments are stored in the direction of the heavy/plus strand of the reference mitogenome of *Helix pomatia*, i.e. some genes have to be reverse-complemented before translation (and *rrnS* is also coded on the light/minus strand!). Similarly, 5' and 3' ends in the list of files below refer to the plus strand, not the direction of the gene if coded on the light strand.

Contents of each alignment file are described below. File names begin with the date of the last major revision.
#### How are the sequences linked to metadata?
Each individual is assigned a unique ID (8 characters, eg. "HE003187"). This ID appears at the beginning of the fasta header of each sequence. It may be followed by another identifier separated by "_" (isolate lab code for original data, GenBank accession number for downloaded data - accession of one gene is usually used for all genes from the same individual and "barcode" if from a barcoding dataset not available from GenBank), species name, country (two-letter ISO code), and locality name (for internal use, not standardized). The ID stands for an individual snail, so it is identical between alignment files where the sequences come from the same snail. The ID provides a link between the sequences and the metadata.
#### How were the sequences aligned?
(to be added...)

### The metadata
(to be added...)

### The files
list of files (only name stem is given in square brackets here, which is preceded by version in the file names):
# alignment files (plain text in fasta format):
[Helicini_COX1]: COX1, mostly incomplete.
[Helicini_TRNV]: complete region between COX1 and 16S alignments; only for individuals with COX1 3' complete and 16S 5' complete.
[Helicini_16S]: 16S rRNA; because the exact start of the gene us uncertain all nucleotides after trnV are included.
[Helicini_TRNL1_TRNA]: complete region between 16S and ND6 alignments, only for individuals with 16S 3' complete and ND6 5' complete.
[Helicini_ND6]: ND6; gene ends with a T preceding the polyA tail even where this T follows after a complete stop codon.
[Helicini_TRNP]: complete region between ND6 and ND5 alignments, only for individuals with ND6 3' complete and ND5 5' complete
[Helicini_ND5]: ND5; ends always with a full stop codon, because ND1 and ND4L are also translated from the same transcript.
[Helicini_ND1_ND4L]: ND1 and ND4L; ND1 is incomplete at the 5' end due to overlap with the preceding ND5; the stop codon of ND1 usually overlaps the start codon of ND4L; would include also any nucleotides between ND4L and CYTB.
[Helicini_CYTB]: CYTB, mostly incomplete.
[Helicini_TRND-TRNF] complete region between CYTB and COX2 alignments, only for individuals with CYTB 3' complete and COX2 5' complete.
[Helicini_COX2]: COX2, mostly incomplete.
[Helicini_TRNY-TRNL2] complete region between COX2 and ATP8 alignments, only for individuals with COX2 3' complete and ATP8 5' complete.
[Helicini_ATP8]: ATP8.
[Helicini_TRNN] complete region between ATP8 and ATP6 alignments, only for individuals with ATP8 3' complete and ATP6 5' complete.
[Helicini_ATP6]: ATP6.
[Helicini_TRNR_TRNE]: complete region between ATP6 and 12S alignments, only for individuals with ATP6 3' complete and 12S alignment 5' complete.
[Helicini_12S_TRNM]: 12S rRNA and TRNM, mostly incomplete; the end position of the gene (5' end of the sequence here) is uncertain and therefore the TRNR_TRNE alignment probably includes part of 12S rRNA; because the 12S amplicon used in PCR-based studies extends into TRNM, TRNM is kept either up to the primer or, if complete, up to the 5' end of ND3.
[Helicini_ND3]: ND3.
[Helicini_TRNS2_TRNT]: complete region between COX3 and ND3 alignments, only for individuals with COX3 3' complete and ND3 5' complete.
[Helicini_COX3]: COX3; only the T (A at the beginning the sequence, because it is in reverse complement as the gene is on the minus strand) preceeding the polyA tail is homologous between species, even though in some species this is followed by AA, seemingly forming a full stop codon.
[Helicini_TRNS1]: complete region between ND3 and ND4 alignments, only for individuals with ND3 3' complete and ND4 5' complete; probably includes the control region in the 5' part.
[Helicini_ND4]: ND4; it covers the whole CDS and all the nucleotides preceeding the polyA tail.
[Helicini_TRNI]: complete region between ND4 and ND2 alignments, only for individuals with ND4 3' complete and ND2 5' complete.
[Helicini_ND2]: ND2.
[Helicini_TRNK]: complete region between ND2 and COX1 alignments, only for individuals with ND2 3' complete and COX1 5' complete.
[thessalica_16S]: partial 16S sequences of Helix thessalica and H. pomatia from a study of their hybrid zone in southern Czechia (Korábek & Hausdorf 2024 Molecular Ecology).

# sequence metadata (Open Document spread sheet)
[Helicini_table_mol]: species identifications, sampling locality data, NCBI accession numbers, shell vouchers, and other relevant metadata relating to the origin of the sequences and their use in publications; detailed explanation of columns inside the file.

# mitogenome annotations (plain text file)
[mitogenome_annotations]: annotation of mitogenomic sequences in the GFF3 format, including the positions where the the individual alignment files align. The current annotation is indicated by "current" in the "source" field.

# change log (plain text file)
[change_log]: record of revisions, in particular to the alignment files.


### The functions
The package provides a set of utilities used to curate and export the data provided with the package and some functions used for analysis of the data in published studies.
The most important function is `representative.alignment`, which helps to export phylogenetically representative (at the level of intraspecific clades, species, or genera) or complete (all individuals) datasets, exporting a concantenate of selected loci or complete mitogenome sequences. As the data provided in the package are aligned, the exported dataset can be radily used for analyses.

### Do you want to contribute?

### Funding



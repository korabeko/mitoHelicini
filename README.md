# mitoHelicini

## An R package containing a collection of geo-referrenced sequences of complete mitogenomes and partial mitochondrial genes of the land snail family Helicidae, tribe Helicini, and fuctions for curating and exporting the data.
maintained by Ondřej Korábek (Department of Zoology, Faculty of Science, Charles University, Praha, Czechia, ondrej.korabek@natur.cuni.cz)

### Installation
(to be added...)

### For the impatient
To export an alignment ready for analysis: run `help(representative.alignment)` and then try the function with different settings. To see the whole data, go to the /inst/extdata folder or use `export.Helicini()` to copy its contents into your working directory.

### The data
#### What are Helicini?
The tribe [Helicini](https://en.wikipedia.org/wiki/Helicini) Rafinesque, 1815 is a clade of land snails from the Western Palearctic (Gastropoda: Heterobranchia: Pneumopulmonata: Stylommatophora: Helicidae: Helicini). It comprises 10 genera of large snails (shell diameter ~2-6 cm). The most diverse genus is [*Helix*](https://en.wikipedia.org/wiki/Helix_(gastropod)) Linnaeus, 1758, which includes [*Helix pomatia*](https://en.wikipedia.org/wiki/Helix_pomatia) Linnaeus, 1758, one of the best studied land snail species.
#### What data are included?
The dataset is an updated and expanded version of that published by Korábek et al. 2022, *Diversity* 14: 24 [https://doi.org/10.3390/d14010024](https://doi.org/10.3390/d14010024). The original data from that paper are available from [Dryad](https://doi.org/10.5061/dryad.pnvx0k6p5). The data consist of nucleotide sequences of partial mitochondrial genes obtained by PCR and assembled sequences of complete and nearly complete mitogenomes. These are data either obtained by O. Korábek or downloaded from public databases (mostly [NCBI's Nucleotide database/GenBank](https://www.ncbi.nlm.nih.gov/nuccore), but also [BOLD](https://v4.boldsystems.org/) and [GBOL](https://data.bolgermany.de/ergebnisse/results). Only sequences accompanied by reasonably precise geographic location are included. Most of the data come from PCR targeting selected parts of the mitochondrial genome and followed by Sanger sequencing. For many individuals, one or more fragments of the mitogenome are thus available with gaps of unknown length in between. Complete mitogenome sequences, mostly assembled from Illumina reads (WGS or RNAseq) are available for a small (but increasing) proportion of the sequenced individuals.

The sequences are provided with complete metadata (locality, date of collection, shell vouchers, the person responsible for identification, etc.). Furthermore, an estimate of the phylogenetic tree of the included samples is provided, although not for all the individuals.
#### Do you want to contribute a sequence?
In the future, the should be a function for that. Until that materializes, write me an email.
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
(to be added...) For a rough idea, see [Korábek et al. 2022](https://doi.org/10.3390/d14010024).

### The metadata
(to be added...) Sometimes the package offers more precise metadata than available from GenBank or the original publication, provided from authors upon my request.

**ID** unique identifier of a sequenced individual

**clade**	“species-level” clade sensu Korábek et al. (2022a)

**subclade** intraspecific clade sensu Korábek et al. (2022a)

**genus**	genus

**species**	species; slash indicates intermediate or hybrid individuals

**identified_by** author of the species-level identification, only the latest revision is given ("authors" indicates that identification comes from the study/database cited as source of the data)

**date_identified** year of identification to species, only the latest revision is given

**basis_for_identification** observations considered for species identification (shell=conchological characters, mtDNA=membership to a mitochondrial lineage, genital=characters on the genital system, origin=the geographic provenance of the sample, SNPs=multilocus data (for methods see the respective publication, for *H. thessalica*/*pomatia* 2024 study 90% cluster membership from STRUCTURE at K=2 was used))

**locality** description of the location of the sampling site

**latitude** northern latitude in degrees, WGS84 datum (negative values south of equator); number of decimal positions roughly reflects the precision of the geographic coordinates

**longitude**	eastern longitude in degrees, WGS84 datum (negative values south of equator); number of decimal positions roughly reflects the precision of the geographic coordinates

**native** does the sample come from within the natural range of the species (Y=yes, N=No, U=unknown/uncertain; see Korábek et al. 2022a for details)

**collected_by** who collected the sample

**date_of_collection** when was the sample collected in the field

**shells_deposited_in**	acronym of the institution or collection, where the shells are deposited (include any shells of the same species collected at the same site on the same occasion, not only the sequenced individual); collection acronyms are listed on a separate sheet “collection_acronyms”; please, contact me if you need information regarding the tissue/DNA vouchers.

**collection_number**	catalogue number, under which the shells are registered in the collection

**isolate_code** isolate codes (lab codes used to identify the individual and DNA extracts made from it); see "comments" to check for additional codes for the same specimen

**16S_acc_nos** GenBank accession for the *rrnL* rRNA gene

**COX1_acc_nos** GenBank accession for the *cox1* gene

**COX2_acc_nos** GenBank accession for the *cox2* gene

**12S_acc_nos**	GenBank accession for the *rrnS+trnM* gene

**CYTB_acc_nos** GenBank accession for the *cytb* gene

**other_genes** semicolon-separated GenBank accessions for other mitochondrial genes/regions (specified in [] as cox1, rrnL, trnS2, etc.)

**genome_acc_nos** GenBank accession for a complete mitogenome

**SRA**	semicolon-separated Sequence Read Archive accessions with character of the data stated (specified in [] as ddRAD, WGS, RNAseq, etc.)

**reference**	the study where the sequences were first reported (full references on a separate sheet “publications”)

**comments**	various comments on nomenclatoric significance of the individual, accuracy/reliability of the geographic coordinates, source of the material, species identification, or sequence reliability separated by semicolon; “GPS only indicative” at the beginning means that coordinates were determined from locality description only (typically by someone else than the collector) or that the precision is uncertain and the sampling point may be in a broader area around that point (used inconsistently, probably too conservative for material from ZMH and HNHM collected after 2010; in the former used for material where source of coordinates was not stated as directly measured, in the latter I lost track how where they obtained for which sample).

### The files
Only name stem is given here, which is preceded by version in the actual file names.
#### alignment files (plain text in FASTA format):
**Helicini_COX1**: *cox1*, mostly incomplete sequences.

**Helicini_TRNV**: complete region between *cox1* and *rrnL*; only for individuals with *cox1* 3' complete and *rrnL* 5' complete.

**Helicini_16S**: *rrnL*; because the exact start of the gene us uncertain, all nucleotides after *trnV* are included.

**Helicini_TRNL1_TRNA**: complete region between *rrnL* and *nd6* alignments, only for individuals with *rrnL* 3' complete and *nd6* 5' complete.

**Helicini_ND6**: *nd6*; gene ends with a T preceding the polyA tail even where this T follows after a complete stop codon.

**Helicini_TRNP**: complete region between *nd6* and *nd5* alignments, only for individuals with *nd6* 3' complete and *nd5* 5' complete

**Helicini_ND5**: *nd5*; ends always with a full stop codon, because *nd1* and *nd4L* are also translated from the same transcript.

**Helicini_ND1_ND4L**: *nd1* and *nd4L*; *nd1* is incomplete at the 5' end due to overlap with the preceding *nd5*; the start codon of *nd4L* usually overlaps the stop codon of *nd1*; would include also any nucleotides between *nd4L* and *cytb*.

**Helicini_CYTB**: *cytb*, mostly incomplete.

**Helicini_TRND-TRNF**: complete region between *cytb* and *cox2* alignments, only for individuals with *cytb* 3' complete and *cox2* 5' complete.

**Helicini_COX2**: COX2, mostly incomplete.

**Helicini_TRNY-TRNL2**: complete region between *cox2* and *atp8* alignments, only for individuals with *cox2* 3' complete and *atp8* 5' complete.

**Helicini_ATP8**: *atp8*.

**Helicini_TRNN** complete region between *atp8* and *atp6* alignments, only for individuals with *atp8* 3' complete and *atp6* 5' complete.

**Helicini_ATP6**: *atp6*.

**Helicini_TRNR_TRNE**: complete region between *atp6* and *rrnS* alignments, only for individuals with *atp6* 3' complete and *rrnS* alignment 5' complete.

**Helicini_12S_TRNM**: *rrnS* and *trnM*, mostly incomplete; the end position of the *rrnS* gene (5' end of the sequence here) is uncertain and therefore the TRNR_TRNE alignment file possibly includes part of *rrnS*; because the *rrnS* amplicon used in PCR-based studies extends into *trnM*, *trnM* is kept either up to the primer or, if complete, up to the 5' end of *nd3*.

**Helicini_ND3**: *nd3*.

**Helicini_TRNS2_TRNT**]: complete region between *cox3* and *nd3* alignments, only for individuals with *cox3* 3' complete and *nd3* 5' complete.

**Helicini_COX3**: *cox3*; only the T (A at the beginning the sequence, because it is in reverse complement as the gene is on the minus strand) preceeding the polyA tail is homologous between species, even though in some species this is followed by AA, seemingly forming a full stop codon.

**Helicini_TRNS1**: complete region between *nd3* and *nd4* alignments, only for individuals with *nd3* 3' complete and *nd4* 5' complete; probably includes the control region in the 5' part.

**Helicini_ND4**: *nd4*; it covers the whole CDS and all the nucleotides preceeding the polyA tail.

**Helicini_TRNI**: complete region between *nd4* and *nd2* alignments, only for individuals with *nd4* 3' complete and *nd2* 5' complete.

**Helicini_ND2**: *nd2*.

**Helicini_TRNK**: complete region between *nd2* and *cox1* alignments, only for individuals with *nd2* 3' complete and *cox1* 5' complete.

#### sequence metadata (Open Document spread sheet)
**Helicini_table_mol**: species identifications, sampling locality data, NCBI accession numbers, shell vouchers, and other relevant metadata relating to the origin of the sequences and their use in publications; detailed explanation of columns inside the file (and above under **"The metadata"**)

#### mitogenome annotations (plain text file)
**mitogenome_annotations**: annotation of mitogenomic sequences in the GFF3 format, including the positions where the the individual alignment files align. The current annotation is indicated by "current" in the "source" field.

#### change log (plain text file)
**change_log**: record of revisions of the above files, in particular of the alignment files.


### The functions
The package provides a set of utilities used to curate and export the data provided with the package and some functions used for analysis of the data in published studies.
The most important function is `representative.alignment`, which helps to export phylogenetically representative (at the level of intraspecific clades, species, or genera) or complete (all individuals) datasets, exporting a concantenate of selected loci or complete mitogenome sequences. As the data provided in the package are aligned, the exported dataset can be radily used for analyses.

### Funding
The core of the R package was developed during my postdoc funded by scholarship from the Alexander-von-Humboldt Stiftung.

### References (sources of the sequence data)
Bouaziz–Yahiatene H, Pfarrer B, Medjdoub-Bensaad F, Neubert E. (2017). Revision of *Massylaea* Möllendorff, 1898 (Stylommatophora, Helicidae). *ZooKeys* 694: 109–133. https://doi.org/10.3897/zookeys.694.15001

Cadahía L, Harl J, Duda M, Sattmann H, Kruckenhauser L, Fehér Z, Zopp L, Haring E. (2014). New data on the phylogeny of Ariantinae (Pulmonata, Helicidae) and the systematic position of *Cylindrus obtusus* based on nuclear and mitochondrial DNA marker sequences. *Journal of Zoological Systematics and Evolutionary Research* 52: 163–169. https://doi.org/10.1111/jzs.12044

Cesaroni D, De Felici S, Riccarducci G, Ciambotta M, Ventura A, Bianchi E, Sbordoni V. (2017). DNA Barcodes of the animal species occurring in Italy under the European “Habitats Directive” (92/43/EEC): a reference library for the Italian National Biodiversity Network. *Biogeographia – The Journal of Integrative Biogeography* 32: 5–23. https://doi.org/10.21426/B632131365

Dahirel M, Olivier E, Guiller A, Martin M–C, Madec L, Ansart A. (2015). Movement propensity and ability correlate with ecological specialization in European land snails: comparative analysis of a dispersal syndrome. *Journal of Animal Ecology* 84: 228–238. https://doi.org/10.1111/1365-2656.12276

Dimzas D, Morelli S, Traversa D, Di Cesare A, Van Bourgonie YR, Breugelmans K, Backeljau T, Frangipane di Regalbono A, Diakou A. (2020). Intermediate gastropod hosts of major feline cardiopulmonary nematodes in an area of wildcat and domestic cat sympatry in Greece. *Parasites & Vectors* 13: 345. https://doi.org/10.1186/s13071-020-04213-z

Fiorentino V, Manganelli G, Giusti F, Ketmaier V. (2016). Recent expansion and relic survival: Phylogeography of the land snail genus *Helix* (Mollusca, Gastropoda) from south to north Europe. *Molecular Phylogenetics and Evolution* 98: 358–372. 
https://doi.org/10.1016/j.ympev.2016.02.017

Gittenberger E, Piel WH, Groenenberg DSJ. (2004). The Pleistocene glaciations and the evolutionary history of the polytypic snail species *Arianta arbustorum* (Gastropoda, Pulmonata, Helicidae). *Molecular Phylogenetics and Evolution* 30: 64–73. https://doi.org/10.1016/S1055-7903(03)00182-9

Giusti F, Fiorentino V, Manganelli G. (2015). A neotype for *Helix cincta* Müller, 1774 (Gastropoda, Pulmonata, Helicidae). *Journal of Conchology* 42: 209–212. https://conchsoc.org/sites/default/files/jconch/42/2/2015-42209.pdf

Groenenberg DSJ, Duijm E. (2019). The complete mitogenome of the Roman snail *Helix pomatia* Linnaeus 1758 (Stylommatophora: Helicidae). *Mitochondrial DNA Part B* 4: 1494–1495. https://doi.org/10.1080/23802359.2019.1601512

Hausdorf B, Parr M, Shappell LJ, Oldeland J, Robinson DG. (2021). The introduction of the European *Caucasotachea vindobonensis* (Gastropoda: Helicidae) in North America, its origin and its potential range. *Biological Invasions* 23: 3281–3289. https://doi.org/10.1007/s10530-021-02579-4

Jaksch K, Eschner A, von Rintelen T, Haring E. (2016). DNA analysis of molluscs from a museum wet collection: a comparison of different extraction methods. *BMC Research Notes* 9: 348. https://doi.org/10.1186/s13104-016-2147-7

Kajtoch Ł, Davison A, Grindon A, Deli T, Sramkó G, Gwardjan M, Kramarenko S, Mierzwa-Szymkowiak D, Ruta R, Tóth JP, Wade C, Kolasa M, Egorov RV, Fehér Z. (2017). Reconstructed historical distribution and phylogeography unravels non-steppic origin of *Caucasotachea vindobonensis* (Gastropoda: Helicidae). *Organisms Diversity & Evolution* 17: 679–692. https://doi.org/10.1007/s13127-017-0337-3

Ketmaier V, Glaubrecht M. (2015). The legacy of the Crusaders: Complex history of colonization and anthropochory in the land snails *Levantina* (Gastropoda, Pulmonata) in the Eastern Mediterranean. *Zoosystematics and Evolution* 91: 81–89. https://doi.org/10.3897/zse.91.4693

Korábek O, Juřičková L, Petrusek A. (2014). Resurrecting *Helix straminea*, a forgotten escargot with trans‐Adriatic distribution: first insights into the genetic variation within the genus Helix (Gastropoda: Pulmonata). *Zoological Journal of the Linnean Society* 171: 72–91. https://doi.org/10.1111/zoj12122

Korábek O, Petrusek A, Neubert E, Juřičková L. (2015). Molecular phylogeny of the genus *Helix* (Pulmonata: Helicidae). *Zoologica Scripta* 44: 263–280. https://doi.org/10.1111/zsc.12101

Korábek O, Juřičková L, Petrusek A. (2016). Splitting the Roman snail *Helix pomatia* Linnaeus, 1758 (Stylommatophora: Helicidae) into two: redescription of the forgotten *Helix thessalica* Boettger, 1886. *Journal of Molluscan Studies* 82: 11–22. https://doi.org/10.1093/mollus/eyv048

Korábek O, Petrusek A, Juřičková L. (2018). Glacial refugia and postglacial spread of an iconic large European land snail, *Helix pomatia* (Pulmonata: Helicidae). *Biological Journal of the Linnean Society* 123: 218–234. https://doi.org/10.1093/biolinnean/blx135

Korábek O, Juřičková L, Balashov I, Petrusek A. (2018). The contribution of ancient and modern anthropogenic introductions to the colonization of Europe by the land snail *Helix lucorum* Linnaeus, 1758 (Helicidae). *Contributions to Zoology* 87: 61–74. https://doi.org/10.1163/18759866-08702001

Korábek O, Petrusek A, Rovatsos M. (2019). Mitogenome of *Helix pomatia* and the basal phylogeny of Helicinae (Gastropoda, Stylommatophora, Helicidae). *ZooKeys* 827: 19–30. https://doi.org/10.3897/zookeys.827.33057

Korábek O, Juřičková L, Petrusek A. (2020). Inferring the sources of postglacial range expansion in two large European land snails. *Journal of Zoological Systematics and Evolutionary Research* 58: 944–956. https://doi.org/10.1111/jzs.12368

Korábek O, Kosová T, Dolejš P, Petrusek A, Neubert E, Juřičková L. (2021). Geographic isolation and human–assisted dispersal in land snails: a Mediterranean story of *Helix borealis* and its relatives (Gastropoda: Stylommatophora: Helicidae). *Zoological Journal of the Linnean Society* 193: 1310–1335. https://doi.org/10.1093/zoolinnean/zlaa186

Korábek O, Juřičková L, Petrusek A. (2022a). Diversity of land snail tribe Helicini (Gastropoda: Stylommatophora: Helicidae): Where do we stand after 20 years of sequencing mitochondrial markers? *Diversity* 14: 24. https://doi.org/10.3390/d14010024

Korábek O, Glaubrecht M, Hausdorf B, Neiber MT. (2022b). Phylogeny of the land snail *Levantina* reveals long-distance dispersal in the Middle East. *Zoologica Scripta* 51: 161–172. https://doi.org/10.1111/zsc.12526

Korábek O, Adamcová T, Proćków M, Petrusek A, Hausdorf B, Juřičková L. (2023a). In both directions: Expansions of European land snails to the north and south from glacial refugia. *Journal of Biogeography* 40: 654–668. https://doi.org/10.1111/jbi.14531

Korábek O, Balashov I, Neiber MT, Walther F, Hausdorf B. (2023b). The Caucasus is neither a cradle nor a museum of diversity of the land snail genus *Helix* (Gastropoda, Stylommatophora, Helicidae), while Crimea is home to an ancient lineage. *Zoosystematics and Evolution* 99: 535–543. https://doi.org/10.3897/zse.99.110610

Korábek O, Hausdorf B. (2023c). Unravelling the double *Helix escherichi*, *with description of *Helix ankae* sp. nov. and a discussion of Maltzanella and the subgenera of Helix (Gastropoda: Helicidae). *Archiv für Molluskenkunde* 152: 239–256. https://doi.org/10.1127/arch.moll/152/239-256

Korábek O, Hausdorf B. (2024). Accelerated mitochondrial evolution and asymmetric viability of hybrids contribute to the persistence of *Helix thessalica* in the *Helix pomatia* range. *Molecular Ecology* 33: e17474. https://doi.org/10.1111/mec.17474

Kotsakiozi P, Parmakelis A, Giokas S, Papanikolaou I, Valakos ED. (2012). Mitochondrial phylogeny and biogeographic history of the Greek endemic land-snail genus *Codringtonia* Kobelt 1898 (Gastropoda, Pulmonata, Helicidae). *Molecular Phylogenetics and Evolution* 62: 681–692. https://doi.org/10.1016/j.ympev.2011.11.012

Layton KKS, Warne CPK, Nicolai A, Ansart A, deWaard JR. (2019). Molecular evidence for multiple introductions of the banded grove snail (*Cepaea nemoralis*) in North America. *Canadian Journal of Zoology* 97: 392–398. https://doi.org/10.1139/cjz-2018-0084

Manganelli G, Salomone N, Giusti F. (2005). A molecular approach to the phylogenetic relationships of the western palaearctic Helicoidea (Gastropoda: Stylommatophora). *Biological Journal of the Linnean Society* 85: 501–512. https://doi.org/10.1111/j.1095-8312.2005.00514.x

Mumladze L, Tarkhnishvili D, Murtskhvaladze M. (2013). Systematics and evolutionary history of large endemic snails from the Caucasus (*Helix buchii* and *H. goderdziana*) (Helicidae). *American Malacological Bulletin* 31: 225–234. https://doi.org/10.4003/006.031.0202

Neiber MT, Hausdorf B. (2015). Molecular phylogeny reveals the polyphyly of the snail genus *Cepaea* (Gastropoda: Helicidae). *Molecular Phylogenetics and Evolution* 93: 143–149. https://doi.org/10.1016/j.ympev.2015.07.022

Neiber MT, Sagorny C, Hausdorf B. (2016a). Increasing the number of molecular markers resolves the phylogenetic relationship of *‘Cepaea’ vindobonensis* (Pfeiffer 1828) with *Caucasotachea* Boettger 1909 (Gastropoda: Pulmonata: Helicidae). *Journal of Zoological Systematics and Evolutionary Research* 54: 40–45. https://doi.org/10.1111/jzs.12116

Neiber MT, Sagorny C, Sauer J, Walther F, Hausdorf B. (2016b). Phylogeographic analyses reveal Transpontic long distance dispersal in land snails belonging to the *Caucasotachea atrolabiata* complex (Gastropoda: Helicidae). *Molecular Phylogenetics and Evolution* 103: 172–183. https://doi.org/10.1016/j.ympev.2016.07.017

Neiber MT, Korábek O, Glaubrecht M, Hausdorf B. (2022). A misinterpreted disjunction: the phylogenetic relationships of the Libyan land snail *Gyrostomella* (Gastropoda: Stylommatophora: Helicidae). *Zoological Journal of the Linnean Society* 194: 1236–1251. https://doi.org/10.1093/zoolinnean/zlab059

Parmakelis A, Kotsakiozi P, Rand D. (2013). Animal mitochondria, positive selection and cyto-nuclear coevolution: insights from Pulmonates. *PLoS ONE* 8: e61970. https://doi.org/10.1371/journal.pone.0061970

Petraccioli A, Niero A, Carandente F, Crovato P, de Vico G, Odierna G, Picariello OLA, Tardy E, Viglietti S, Guarino FM, Maio N. (2021). *Helix straminea* Briganti, 1825 in Italy (Gastropoda: Pulmonata): taxonomic history, morphology, biology, distribution and phylogeny. *The European Zoological Journal* 88: 390–416. https://doi.org/10.1080/24750263.2021.1892217

Petraccioli A, Crovato P, Guarino FM, Mezzasalma M, Odierna G, Picariello O, Maio N. (2021). Chromosome diversity and evolution in Helicoidea (Gastropoda: Stylommatophora): A synthesis from original and literature data. *Animals* 11: 2551. https://doi.org/10.3390/ani11092551

Psonis N, Vardinoyannis K, Mylonas M, Poulakakis N. (2015). Evaluation of the taxonomy of *Helix cincta* (Muller, 1774) and *Helix nucula* (Mousson, 1854); insights using mitochondrial DNA sequence data. *Journal of Natural History* 49: 383–392. https://doi.org/10.1080/00222933.2013.825023

Psonis N, Vardinoyannis K, Poulakakis N. (2022). High-throughput degraded DNA sequencing of subfossil shells of a critically endangered stenoendemic land snail in the Aegean. *Molecular Phylogenetics and Evolution* 175: 107561. https://doi.org/10.1016/j.ympev.2022.107561

Razkin O, Gómez-Moliner BJ, Prieto CE, Martínez-Ortí A, Arrébola JR, Muñoz B, Chueca LJ, Madeira MJ. (2015). Molecular phylogeny of the western Palaearctic Helicoidea (Gastropoda, Stylommatophora). *Molecular Phylogenetics and Evolution* 83: 99–117. https://doi.org/10.1016/j.ympev.2014.11.014

Schmera D, Pizá J, Reinartz E, Ursenbacher S, Baur B. (2016). Breeding system, shell size and age at sexual maturity affect sperm length in stylommatophoran gastropods. *BMC Evolutionary Biology* 16: 89. https://doi.org/10.1186/s12862-016-0661-9

Sei M, Robinson DG, Geneva AJ, Rosenberg G. (2017). Doubled helix: Sagdoidea is the overlooked sister group of Helicoidea (Mollusca: Gastropoda: Pulmonata). *Biological Journal of the Linnean Society* 122: 697–728.	https://doi.org/10.1093/biolinnean/blx082

Zając KS, Turek J, Boroń A. (2023). *Helix lucorum* Linnaeus, 1758 (Gastropoda: Helicidae) – the morphological and molecular analysis of a new species to the Polish malacofauna. *Folia Biologica (Kraków)* 71: 137–145.	https://doi.org/10.3409/fb_71-3.14

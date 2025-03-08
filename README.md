# mitoHelicini

## An R package containing a collection of geo-referrenced sequences of complete mitogenomes and partial mitochondrial genes of the land snail family Helicidae, tribe Helicini, and fuctions for curating and exporting the data.
maintained by Ondřej Korábek (Department of Zoology, Faculty of Science, Charles University, Praha, Czechia)

### The data
#### What are Helicini?
The tribe [Helicini](https://en.wikipedia.org/wiki/Helicini) is a clade of land snails from the Western Palearctic (Gastropoda: Heterobranchia: Pneumopulmonata: Stylommatophora: Helicidae: Helicini). It comprises 10 genera of large snails (shell diameter ~2-6 cm). The most diverse genus is [*Helix*](https://en.wikipedia.org/wiki/Helix_(gastropod)), which includes [*Helix pomatia*](https://en.wikipedia.org/wiki/Helix_pomatia) Linnaeus, 1758, one of the best studied land snail species.
#### What data are included?
The dataset is an updated and expanded version of that published by Korábek et al. 2022, *Diversity* 14: 24 [https://doi.org/10.3390/d14010024](https://doi.org/10.3390/d14010024). The original data from that paper are available from [Dryad](https://doi.org/10.5061/dryad.pnvx0k6p5).


### The functions
The package provides a set of utilities used to curate and export the data provided with the package and some functions used for analysis of the data in published studies.
The most important function is `representative.alignment`, which helps to export phylogenetically representative (at the level of intraspecific clades, species, or genera) or complete (all individuals) datasets, exporting a concantenate of selected loci or complete mitogenome sequences. As the data provided in the package are aligned, the exported dataset can be radily used for analyses.




## For the Human 16S rRNA Gene Sequencing
The Required Input Data
* Microbe Abundance (no-miss-table-dada2.qza):An abundance table with microbial species in rows and samples in columns.
* Taxonomy Table (taxonomy.qza): Table of the taxonomic annotation
* MetaData (BX_Sputum_convert.txt): A spreadsheet with samples in rows and metadata in columns
* Tree File (rooted-tree_quality.qza): 

1. Open R
1. Follow Code Bronchiectasis.Neutrophil.Elastase.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```

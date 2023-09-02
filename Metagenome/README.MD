*Metagenome Analysis*
## Running the Analysis
The Required Input Data
* Microbe Raw Abundance (counts_taxa_metagenome_cat.txt):An abundance table with microbial species in rows and samples in columns.
* MetaData (COPD_MAP_Final_git.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COPD.edgeR.Function.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```

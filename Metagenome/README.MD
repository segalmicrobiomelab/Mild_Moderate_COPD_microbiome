*Metagenome Analysis*
## Running the Analysis
The Required Input Data
* Microbe Raw Abundance (taxonomy_bacteria_metatranscriptome_G.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_bacteria_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201104.WGS.UPDATE.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19_Functional_Trial.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```

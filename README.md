## Integrative analysis reveals mouse strain-dependent responses to acute ozone exposure associated with airway macrophage transcriptional activity

This repository contains code for analyses and figures contained in this manuscript, which were most recently tested with R 4.1.1 on September 26, 2021.

The analysis is separated into several files, all sourcing two files: `0_functions.R`, which contains the environment and accessory functions and `1_data_cleaning.R`, which prepares raw data. The subsequent files should be run separately, and you should remember to set your working directory back to the analysis directory after running each with `setwd()`. A separate folder contains a raw data spreadsheet which is also available on the manuscript's figshare repository. A script within `1_data_cleaning.R` will download the gene expression counts from GEO.

This analysis is open source under an MIT license.
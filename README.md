## Integrative analysis reveals mouse strain-dependent responses to acute ozone exposure associated with airway macrophage transcriptional activity

This repository contains code for analyses and figures contained in this manuscript, which were most recently tested with R 4.1.1 on September 26, 2021.

The analysis is separated into several files, all sourcing two scripts: `0_functions.R`, which contains the environment and accessory functions and `1_data_cleaning.R`, which prepares raw data. The subsequent scripts should be run in order, and you should remember to set your working directory back to the analysis directory after running each with `setwd()`. A separate folder contains a raw data spreadsheet which is also available on the manuscript's figshare repository. A function within `1_data_cleaning.R` will download the gene expression counts from GEO.

Code was written by [Wes Crouse](https://github.com/wesleycrouse) and [Adelaide Tovar](https://github.com/adelaidetovar), with specific contributions noted in each script.

This analysis is open source under an MIT license.

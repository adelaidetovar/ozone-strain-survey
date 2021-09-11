source("0_functions.R")

setwd("../../data")

####################
#### Data Input ####
####################

raw_data <- read.xlsx('raw_data.xlsx', sheet = 2, startRow = 2)
raw_data[,c(2:4)] <- sapply(raw_data[,c(2:4)], as.factor)
raw_data$strain <- factor(raw_data$strain, levels = c("C57BL/6J", "CC003", "CC017", "CC025", "CC039", "CC059"))

# BAL
bal <- raw_data[,c(1:14)]
bal[,c(6:14)] <- sapply(bal[,c(6:14)], as.numeric)

# protein
protein <- raw_data[,c(1:5, 15)]
protein$protein_conc <- as.numeric(protein$protein_conc)

# albumin
albumin <- raw_data[,c(1:5, 16)]
# remove samples that weren't measured
albumin <- albumin[complete.cases(albumin),]
# not detected samples become NA
albumin$albumin_conc <- as.numeric(albumin$albumin_conc)

# cytokines
cytokines <- raw_data[,c(1:5, 17:31)]
# remove samples that weren't measured
cytokines <- cytokines[complete.cases(cytokines),]
# not detected samples become NA
cytokines[,c(6:20)] <- sapply(cytokines[,c(6:20)], as.numeric)


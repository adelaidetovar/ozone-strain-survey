source("0_functions.R")

setwd("../data")

##############################
#### Phenotype Data Input ####
##############################

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

####################################
#### Gene Expression Data Input ####
####################################

#load subject information
info <- data.frame(read_xlsx("rnaseq_ids.xlsx"))

rownames(info) <- paste0("X", info$mouse_no)
info$strain[info$strain=="C57BL/6J"] <- "B6"

# Download expression data from GEO
getGEOSuppFiles("GSE174205", filter_regex = ".txt.gz")
df <- read.table(gzfile("../GSE174205/GSE174205_summary_count_matrix.txt.gz"), head=T, row.names=NULL)
df <- df[!duplicated(df$Genes),]
rownames(df) <- NULL
df <- df %>% column_to_rownames(var = "Genes")

df <- t(df)
rownames(df) <- sapply(rownames(df), function(x){unlist(strsplit(x, split="_"))[1]})

# Make phenotype data for Wes's gene expression analyses

ids <- rownames(info)
ids <- gsub("X", "", ids)
set <- c("mouse_no","eotaxin", "gcsf", "gmcsf", "il10", "il12p70", "il6",
         "ip10", "kc", "lix", "mcp1", "mip1a", "mip1b", "mip2")

pheno <- full_join(cytokines[cytokines$mouse_no%in%ids,set],
               bal[bal$mouse_no%in%ids,c("mouse_no","per_neu", "no_neu", "per_macs", "no_macs")],
               by = "mouse_no")
pheno <- full_join(pheno,
                   protein[protein$mouse_no%in%ids, c("mouse_no","protein_conc")],
                   by = "mouse_no")
pheno$mouse_no <- paste0("X",pheno$mouse_no)
pheno <- pheno %>% column_to_rownames("mouse_no")

####################
#data processing and formatting

#drop subjects without pairs
info <- info[info$pair!="x",]
df <- df[rownames(info),]
pheno <- pheno[rownames(info),]

#drop genes with fewer counts than samples
df <- df[,colSums(df) >= nrow(df)]

#ensure genes are integer counts
mode(df) <- "integer"

#drop genes with identical counts
id <- apply(df, 2, paste, collapse=",")
id_table <- table(id)
id_table <- id_table[id_table!=1]

df <- df[,!(colnames(df) %in% unlist(lapply(names(id_table), function(x){names(which(id==x))})))]
save(df, file="df.RData")

#create unique pair ID
info$pair_ID <- apply(info, 1, function(x){paste(x["strain"], x["pair"], sep="_") })
info$pair_ID <- as.factor(info$pair_ID)

#format variables
info <- info[,-(1:2)]
info$strain <- as.factor(info$strain)
info$rx <- as.factor(info$rx)
info$sex <- as.factor(info$sex)
info$pair <- as.factor(info$pair)

info <- info[,!(c("mouse_no", "file") %in% rownames(info))]
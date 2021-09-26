source("0_functions.R")
source("1_data_cleaning.R")

output_dir <- file.path("../output/gene_expr/")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

setwd(output_dir)
rm(output_dir)

##########################################
#### Differential Expression Analysis ####
##########################################

#initialize object for paired samples with G, T effect
dds <- DESeqDataSetFromMatrix(countData = t(df),
                              colData = info,
                              design= ~ pair_ID + rx)

#design matrix for paired samples with G, T effects; cell-means for strain and sum-to-zero contrasts for strain:pair
#this makes the G effect average over pairs
X_gt <- model.matrix(~0 + strain, info)

pair_contrast <- model.matrix(~0 + pair:strain, info)
pair_contrast <- pair_contrast[,colSums(pair_contrast)!=0]

for (i in levels(info$strain)){
  pair_contrast_current <- pair_contrast[,grep(i, colnames(pair_contrast))]
  pair_contrast_current <- pair_contrast_current%*%sumtozero_contrast(ncol(pair_contrast_current))
  colnames(pair_contrast_current) <- sapply(1:ncol(pair_contrast_current), function(x){paste("pair", i, x, sep="_")})

  X_gt <- cbind(X_gt, pair_contrast_current)
}

rm(pair_contrast, pair_contrast_current)

X_gt <- cbind(X_gt, model.matrix(~rx, info)[,2,drop=F])

#call DESeq for G,T vs G
dds <- DESeq(dds, full=X_gt)
save(dds, file="dds_Wald_gt_g.RData")

#shrunken log fold change for T using ashr
resLFC <- lfcShrink(dds, coef="rxO3", type="ashr")
save(resLFC, file="resLFC_LRT_gt_g.RData")

#output list of differentially expressed genes
table <- data.frame(resLFC)

write.table(rownames(table)[!is.na(table$padj) & table$padj < 0.05 & table$log2FoldChange>0], file="DE_positive.txt", quote=F, row.names=F, col.names=F)
write.table(rownames(table)[!is.na(table$padj) & table$padj < 0.05 & table$log2FoldChange<0], file="DE_negative.txt", quote=F, row.names=F, col.names=F)
write.table(rownames(table)[!is.na(table$padj) & table$padj < 0.05], file="DE_all.txt", quote=F, row.names=F, col.names=F)

#store output table
output <- resLFC

#---------------

#collect G coefficients for each strain
G_coef <- matrix(NA, ncol(df), nlevels(info$strain))

colnames(G_coef) <- levels(info$strain)
rownames(G_coef) <- colnames(df)

G_pval <- G_coef
colnames(G_pval) <- sapply(levels(info$strain), paste0, "_pval")

G_SE <- G_coef
colnames(G_SE) <- sapply(colnames(G_pval), paste0, "_SE")

# #collect G coefficients for each strain
for (i in 1:ncol(G_coef)){
  c_vector <- rep(-1/(ncol(G_coef)-1), ncol(G_coef))
  c_vector[i] <- 1
  c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-ncol(G_coef)))

  dds_current <- results(dds, c_vector)

  G_coef[,i] <- dds_current$log2FoldChange
  G_pval[,i] <- dds_current$pvalue
  G_SE[,i] <- dds_current$lfcSE

#   #assign NA for G_pval if NA for DESeq2 padj
  G_pval[is.na(dds_current$padj),i] <- NA
}

G_coef_vec <- as.vector(G_coef)
G_SE_vec <- as.vector(G_SE)

# #jointly shrink independent G effects
G_coef_ash <- ash(G_coef_vec, G_SE_vec, mixcompdist = "normal", method = "shrink")

G_coef_ash <- matrix(get_pm(G_coef_ash), ncol=ncol(G_coef))
colnames(G_coef_ash) <- colnames(G_coef)
rownames(G_coef_ash) <- rownames(G_coef)

G_coef <- G_coef_ash
rm(G_coef_ash)

# #FDR correction for G p-values
G_padj <- matrix(p.adjust(as.vector(G_pval), method="fdr"), ncol=ncol(G_pval))
rownames(G_padj) <- rownames(G_pval)
colnames(G_padj) <- sapply(colnames(G_coef), function(x){paste0(x, "_padj")})

# #G contrasts for specific pairs of strains
c_vector <- rep(0, 6)
c_vector[3] <- 1
c_vector[2] <- -1
c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-ncol(G_coef)))

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

G_coef <- cbind(G_coef, resLFC_current$log2FoldChange)
G_padj <- cbind(G_padj, resLFC_current$padj)

colnames(G_coef)[7] <- "CC017_vs_CC003"
colnames(G_padj)[7] <- "CC017_vs_CC003_padj"

c_vector <- rep(0, 6)
c_vector[5] <- 1
c_vector[2] <- -1
c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-6))

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

G_coef <- cbind(G_coef, resLFC_current$log2FoldChange)
G_padj <- cbind(G_padj, resLFC_current$padj)

colnames(G_coef)[8] <- "CC039_vs_CC003"
colnames(G_padj)[8] <- "CC039_vs_CC003_padj"

c_vector <- rep(0, 6)
c_vector[5] <- 1
c_vector[3] <- -1
c_vector <- c(rep(0,length(resultsNames(dds))-6), c_vector)

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

G_coef <- cbind(G_coef, resLFC_current$log2FoldChange)
G_padj <- cbind(G_padj, resLFC_current$padj)

colnames(G_coef)[9] <- "CC039_vs_CC017"
colnames(G_padj)[9] <- "CC039_vs_CC017_padj"

save(G_coef, file="G_coef.RData")
save(G_pval, file="G_pval.RData")
save(G_padj, file="G_padj.RData")

#CC003 vs CC017 G genes with coefficient
write.csv(as.matrix(G_coef[sapply(G_padj[,"CC017_vs_CC003_padj"]<0.05, isTRUE),"CC017_vs_CC003"]),
          file="CC017_vs_CC003_G_with_coefs.csv", quote=F)

#---------------

#design matrix for paired samples with G, T, X effects; sum-to-zero contrasts
X_gtx <- model.matrix(~0 + strain:rx, info)
X_gtx <- X_gtx[,tail(colnames(X_gtx), nlevels(info$strain))]
X_gtx <- cbind(X_gt, X_gtx)
X_gtx <- X_gtx[,-ncol(X_gt)]

#significant GxT
#call DESeq for G,T,X vs G,T
dds <- DESeq(dds, full=X_gtx, reduced=X_gt, test="LRT")
save(dds, file="dds_LRT_gtx_gt.RData")

#number of significant GxT genes
sum(results(dds)$padj<0.05, na.rm=T)

#save output table
results_df <- results(dds)
output$padj_GTX <- results_df$padj
write.csv(output, file="ozone_DE_GT_GTX.csv")

#VST
df_vst <- t(assay(vst(dds, blind=FALSE)))
save(df_vst, file="ozone_vst.RData")

#call DESeq for G,T,X vs G,T using Wald test
dds <- DESeq(dds, full=X_gtx)
save(dds, file="dds_Wald_gtx.RData")

#---------------

#collect G coefficients for each strain
G_coef_in_GTX <- matrix(NA, ncol(df), nlevels(info$strain))

colnames(G_coef_in_GTX) <- sapply(levels(info$strain), paste0, "_G_in_GTX")
rownames(G_coef_in_GTX) <- colnames(df)
 
G_pval_in_GTX <- G_coef_in_GTX
colnames(G_pval_in_GTX) <- sapply(levels(info$strain), paste0, "_G_in_GTX_pval")
 
G_SE_in_GTX <- G_coef_in_GTX
colnames(G_SE_in_GTX) <- sapply(colnames(G_pval_in_GTX), paste0, "_G_in_GTX_SE")
# 
# # #collect G coefficients for each strain
for (i in 1:ncol(G_coef_in_GTX)){
  c_vector <- rep(-1/(ncol(G_coef_in_GTX)-1), ncol(G_coef_in_GTX))
  c_vector[i] <- 1
  c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-ncol(G_coef_in_GTX)))

  dds_current <- results(dds, c_vector)

  G_coef_in_GTX[,i] <- dds_current$log2FoldChange
  G_pval_in_GTX[,i] <- dds_current$pvalue
  G_SE_in_GTX[,i] <- dds_current$lfcSE

  #assign NA for GTX_pval if NA for DESeq2 padj
  G_pval_in_GTX[is.na(dds_current$padj),i] <- NA
}

G_coef_vec_in_GTX <- as.vector(G_coef_in_GTX)
G_SE_vec_in_GTX <- as.vector(G_SE_in_GTX)

#jointly shrink independent G effects
G_coef_ash_in_GTX <- ash(G_coef_vec_in_GTX, G_SE_vec_in_GTX, mixcompdist = "normal", method = "shrink")

G_coef_ash_in_GTX <- matrix(get_pm(G_coef_ash_in_GTX), ncol=ncol(G_coef_in_GTX))
colnames(G_coef_ash_in_GTX) <- colnames(G_coef_in_GTX)
rownames(G_coef_ash_in_GTX) <- rownames(G_coef_in_GTX)

G_coef_in_GTX <- G_coef_ash_in_GTX
rm(G_coef_ash_in_GTX)

#FDR correction for G p-values
G_padj_in_GTX <- matrix(p.adjust(as.vector(G_pval_in_GTX), method="fdr"), ncol=ncol(G_pval_in_GTX))
rownames(G_padj_in_GTX) <- rownames(G_pval_in_GTX)
colnames(G_padj_in_GTX) <- sapply(colnames(G_coef_in_GTX), function(x){paste0(x, "_padj")})

#G contrasts for specific pairs of strains
c_vector <- rep(0, 6)
c_vector[3] <- 1
c_vector[2] <- -1
c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-6))

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

G_coef_in_GTX <- cbind(G_coef_in_GTX, resLFC_current$log2FoldChange)
G_padj_in_GTX <- cbind(G_padj_in_GTX, resLFC_current$padj)

colnames(G_coef_in_GTX)[7] <- "CC017_vs_CC003_G_in_GTX"
colnames(G_padj_in_GTX)[7] <- "CC017_vs_CC003_G_in_GTX_padj"

c_vector <- rep(0, 6)
c_vector[5] <- 1
c_vector[2] <- -1
c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-6))

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

G_coef_in_GTX <- cbind(G_coef_in_GTX, resLFC_current$log2FoldChange)
G_padj_in_GTX <- cbind(G_padj_in_GTX, resLFC_current$padj)

colnames(G_coef_in_GTX)[8] <- "CC039_vs_CC003_G_in_GTX"
colnames(G_padj_in_GTX)[8] <- "CC039_vs_CC003_G_in_GTX_padj"

c_vector <- rep(0, 6)
c_vector[5] <- 1
c_vector[3] <- -1
c_vector <- c(c_vector, rep(0,length(resultsNames(dds))-6))

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

G_coef_in_GTX <- cbind(G_coef_in_GTX, resLFC_current$log2FoldChange)
G_padj_in_GTX <- cbind(G_padj_in_GTX, resLFC_current$padj)

colnames(G_coef_in_GTX)[9] <- "CC039_vs_CC017_G_in_GTX"
colnames(G_padj_in_GTX)[9] <- "CC039_vs_CC017_G_in_GTX_padj"

save(G_coef_in_GTX, file="G_coef_in_GTX.RData")
save(G_pval_in_GTX, file="G_pval_in_GTX.RData")
save(G_padj_in_GTX, file="G_padj_in_GTX.RData")

sum(apply(G_padj_in_GTX, 1, function(x){any(x<0.05, na.rm=T)}))

#---------------

#collect GxT coefficients for each strain
GTX_coef <- matrix(NA, ncol(df), nlevels(info$strain))

colnames(GTX_coef) <- levels(info$strain)
rownames(GTX_coef) <- colnames(df)

GTX_pval <- GTX_coef
colnames(GTX_pval) <- sapply(levels(info$strain), paste0, "_pval")

GTX_SE <- GTX_coef
colnames(GTX_SE) <- sapply(colnames(GTX_pval), paste0, "_SE")

#collect GxT coefficients for each strain
for (i in 1:ncol(GTX_coef)){
  c_vector <- rep(-1/(ncol(GTX_coef)-1), ncol(GTX_coef))
  c_vector[i] <- 1
  c_vector <- c(rep(0,length(resultsNames(dds))-ncol(GTX_coef)), c_vector)

  dds_current <- results(dds, c_vector)

  GTX_coef[,i] <- dds_current$log2FoldChange
  GTX_pval[,i] <- dds_current$pvalue
  GTX_SE[,i] <- dds_current$lfcSE

  #assign NA for GTX_pval if NA for DESeq2 padj
  GTX_pval[is.na(dds_current$padj),i] <- NA
}

GTX_coef_vec <- as.vector(GTX_coef)
GTX_SE_vec <- as.vector(GTX_SE)

#jointly shrink independent GxT effects
GTX_coef_ash <- ash(GTX_coef_vec, GTX_SE_vec, mixcompdist = "normal", method = "shrink")

GTX_coef_ash <- matrix(get_pm(GTX_coef_ash), ncol=ncol(GTX_coef))
colnames(GTX_coef_ash) <- colnames(GTX_coef)
rownames(GTX_coef_ash) <- rownames(GTX_coef)

GTX_coef <- GTX_coef_ash
rm(GTX_coef_ash)

#FDR correction for GxT p-values
GTX_padj <- matrix(p.adjust(as.vector(GTX_pval), method="fdr"), ncol=ncol(GTX_pval))
rownames(GTX_padj) <- rownames(GTX_pval)
colnames(GTX_padj) <- sapply(colnames(GTX_coef), function(x){paste0(x, "_padj")})

#GxT contrasts for specific pairs of strains
c_vector <- rep(0, 6)
c_vector[3] <- 1
c_vector[2] <- -1
c_vector <- c(rep(0,length(resultsNames(dds))-6), c_vector)

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

GTX_coef <- cbind(GTX_coef, resLFC_current$log2FoldChange)
GTX_padj <- cbind(GTX_padj, resLFC_current$padj)

colnames(GTX_coef)[7] <- "CC017_vs_CC003"
colnames(GTX_padj)[7] <- "CC017_vs_CC003_padj"

c_vector <- rep(0, 6)
c_vector[5] <- 1
c_vector[2] <- -1
c_vector <- c(rep(0,length(resultsNames(dds))-6), c_vector)

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

GTX_coef <- cbind(GTX_coef, resLFC_current$log2FoldChange)
GTX_padj <- cbind(GTX_padj, resLFC_current$padj)

colnames(GTX_coef)[8] <- "CC039_vs_CC003"
colnames(GTX_padj)[8] <- "CC039_vs_CC003_padj"

c_vector <- rep(0, 6)
c_vector[5] <- 1
c_vector[3] <- -1
c_vector <- c(rep(0,length(resultsNames(dds))-6), c_vector)

dds_current <- results(dds, c_vector)
resLFC_current <- lfcShrink(dds, contrast=c_vector, type="ashr")

GTX_coef <- cbind(GTX_coef, resLFC_current$log2FoldChange)
GTX_padj <- cbind(GTX_padj, resLFC_current$padj)

colnames(GTX_coef)[9] <- "CC039_vs_CC017"
colnames(GTX_padj)[9] <- "CC039_vs_CC017_padj"

save(GTX_coef, file="GTX_coef.RData")
save(GTX_pval, file="GTX_pval.RData")
save(GTX_padj, file="GTX_padj.RData")

#number of genes with any significant GTX coefficient; plot baseline expression for these
sum(apply(GTX_padj[,1:6], 1, function(x){any(x<0.05, na.rm=T)}))

#statistics on GTX genes
apply(GTX_padj<0.05, 2, sum, na.rm=T)
apply(GTX_padj<0.05 & GTX_coef>0, 2, sum, na.rm=T)
apply(GTX_padj<0.05 & GTX_coef<0, 2, sum, na.rm=T)

#save gene lists for GTX genes
for (i in 1:ncol(GTX_coef)){
  current <- GTX_padj[,i]
  current <- current[!is.na(current)]
  names(current)[current<0.05]
  
  write.table(names(current)[current<0.05], file=paste0("GTX_", colnames(GTX_coef)[i], "_sig.txt"), quote=F, row.names=F, col.names=F)
  write.table(names(which(abs(GTX_coef[,i])>quantile(abs(GTX_coef[,i]), 0.95))), file=paste0("GTX_", colnames(GTX_coef)[i], "_top05.txt"), quote=F, row.names=F, col.names=F)
}

contrast_list <- colnames(GTX_coef)[-(1:6)]

for (i in contrast_list){
  current_padj <- GTX_padj[,paste0(i,"_padj")]
  current_coef <- GTX_coef[!is.na(current_padj),i]
  current_padj <- current_padj[!is.na(current_padj)]
  
  write.table(names(current_coef[current_padj < 0.05 & current_coef > 0]), file=paste0("GTX_", i, "_positive_sig.txt"), quote=F, row.names=F, col.names=F)
  write.table(names(current_coef[current_padj < 0.05 & current_coef < 0]), file=paste0("GTX_", i, "_negative_sig.txt"), quote=F, row.names=F, col.names=F)

  write.table(names(which(GTX_coef[,i]<quantile(GTX_coef[,i], 0.01))), file=paste0("GTX_", i, "_negative_top01.txt"), quote=F, row.names=F, col.names=F)
  write.table(names(which(GTX_coef[,i]>quantile(GTX_coef[,i], 0.99))), file=paste0("GTX_", i, "_positive_top01.txt"), quote=F, row.names=F, col.names=F)
}

#CC003 vs CC017 GxT genes with coefficient
write.csv(as.matrix(GTX_coef[sapply(GTX_padj[,"CC017_vs_CC003_padj"]<0.05, isTRUE),"CC017_vs_CC003"]),
          file="CC017_vs_CC003_GTX_with_coefs.csv", quote=F)

#---------------

#GT model without pair information (but including pair as batch covariate, with sex)

dds_nopair <- DESeqDataSetFromMatrix(countData = t(df),
                                     colData = info,
                                     design= ~ sex + pair + strain + rx)

X_gt_nopair <- X_gt[,c(1:6,24),drop=F]
X_t_nopair <- cbind(1,X_gt_nopair[,7,drop=F])
X_g_nopair <- X_gt_nopair[,-7,drop=F]

Z <- model.matrix(~1 + sex + pair, info)[,-1]

X_gt_nopair <- cbind(X_gt_nopair, Z)
X_t_nopair <- cbind(X_t_nopair, Z)
X_g_nopair <- cbind(X_g_nopair, Z)

dds_nopair <- DESeq(dds_nopair, full=X_gt_nopair, reduced=X_t_nopair, test="LRT")

save(dds_nopair, file="dds_LRT_gt_t_nopair.RData")

#number of significant G genes using different models
sum(apply(G_padj, 1, function(x){any(x<0.05, na.rm=T)}))
sum(apply(G_padj_in_GTX, 1, function(x){any(x<0.05, na.rm=T)}))
sum(results(dds_nopair)$padj<0.05, na.rm=T)

#number of significant G genes in all models
sum(apply(G_padj_in_GTX, 1, function(x){any(x<0.05, na.rm=T)}) &
      apply(G_padj, 1, function(x){any(x<0.05, na.rm=T)}) & 
      results(dds_nopair)$padj<0.05)

#---------------

colnames(G_coef) <- sapply(colnames(G_coef), function(x){paste0(x, "_G")})
colnames(G_padj) <- sapply(colnames(G_padj), function(x){paste0(x, "_G")})

colnames(GTX_coef) <- sapply(colnames(GTX_coef), function(x){paste0(x, "_GTX")})
colnames(GTX_padj) <- sapply(colnames(GTX_padj), function(x){paste0(x, "_GTX")})

output <- cbind(output, G_coef, G_padj, GTX_coef, GTX_padj, G_coef_in_GTX, G_padj_in_GTX)

output <- cbind(output[,1:6], results(dds_nopair)$padj, output[,-c(1:6)])
colnames(output)[7] <- "padj_GT_vs_T_unpaired"

#---------------

####################
# Correlation with all data (FA and O3)
pheno_corr <- round(cor(pheno, use="pairwise.complete.obs", method="spearman"), 2)

gene_pheno_corr <- matrix(NA, ncol(df_vst), ncol(pheno))
rownames(gene_pheno_corr) <- colnames(df_vst)
colnames(gene_pheno_corr) <- colnames(pheno)

counter <- 0
for (i in rownames(gene_pheno_corr)){
  counter <- counter + 1
  print(counter)
  for (j in colnames(gene_pheno_corr)){
    gene_pheno_corr[i,j] <- cor(df_vst[,i], pheno[,j], use="pairwise.complete.obs", method="spearman")
  }
}

save(gene_pheno_corr, file="ozone_gene_pheno_corr.RData")

output <- cbind(output, gene_pheno_corr)

##########
# Correlation with baseline phenotypes and gene expression only
df_vst_FA <- df_vst[info$rx=="FA",]
rownames(df_vst_FA) <- info$pair_ID[info$rx=="FA"]

pair_link <- sapply(rownames(df_vst_FA), function(x){which(info$pair_ID==x & info$rx=="O3")})
pheno_O3 <- pheno[pair_link,]
rownames(pheno_O3) <- names(pair_link)

gene_pheno_corr <- matrix(NA, ncol(df_vst_FA), ncol(pheno_03))
rownames(gene_pheno_corr) <- colnames(df_vst_FA)
colnames(gene_pheno_corr) <- colnames(pheno)

counter <- 0
for (i in rownames(gene_pheno_corr)){
  counter <- counter + 1
  print(counter)
  for (j in colnames(gene_pheno_corr)){
    gene_pheno_corr[i,j] <- cor(df_vst_FA[,i], pheno_O3[,j], use="pairwise.complete.obs", method="spearman")
  }
}

colnames(gene_pheno_corr) <- paste0(colnames(pheno), "_FA")

save(gene_pheno_corr, file="ozone_gene_pheno_corr_FA.RData")

output <- cbind(output, gene_pheno_corr)

##########
#correlations between treated genes and phenotypes
df_vst_O3 <- df_vst[info$rx=="O3",]
rownames(df_vst_O3) <- info$pair_ID[info$rx=="O3"]

pair_link <- sapply(rownames(df_vst_O3), function(x){which(info$pair_ID==x & info$rx=="O3")})
pheno_O3 <- pheno[pair_link,]
rownames(pheno_O3) <- names(pair_link)

gene_pheno_corr <- matrix(NA, ncol(df_vst_O3), ncol(pheno_O3))
rownames(gene_pheno_corr) <- colnames(df_vst_O3)
colnames(gene_pheno_corr) <- colnames(pheno)

counter <- 0
for (i in rownames(gene_pheno_corr)){
  counter <- counter + 1
  print(counter)
  for (j in colnames(gene_pheno_corr)){
    gene_pheno_corr[i,j] <- cor(df_vst_O3[,i], pheno_O3[,j], use="pairwise.complete.obs", method="spearman")
  }
}

colnames(gene_pheno_corr) <- paste0(colnames(pheno), "_O3")

save(gene_pheno_corr, file="ozone_gene_pheno_corr_O3.RData")

output <- cbind(output, gene_pheno_corr)

#---------------

write.csv(output, file="ozone_DE_corr_merged.csv")

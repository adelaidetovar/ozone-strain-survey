source("0_functions.R")
source("1_data_cleaning.R")

output_dir <- file.path("../output/wgcna/")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

setwd(output_dir)
rm(output_dir)

#if using R < 4.0
options(stringsAsFactors = FALSE)

#allow multi-threading within WGCNA
enableWGCNAThreads()

load("../gene_expr/ozone_vst.RData")

###############
#### WGCNA ####
###############

datExpr0 <- df_vst

#check for excessive missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#sample clustering to detect outliers
sampleTree = hclust(dist(datExpr0), method = "average")
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=F)
plot(sampleTree)

#PCA to detect outliers
pca <- prcomp(datExpr0)
factoextra::fviz_pca_biplot(pca, geom="point")

#cleaned data (no change)
datExpr <- datExpr0

####################
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

genes <- colnames(df_vst)
G_list <- getBM(filters= "external_gene_name", attributes= c("ensembl_gene_id", "external_gene_name", "description", "entrezgene_id"), values=genes, mart=mart)

#this step does not correctly handle one to many matches, just picks the first one
G_list <- G_list[sapply(genes, match, table=G_list$external_gene_name),]

#---------------------

#call the network topology analysis function
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

sft = pickSoftThreshold(datExpr, verbose = 5,
                        networkType = "signed",
                        corFnc = bicor, corOptions = list(maxPOutliers = 0.05),
                        powerVector = powers)

#plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

for (selected.power in c(9,12)){
  print(selected.power)
  #one-step network construction and module detection
  net = blockwiseModules(datExpr, power = selected.power,
                         corType = "bicor",
                         maxPOutliers = 0.05,
                         TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = paste0("TOM_ozone_", selected.power),
                         verbose = 3,
                         maxBlockSize = 20000)

  save(net, file=paste0("net_", selected.power, ".RData"))

  #open a graphics window
  sizeGrWindow(12, 9)
  #convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  #plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)

  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]]
  save(MEs, moduleLabels, moduleColors, geneTree, file = paste0("WGCNA_output_ozone_",selected.power,".RData"))

  #compute module eigengenes
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)

  save(MEs, file=paste0("module_eigengenes_", selected.power, ".RData"))
}


selected.power <- 12

load(paste0("net_", selected.power, ".RData"))
load(paste0("WGCNA_output_ozone_",selected.power,".RData"))
load(paste0("module_eigengenes_", selected.power, ".RData"))

#open a graphics window
sizeGrWindow(12, 9)
#convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
#plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

print(-sort(-table(moduleColors)))
print(length(unique(moduleColors)))

#correlations between eigengenes and phenotypes
datTraits <- pheno
datTraits <- datTraits[,c(1,4,6:9,12,14,15,17,18)]
datTraits$mouse_no <- rownames(datTraits)
datTraits$mouse_no <- sapply(datTraits$mouse_no, 
                             function(x){unlist(strsplit(x, split="X"))[2]})
datTraits <- merge(bal[,c(1,14)], datTraits, by="mouse_no")
datTraits$mouse_no <- paste0("X",datTraits$mouse_no)
rownames(datTraits) <- datTraits$mouse_no
datTraits <- datTraits[,c(3:9,13,10,11,2)]
datTraits <- datTraits[rownames(MEs),]

#Quantifying module-trait associations
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

traitNames <- c("Eotaxin", "IL-10", "IL-6", "IP-10", "KC", "LIX", "MIP-1b",
                "Protein Conc.", "Percent Neuts", "Num. Neuts", "Total Cells")

pdf(file=paste0("module_correlations_", selected.power, ".pdf"), 
    height = 9, width = 9)

# Display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = traitNames,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.75,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

####################
#Gene relationship to traits and important modules: Gene Significance and Module Membership

# names (colors) of the modules
modNames = substring(names(MEs), 3)

#calculate correlations with eigengenes
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#calculate correlations with traits
trait.correlations <- cor(datExpr, datTraits, use="pairwise.complete.obs", method="spearman")

#output module membership and trait correlations
geneInfo <- data.frame(G_list, moduleColor=moduleColors, moduleMembership=NA, trait.correlations)

for (i in 1:nrow(geneInfo)){
  geneInfo[i, "moduleMembership"] <- geneModuleMembership[i,paste0("MM", geneInfo[i, "moduleColor"])]
}

geneInfo$external_gene_name <- genes

save(geneInfo, file=paste0("WGCNA_table_", selected.power, ".RData"))

geneInfo <- geneInfo[order(geneInfo$moduleColor, -abs(geneInfo$moduleMembership)),]

write.csv(geneInfo, file=paste0("module_membership_ozone_", selected.power, ".csv"), row.names=F)

#---------------------

# Enrichment analysis for selected modules

dbs <- c("KEGG_2019_Mouse", "MGI_Mammalian_Phenotype_Level_4_2019", "MSigDB_Hallmark_2020",
         "ProteomicsDB_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2018", 
         "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "Genome_Browser_PWMs",
         "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs")

# pink module
write.xlsx(enrichr(genes = geneInfo$external_gene_name[geneInfo$moduleColor=="pink"], dbs),"enrichR_pink.xlsx")
# red module
write.xlsx(enrichr(genes = geneInfo$external_gene_name[geneInfo$moduleColor=="red"], dbs),"enrichR_red.xlsx")
# blue module
write.xlsx(enrichr(genes = geneInfo$external_gene_name[geneInfo$moduleColor=="blue"], dbs),"enrichR_blue.xlsx")
# midnightblue module
write.xlsx(enrichr(genes = geneInfo$external_gene_name[geneInfo$moduleColor=="midnightblue"], dbs),"enrichR_midnightblue.xlsx")
# brown module
write.xlsx(enrichr(genes = geneInfo$external_gene_name[geneInfo$moduleColor=="brown"], dbs),"enrichR_brown.xlsx")

#---------------------

# Hub gene analysis

## intramodular connectivity

MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)

## module membership

datKME=signedKME(datExpr, MEs, outputColumnName = "MM.")
colorlevels = unique(moduleColors)

## gene significance/correlation for protein
GS1_protein=as.numeric(WGCNA::cor(datTraits$protein_conc,datExpr,use="p"),2)
GS1_protein=abs(GS1_protein)
ModuleSignificance_protein=tapply(GS1_protein, moduleColors, mean, na.rm=T)

## gene significance/correlation for percent neutrophils
GS1_perneuts=as.numeric(WGCNA::cor(datTraits$per_neu,datExpr,use="p"),2)
GS1_perneuts=abs(GS1_perneuts)
ModuleSignificance_perneuts=tapply(GS1_perneuts, moduleColors, mean, na.rm=T)

## gene significance/correlation for IL-10
GS1_il10=as.numeric(WGCNA::cor(datTraits$il10,datExpr,use="p"),2)
GS1_il10=abs(GS1_il10)
ModuleSignificance_il10=tapply(GS1_il10, moduleColors, mean, na.rm=T)

## gene significance/correlation for KC
GS1_kc=as.numeric(WGCNA::cor(datTraits$kc,datExpr,use="p"),2)
GS1_kc=abs(GS1_kc)
ModuleSignificance_kc=tapply(GS1_kc, moduleColors, mean, na.rm=T)

##########
#### neuts ----------
##########
# gene-eigengene (kME) and gene-trait (GS1n) relationship
GS1p <- as.data.frame(GS1_protein)
rownames(GS1p) <- GS1p$gene <- colnames(datExpr)
GS1p$moduleColors <- moduleColors
GTp = merge(datKME,GS1p,by=c("gene","moduleColors"))

GTp <- GTp %>%
  mutate(red.label = ifelse(gene %in% red_genes$gene, gene, ""),
         brown.label = ifelse(gene %in% brown_genes$gene, gene, ""),
         blue.label = ifelse(gene %in% blue_genes$gene, gene, ""),
         pink.label = ifelse(gene %in% pink_genes$gene, gene, ""),
         midnightblue.label = ifelse(gene %in% midnightblue_genes$gene,gene, ""),
         red.color = ifelse(gene %in% red_genes$gene, TRUE, FALSE),
         brown.color = ifelse(gene %in% brown_genes$gene, TRUE, FALSE),
         blue.color = ifelse(gene %in% blue_genes$gene, TRUE, FALSE),
         pink.color = ifelse(gene %in% pink_genes$gene, TRUE, FALSE),
         midnightblue.color = ifelse(gene %in% midnightblue_genes$gene, TRUE, FALSE))

brown_protein <- wgcnaCorPlot(df = GTp, modCol = "brown",
                              modMem = MM.brown, pheno = GS1_protein,
                              label = brown.label, color = brown.color,
                              yaxis = "Correlation with BAL Protein",
                              title = "Brown Module")


ggsave(brown_protein, filename = "brown_corprot.png", units = "in", width = 6, height = 6)

pink_protein <- wgcnaCorPlot(df = GTp, modCol = "pink",
                             modMem = MM.pink, pheno = GS1_protein,
                             label = pink.label, color = pink.color,
                             yaxis = "Correlation with BAL Protein",
                             title = "Pink Module")

ggsave(pink_protein, filename = "pink_corprot.png", units = "in", height = 6, width = 6)

##########
#### neuts ----------
##########
# gene-eigengene (kME) and gene-trait (GS1n) relationship
GS1n <- as.data.frame(GS1_perneuts)
rownames(GS1n) <- GS1n$gene <- colnames(datExpr)
GS1n$moduleColors <- moduleColors
GTn = merge(datKME,GS1n,by=c("gene","moduleColors"))

GTn <- GTn %>%
  mutate(red.label = ifelse(gene %in% red_genes$gene, gene, ""),
         brown.label = ifelse(gene %in% brown_genes$gene, gene, ""),
         blue.label = ifelse(gene %in% blue_genes$gene, gene, ""),
         pink.label = ifelse(gene %in% pink_genes$gene, gene, ""),
         midnightblue.label = ifelse(gene %in% midnightblue_genes$gene,gene, ""),
         red.color = ifelse(gene %in% red_genes$gene, TRUE, FALSE),
         brown.color = ifelse(gene %in% brown_genes$gene, TRUE, FALSE),
         blue.color = ifelse(gene %in% blue_genes$gene, TRUE, FALSE),
         pink.color = ifelse(gene %in% pink_genes$gene, TRUE, FALSE),
         midnightblue.color = ifelse(gene %in% midnightblue_genes$gene, TRUE, FALSE))

pink_neut <- wgcnaCorPlot(df = GTn, modCol = "pink",
                          modMem = MM.pink, pheno = GS1_perneuts,
                          label = pink.label, color = pink.color,
                          yaxis = "Correlation with BAL Percent Neutrophils",
                          title = "Pink Module")

ggsave(pink_neut, filename = "pink_corneut.png", units = "in", height = 6, width = 6)

red_neut <- wgcnaCorPlot(df = GTn, modCol = "red",
                         modMem = MM.red, pheno = GS1_perneuts,
                         label = red.label, color = red.color,
                         yaxis = "Correlation with BAL Percent Neutrophils",
                         title = "Red Module")
  
ggsave(red_neut, filename = "red_corneut.png", units = "in", height = 6, width = 6)

midnightblue_neut <- wgcnaCorPlot(df = GTn, modCol = "midnightblue",
                                  modMem = MM.midnightblue, pheno = GS1_perneuts,
                                  label = midnightblue.label, color = midnightblue.color,
                                  yaxis = "Correlation with BAL Percent Neutrophils",
                                  title = "Midnightblue Module")
  
ggsave(midnightblue_neut, filename = "midnightblue_corneut.png", units = "in", width = 6, height = 6)

##########
#### KC ----------
##########
# gene-eigengene (kME) and gene-trait (GS1n) relationship
GS1k <- as.data.frame(GS1_kc)
rownames(GS1k) <- GS1k$gene <- colnames(datExpr)
GS1k$moduleColors <- moduleColors
datKME$gene <- rownames(datKME)
datKME$moduleColors <- moduleColors
GTk = merge(datKME,GS1k,by=c("gene","moduleColors"))

# pick filter values
GTk %>% filter(moduleColors == "pink") %>%
  ggplot(aes(x = MM.pink, y = GS1_kc)) + geom_point()

GTk <- GTk %>%
  mutate(red.color = ifelse(gene %in% red_genes$gene, TRUE, FALSE),
         brown.color = ifelse(gene %in% brown_genes$gene, TRUE, FALSE),
         blue.color = ifelse(gene %in% blue_genes$gene, TRUE, FALSE),
         pink.color = ifelse(gene %in% pink_genes$gene, TRUE, FALSE),
         midnightblue.color = ifelse(gene %in% midnightblue_genes$gene, TRUE, FALSE),
         red.label = ifelse(gene %in% red_genes$gene, gene, ""),
         brown.label = ifelse(gene %in% brown_genes$gene, gene, ""),
         blue.label = ifelse(gene %in% blue_genes$gene, gene, ""),
         pink.label = ifelse(gene %in% pink_genes$gene, gene, ""),
         midnightblue.label = ifelse(gene %in% midnightblue_genes$gene, gene, ""))

pink_kc <- wgcnaCorPlot(df = GTk, modCol = "pink",
                        modMem = MM.pink, pheno = GS1_kc,
                        label = pink.label, color = pink.color,
                        yaxis = "Correlation with BAL KC",
                        title = "Pink Module")

ggsave(pink_kc, filename = "pink_corkc.png", units = "in", height = 6, width = 6)

##########
#### IL-10 ----------
##########
# gene-eigengene (kME) and gene-trait (GS1n) relationship
GS1t <- as.data.frame(GS1_il10)
rownames(GS1t) <- GS1t$gene <- colnames(datExpr)
GS1t$moduleColors <- moduleColors
datKME$gene <- rownames(datKME)
datKME$moduleColors <- moduleColors
GTt = merge(datKME,GS1t,by=c("gene","moduleColors"))

GTt <- GTt %>%
  mutate(red.color = ifelse(gene %in% red_genes$gene, TRUE, FALSE),
         brown.color = ifelse(gene %in% brown_genes$gene, TRUE, FALSE),
         blue.color = ifelse(gene %in% blue_genes$gene, TRUE, FALSE),
         pink.color = ifelse(gene %in% pink_genes$gene, TRUE, FALSE),
         midnightblue.color = ifelse(gene %in% midnightblue_genes$gene, TRUE, FALSE),
         red.label = ifelse(gene %in% red_genes$gene, gene, ""),
         brown.label = ifelse(gene %in% brown_genes$gene, gene, ""),
         blue.label = ifelse(gene %in% blue_genes$gene, gene, ""),
         pink.label = ifelse(gene %in% pink_genes$gene, gene, ""),
         midnightblue.label = ifelse(gene %in% midnightblue_genes$gene, gene, ""))

blue_il10 <- wgcnaCorPlot(df = GTt, modCol = "blue",
                          modMem = MM.blue, pheno = GS1_il10,
                          label = blue.label, color = blue.color,
                          yaxis = "Correlation with BAL IL-10",
                          title = "Blue Module")

ggsave(blue_il10, filename = "blue_coril10.png", units = "in", height = 6, width = 6)

# module connectivity for modules discussed in paper
## red, brown, blue, pink, midnightblue 

# filter genes for each based on intramodular connectivity and module membership, use for labels
Alldegrees1$gene <- rownames(Alldegrees1)
datKME$gene <- rownames(datKME)
Alldegrees1$moduleColors <- moduleColors
datKME$moduleColors <- moduleColors

FilterGenes = merge(Alldegrees1, datKME, by =c("gene","moduleColors"))

red_genes <- FilterGenes %>% filter(moduleColors == "red" & MM.red > .75) %>% arrange(desc(kWithin, MM.red)) %>% top_n(10, wt = kWithin) %>% dplyr::select(gene)
blue_genes <- FilterGenes %>% filter(moduleColors == "blue" & MM.blue > .75) %>% arrange(desc(kWithin, MM.blue)) %>% top_n(20, wt = kWithin) %>% select(gene)
brown_genes <- FilterGenes %>% filter(moduleColors == "brown" & MM.brown>.75) %>% arrange(desc(kWithin)) %>% top_n(10, wt = kWithin) %>% select(gene)
pink_genes <- FilterGenes %>% filter(moduleColors == "pink" & MM.pink>.75) %>% arrange(desc(kWithin)) %>% top_n(10, wt = kWithin) %>% select(gene)
midnightblue_genes <- FilterGenes %>% filter(moduleColors == "midnightblue" & MM.midnightblue>.5) %>% arrange(desc(kWithin)) %>% top_n(10, wt = kWithin) %>% select(gene)

FilterGenes <- FilterGenes %>%
  mutate(red.color = ifelse(gene %in% red_genes$gene, TRUE, FALSE),
         brown.color = ifelse(gene %in% brown_genes$gene, TRUE, FALSE),
         blue.color = ifelse(gene %in% blue_genes$gene, TRUE, FALSE),
         pink.color = ifelse(gene %in% pink_genes$gene, TRUE, FALSE),
         midnightblue.color = ifelse(gene %in% midnightblue_genes$gene, TRUE, FALSE),
         red.label = ifelse(gene %in% red_genes$gene, gene, ""),
         brown.label = ifelse(gene %in% brown_genes$gene, gene, ""),
         blue.label = ifelse(gene %in% blue_genes$gene, gene, ""),
         pink.label = ifelse(gene %in% pink_genes$gene, gene, ""),
         midnightblue.label = ifelse(gene %in% midnightblue_genes$gene, gene, ""))

hubGenePlot <- function(df, modCol, modMem, label, color, title){
  df %>% filter(moduleColors == {{modCol}} & {{modMem}} > 0) %>%
    ggplot(aes(x = kWithin, y = {{modMem}}, label = {{label}})) +
    geom_point(aes(color = {{color}})) +
    geom_text_repel(aes(fontface="italic"), ylim=c(0.7, 1), size = 4, max.overlaps = 50) +
    theme_linedraw(base_size = 12) + scale_color_manual(values = c("black", "red")) +
    labs(y = "Module Membership", x = "Intramodule Connectivity") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) + ggtitle(title)
}

red_hubgene <- hubGenePlot(df = FilterGenes, modCol = "red",
                           modMem = MM.red, label = red.label,
                           color = red.color, title = "Red Module")

ggsave(red_hubgene, file = "red_hubgene.png", units = "in", height = 6, width = 6)

blue_hubgene <- hubGenePlot(df = FilterGenes, modCol = "blue",
                           modMem = MM.blue, label = blue.label,
                           color = blue.color, title = "Blue Module")

ggsave(blue_hubgene, file = "blue_hubgene.png", units = "in", height = 6, width = 6)

midnightblue_hubgene <- hubGenePlot(df = FilterGenes, modCol = "midnightblue",
                           modMem = MM.midnightblue, label = midnightblue.label,
                           color = midnightblue.color, title = "Midnightblue Module")

ggsave(midnightblue_hubgene, file = "midnightblue_hubgene.png", units = "in", height = 6, width = 6)

pink_hubgene <- hubGenePlot(df = FilterGenes, modCol = "pink",
                           modMem = MM.pink, label = pink.label,
                           color = pink.color, title = "Pink Module")

ggsave(pink_hubgene, file = "pink_hubgene.png", units = "in", height = 6, width = 6)

brown_hubgene <- hubGenePlot(df = FilterGenes, modCol = "brown",
                           modMem = MM.brown, label = brown.label,
                           color = brown.color, title = "Brown Module")

ggsave(brown_hubgene, file = "brown_hubgene.png", units = "in", height = 6, width = 6)
# contributed by Adelaide Tovar

source("0_functions.R")
source("1_data_cleaning.R")

output_dir <- file.path("../output/plots/")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

setwd(output_dir)
rm(output_dir)

##################
#### Plotting ####
##################

# Albumin plot

albumin_plot <- phenoPlot(albumin, albumin_conc, ylab = "Concentration (ng/ml)", title = "BAL Albumin")  +
  scale_y_continuous(breaks=c(seq(0,1300,by=150))) + 
  coord_cartesian(ylim=c(0,1300), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -60, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -120, label = me93_strains, size= 7, fontface=2)

ggsave(plot=albumin_plot,'albumin_plot.png',
       dpi=300, height=7.5, width=10, units="in")

# Total protein plot

protein_plot <- phenoPlot(protein, protein_conc, ylab = "Concentration (ug/ml)", title = NULL) +
  scale_y_continuous(breaks=c(seq(0,3500,by=500))) +
  coord_cartesian(ylim=c(0,3500), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.04)*3500, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.09)*3500, label = me93_strains, size= 7, fontface=2)

ggsave(plot = protein_plot, 'total_protein_plot.png',
       dpi = 300, height = 7.5, width = 10, units = "in")

# Total Cell Number Plot

tot_num_plot <- phenoPlot(bal, total_cells, ylab = parse(text=paste0('"# of cells"','~ (10^6)')), title = NULL) +
  scale_y_continuous(breaks=c(seq(0,1.2e6,by=1e5)),labels=c(seq(0,1.2,by=0.1))) +
  coord_cartesian(ylim=c(-2e4,1.2e6), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.06)*1.22e6, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.105)*1.22e6, label = me93_strains, size= 7, fontface=2)

ggsave(plot = tot_num_plot, 'tot_num_plot.png',
       dpi = 300, height = 7.5, width = 10, units = "in")

# Neutrophil Number Plot

neut_num_plot <- phenoPlot(bal, no_neu, ylab = parse(text=paste0('"# of neutrophils"','~ (10^5)')), title = NULL) +
  scale_y_continuous(breaks=c(seq(0,6e5,by=1e5)),labels=c(0,1,2,3,4,5,6)) +
  coord_cartesian(ylim=c(-1e4,6e5), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.06)*6e5, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.105)*6e5, label = me93_strains, size= 7, fontface=2)

ggsave(plot = neu_num_plot, 'neut_num_plot.png',
       dpi = 300, height = 7.5, width = 10, units = "in")

# Neutrophil Percentage Plot

neut_per_plot <- phenoPlot(bal, per_neu, ylab = "% neutrophils in total BAL cells",
                           title = "% neutrophils in BAL") +
  scale_y_continuous(breaks=c(seq(0,50,by=10))) +
  coord_cartesian(ylim=c(-2,53),clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -4.2, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -6.9, label = me93_strains, size= 7, fontface=2)

ggsave(plot = neu_per_plot, 'neut_per_plot.png',
       dpi = 300, height = 7.5, width = 10, units = "in")

# Macrophage percentage plot

mac_per_plot <- phenoPlot(bal, per_macs, ylab = "% macrophages in total BAL cells",
                          title = "% macrophages in BAL") +
  scale_y_continuous(breaks=c(seq(50,100,by=10)))  +
  guides(shape = "none",
         fill = "none",
         color = guide_legend(override.aes = list(shape = 21, stroke = 1)))+
  coord_cartesian(ylim=c(45,102),clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = 43, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = 40, label = me93_strains, size= 7, fontface=2)

ggsave(plot = mac_per_plot, 'mac_per_plot.png',
       dpi = 300, height = 7.5, width = 10, units = "in")

# Macrophage number plot

mac_num_plot <- phenoPlot(bal, no_macs, ylab = parse(text=paste0('"# of macrophages"','~ (10^5)')),
                          title = "Number of macrophages in BAL") +
  scale_y_continuous(breaks=c(seq(0,9e5,by=1e5)),labels=c(0,1,2,3,4,5,6,7,8,9)) +
  coord_cartesian(ylim=c(-500,9e5),clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -35500, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -80500, label = me93_strains, size= 7, fontface=2)

ggsave(plot = mac_num_plot, 'mac_num_plot.png',
       dpi = 300, height = 7.5, width = 10, units = "in")

# IL-6 plot

ylab =  "Concentration (pg/mL)"

il6_plot <- phenoPlot(cytokines, il6, ylab = ylab, title = NULL) +
  coord_cartesian(ylim=c(0,800), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*800, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*800, label = me93_strains, size= 7, fontface=2)

ggsave(plot = il6_plot, 'il6_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

# KC (CXCL1) plot

kc_plot <- phenoPlot(cytokines, kc, ylab =  ylab, title = NULL) +
  coord_cartesian(ylim=c(0,190), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*190, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*190, label = me93_strains, size= 7, fontface=2)

ggsave(plot = kc_plot, 'cxcl1_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

# IL-10 plot

il10_plot <- phenoPlot(cytokines, il10, ylab = ylab, title = NULL) +
  coord_cartesian(ylim=c(0,250), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*250, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*250, label = me93_strains, size= 7, fontface=2)

ggsave(plot = il10_plot, 'il10_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

# IP-10 plot

ip10_plot <- phenoPlot(cytokines, ip10, ylab = ylab, title = NULL) +
  coord_cartesian(ylim=c(0,60), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*60, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*60, label = me93_strains, size= 7, fontface=2)

ggsave(plot = ip10_plot, 'ip10_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

# eotaxin (CCL11) plot

eotaxin_plot <- phenoPlot(cytokines, eotaxin, ylab = ylab, title = NULL) +
  guides(shape = "none",
         fill = "none",
         color = guide_legend(override.aes = list(shape = c(21), stroke = 1))) +
  coord_cartesian(ylim=c(0,280), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*280, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*280, label = me93_strains, size= 7, fontface=2)

ggsave(plot = eotaxin_plot, 'eotaxin_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

# LIX (CXCL5) plot

lix_plot <- phenoPlot(cytokines, lix, ylab = ylab, title = NULL) +
  coord_cartesian(ylim=c(0,200), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*200, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*200, label = me93_strains, size= 7, fontface=2)

ggsave(plot = lix_plot, 'lix_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

# MIP-1B

mip1b_plot <- phenoPlot(cytokines, mip1b, ylab = ylab, title = NULL) +
  coord_cartesian(ylim=c(0,55), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -(0.05)*55, label = me93_labels, size = 5.5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -(0.1083)*55, label = me93_strains, size= 7, fontface=2)

ggsave(plot = mip1b_plot, 'mip1b_plot.png',
       dpi = 300, height = 6, width = 8, units = "in")

############################################
#### Correlation of cells and cytokines ####
############################################

# Select cytokines of interest: eotaxin, il10, il6, ip10, kc, lix, mip1b
cytokines_corr <- cytokines[,c("mouse_no", "strain", "tx", "sex", "batch",
                               "eotaxin", "il10", "il6", "ip10", "kc", "lix", "mip1b")]

# Merge with BAL data
cytokines_corr <- left_join(cytokines_corr, bal[,c(1, 8, 12, 14)], by = "mouse_no")

# Merge with protein data
cytokines_corr <- left_join(cytokines_corr, protein[,c(1, 6)], by = "mouse_no")

# Change column names for plotting
colnames(cytokines_corr) <- c("mouse_no","strain","tx","sex","batch","Eotaxin",
                              "IL-10","IL-6","IP-10","CXCL1","LIX","MIP-1b",
                              "Percent Neuts","Num. Neuts","Num. Cells","Protein Conc.")

# Compute correlation matrices, both complete observations and with all complete pairwise observations
phenos_corr <- cor_mat(cytokines_corr[,c(6:16)],method="spearman")
phenos_pcorr <- cor_pmat(cytokines_corr[,c(6:16)],method="spearman",exact=FALSE)

# Correlation plot
phenos_corr_plot <- ggcorrplot(phenos_corr, outline.col = "white",
                               ggtheme = ggplot2::theme_linedraw(base_size=18),
                               type = "lower",
                               colors = c("#00b9ff", "white", "#f23f00"),
                               show.diag=TRUE, lab=TRUE, lab_size=4,tl.srt = 40, insig="blank") +
  theme(axis.text=element_text(size=12, color="black"),
        plot.title=element_text(hjust=0.5,face="bold")) +
  ggtitle("Correlation between Inflammation/Injury Phenotypes\nand Cytokine Measurements")

ggsave(phenos_corr_plot, file = 'phenos_corr_plot.png', dpi = 300, units="in", height = 8, width = 9)

#########################################################
#### Principal Components Analysis for Top 500 Genes ####
#########################################################

cts <- df
coldata <- info

# Data normalization
dds <- DESeqNormalize(counts = cts, coldata = coldata, designFormula = ~ strain*rx)

# VSD, select top 500 most variable genes
vsd <- assay(vsd)
vsd.var <- varFilter(vsd, variableGenes = 500)
vsd.thres.var <- varFilter(vsd.var, variableGenes = 500)

pca = prcomp(t(vsd.thres.var))
print(summary(pca))

pcaData = as.data.frame(pca$x)
pcaData$mouse_no=rownames(pcaData)
pcaData=merge(pcaData, coldata)
percentVar = round(100 * (pca$sdev^2 / sum(pca$sdev^2)),digits = 2)

pcaPlot <- pcaData %>% ggplot(aes(x = PC1, y = PC2, color = factor(strain), shape = factor(rx))) +
  geom_point(size = 4,
             aes(fill = factor(strain), alpha = sex)) +
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values=c(21,24),
                     name="Treatment",
                     labels=c("FA",expression("O"[3]))) +
  scale_color_manual(values = me93_colors,
                     name = "Strain") +
  scale_fill_manual(values = me93_colors) +
  scale_alpha_manual(values=c("F"=1, "M"=0), name = "Sex") +
  theme_linedraw(base_size = 20) +
  theme(plot.margin=unit(c(1,1,1,1),"lines")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  guides(color = guide_legend(override.aes = list(color = me93_colors)),
         shape = guide_legend(override.aes = list(shape = c(16,17))),
         fill = "none",
         alpha = guide_legend(override.aes = list(shape = c(15, 0), alpha = 1, stroke = 1)))

ggsave(plot = pcaPlot, 'rnaseq_pca.png', dpi = 300, height = 8, width = 10, units = "in")

#################
#### MA Plot ####
#################

load(file = "../gene_expr/resLFC_LRT_gt_g.RData")
table <- data.frame(resLFC)

# all diff expr genes
diffexp <- table[which(table$padj <= 0.05),]
diffexp$gene_name <- rownames(diffexp)

# all diff expr genes with LFC cutoff
down <- rownames(diffexp[diffexp$log2FoldChange <= -1, ])
up <- rownames(diffexp[diffexp$log2FoldChange >= 1, ])
lg.genes <- union(down, up)

# only padj and not LFC diff expr genes
sm.genes <- setdiff(diffexp$gene_name, lg.genes)

# highest DE genes
top.diffexp <- diffexp[order(abs(diffexp$log2FoldChange), decreasing = TRUE)[1:10], ]
top.diffexp <- rownames(top.diffexp)

x <- log10(table$baseMean)
y <- table$log2FoldChange
gene_name <- rownames(table)
d <- densCols(x, y, nbin = 100,
              colramp = colorRampPalette((brewer.pal(9,"Greys")[-c(1:5)])))
df <- data.frame(x, y, d, gene_name)

p <- ggplot(df, aes(x = x, y = y, label = gene_name)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1)) +
  scale_color_identity() + 
  labs(x = "Mean of normalized counts", y = expression(paste("log"[2]," fold change"))) +
  scale_x_continuous(labels = parse(text=c("10", "100", "10^4", "10^6"))) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8))

maplot.basic <- p +
  geom_point(aes(x, y, col = d), size = 1.3, shape = 16) +
  geom_abline(aes(intercept = -1, slope = 0), size = 0.8, linetype = 3) +
  geom_hline(yintercept = 0, size = 0.8) +
  geom_abline(aes(intercept = 1, slope = 0), size = 0.8, linetype = 3)

maplot.highlight <- p +
  geom_point(data = df[!(df$gene_name %in% diffexp$gene_name),],
             aes(x, y, col = "gray60"), size = 1.3, shape = 16, alpha = 0.25) +
  geom_point(data = df[df$gene_name %in% sm.genes,],
             aes(x, y, col = "#9999ff"), size = 1.3, shape = 16, alpha = 0.65) +
  geom_point(data = df[df$gene_name %in% lg.genes,],
             aes(x, y, col = "blue"), size = 1.3, shape = 16, alpha = 0.65) +
  geom_hline(yintercept = 0, size = 0.8) +
  geom_label_repel(data = df[df$gene_name %in% top.diffexp,],
                   point.padding = 0.2,
                   nudge_x = 0.2,
                   nudge_y = 0.5,
                   min.segment.length = 0,
                   max.overlaps = 10,
                   fill = "white",
                   xlim = c(-5,6),
                   ylim = c(0,7.5),
                   fontface = "italic")

ggsave(plot = maplot.highlight, "maplot_labels.png", dpi = 300, height = 5, width = 7, units = "in")

########################################
#### PCA biplot for genes x strains ####
########################################

load('../gene_expr/GTX_coef.RData')
gxt_coef <- GTX_coef
gxt_coef <- as.data.frame(gxt_coef[,c(1:6)])
load('../gene_expr/GTX_padj.RData')
gxt_padj <- GTX_padj
gxt_padj <- as.data.frame(gxt_padj[,c(1:6)])
gxt_padj_sig <- filter_all(gxt_padj, any_vars(. <= 0.05))

# Generate PCA
pca = prcomp(gxt_coef)
# Rotate coordinates
pca$x <- pca$x * -1
pca$rotation <- pca$rotation * -1

# Get genes with highest loadings
pc_genes <- pca$x
pc_genes <- pc_genes %>%
  as_tibble(rownames = "gene")

top_genes <- pc_genes %>%
  dplyr::select(gene, PC1, PC2) %>%
  pivot_longer(matches("PC"), names_to="PC", values_to="loading") %>%
  group_by(PC) %>%
  arrange(desc(abs(loading))) %>%
  slice(1:17) %>%
  pull(gene) %>%
  unique()

top_loadings <- pc_loadings %>%
  filter(gene %in% top_genes)

# Calculate percent variance explained
percentVar = round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
pc_eigenvalues <- pca$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)),
                         variance = pc_eigenvalues) %>%
  mutate(pct = round(100*(variance/sum(variance))),
         pct_cum = cumsum(pct))

text_genes <- c("Slc40a1", "Plxdc2", "Pla2g7", "Vcan", "Cd74",
                "Card11", "Rxra", "Lrch4", "Hif1a", "Mmp12", "Clec4a2",
                "Ccl9", "Arg1", "Cd86", "Osbpl3", "Marco")

top_genes <- c(top_genes, text_genes)
top_genes <- top_genes[!duplicated(top_genes)]
top_genes <- top_genes[c(1:23, 25:41)]

pc_scores <- pca$x

pc_scores <- pc_scores %>%
  as_tibble(rownames = "gene")

pc_scores <- pc_scores %>%
  mutate(label = case_when(gene %in% top_genes ~ gene,
                           TRUE ~ ""),
         in_text = case_when(gene %in% text_genes ~ TRUE,
                             TRUE ~ FALSE))
pc_loadings <- pca$rotation
pc_loadings <- pc_loadings %>%
  as_tibble(rownames = "strain")

dgxt <- ggplot_build(gxt_pca)

gxt_loadings <- dgxt$data[[2]]
gxt_loadings$label <- dgxt$data[[3]]

gxt_loadings_labels <- gxt_loadings

# PCA Plot
gxt_pca <- ggplot(data = pc_scores, aes(x = PC1, y= PC2)) +
  geom_point(size=1, alpha = 0.5) +
  geom_text_repel(aes(label = label, color = in_text, fontface="italic", segment.color = "black"), min.segment.length = 0, max.overlaps=Inf) +
  scale_color_manual(values = c("#9999ff","blue")) +
  geom_segment(data = gxt_loadings, aes(x = x, y = x, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.1, "in")),
               color = "#FF4A3F") +
  geom_text_repel(data = gxt_loadings, aes(x = xend, y = yend, label = label),
                  size = 6, color = "#FF4A3F") +
  theme_linedraw(base_size=18) +
  theme(legend.position="none") +
  xlab(paste0("PC1: ", round(pc_eigenvalues$pct[1],2), "% variance")) +
  ylab(paste0("PC2: ", round(pc_eigenvalues$pct[2],2), "% variance")) +
  guides(fill = guide_legend(override.aes=list(shape=22,color=NA))) +
  scale_y_continuous(limits = c(-3, 6), breaks = seq(-3, 6, by = 2)) +
  scale_x_continuous(limits = c(-5, 3), breaks = seq(-5, 3, by = 2))

ggsave(plot = gxt_pca, 'gxt_pca.png', dpi = 300, height = 6, width = 6, units = "in")

##########################
#### GxT RNA-seq Plot ####
##########################

#plot genes of interest
plot.df <- info

# Ccl6
plot.df$y <- df[,"Ccl6"]
ccl6_plot <- gxtPlot(plot.df, gene = "Ccl6") +
  coord_cartesian(ylim=c(0,43250), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -0.05*43250, label = me93_labels, size = 5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -0.1083*43250, label = me93_strains, size= 6.5, fontface=2)

# Cxcl3
plot.df$y <- df[,"Cxcl3"]
cxcl3_plot <- gxtPlot(plot.df, gene = "Cxcl3") +
  coord_cartesian(ylim=c(0,875), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -0.05*875, label = me93_labels, size = 5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -0.1083*875, label = me93_strains, size= 6.5, fontface=2)

# Slpi
plot.df$y <- df[,"Slpi"]
slpi_plot <- gxtPlot(plot.df, gene = "Slpi") +
  coord_cartesian(ylim=c(0,2800), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -0.05*2800, label = me93_labels, size = 5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -0.1083*2800, label = me93_strains, size= 6.5, fontface=2)

# Cd36
plot.df$y <- df[,"Cd36"]
cd36_plot <- gxtPlot(plot.df, gene = "Cd36") +
  coord_cartesian(ylim=c(0,22750), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -0.05*22750, label = me93_labels, size = 5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -0.1083*22750, label = me93_strains, size= 6.5, fontface=2)

# Marco
plot.df$y <- df[,"Marco"]
marco_plot <- gxtPlot(plot.df, gene = "Marco") +
  coord_cartesian(ylim=c(0,12000), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -600, label = me93_labels, size = 5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -1300, label = me93_strains, size= 6.5, fontface=2)

# Msr1
plot.df$y <- df[,"Msr1"]
msr1_plot <- gxtPlot(plot.df, gene = "Msr1") +
  coord_cartesian(ylim=c(0,2750), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,12,by=1), y = -0.05*2750, label = me93_labels, size = 5) +
  annotate(geom = "text", x = 1.5 + 2 * (0:5),y = -0.1083*2750, label = me93_strains, size= 6.5, fontface=2)

gxt_plot <- (ccl6_plot + cxcl3_plot + slpi_plot) / (cd36_plot + marco_plot + msr1_plot)

ggsave(gxt_plot, "gxt_plot.png", dpi = 300, units = "in", width = 24, height = 12)
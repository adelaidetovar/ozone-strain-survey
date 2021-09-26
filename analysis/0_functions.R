#####################
#### Environment ####
#####################

#### USING R version 4.1.1

using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

BiocManager::install("GEOquery")

tc_packages <- c("openxlsx", "ggfortify", "ggcorrplot", "rstatix", "dplyr",
                 "tidyverse", "ggrepel", "reshape2", "patchwork", "WGCNA", "enrichR",
                 "ashr", "factoextra", "biomaRt","GEOquery", "readxl", "DESeq2")
using(tc_packages)

tc_colors <- c("#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00")

tc_labels <- c("FA.C57BL/6J"="FA",
                 "O3.C57BL/6J"=expression("O"[3]),
                 "FA.CC059"="FA",
                 "O3.CC059"=expression("O"[3]),
                 "FA.CC017"="FA",
                 "O3.CC017"=expression("O"[3]),
                 "FA.CC025"="FA",
                 "O3.CC025"=expression("O"[3]),
                 "FA.CC039"="FA",
                 "O3.CC039"=expression("O"[3]),
                 "FA.CC003"="FA",
                 "O3.CC003"=expression("O"[3]))

tc_strains <- c("C57BL/6J","CC003","CC017","CC025","CC039","CC059")

output_dir <- file.path("../output/")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

rm(output_dir)

############################
#### Plotting Functions ####
############################

phenoPlot <- function(df, pheno, ylab, title) { ggplot(df,
                       aes(x=interaction(tx,strain),
                           y={{pheno}},
                           fill=strain,
                           group=interaction(tx,strain))) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=interaction(tx,sex), stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(16,17,21,24)) +
  scale_fill_manual(values=tc_labels) +
  theme_linedraw(base_size = 20) +
    labs(x = NULL, y = ylab) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="none",
        plot.margin=unit(c(1,1,3.5,1),"lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle(label=title) +
  scale_x_discrete(labels=tc_labels)
}

gxtPlot <- function(df, gene) {ggplot(df,
                             aes(x=interaction(rx,strain),
                                 y=y,
                                 fill=strain)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=interaction(rx,sex),stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(16,17,21,24)) +
  scale_fill_manual(values=tc_colors) +
  theme_linedraw(base_size=20) +
  xlab(NULL) +
  labs(y="Normalized counts") +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="none",
        plot.margin=unit(c(1,1,3.5,1),"lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_discrete(labels=tc_labels) +
  ggtitle(label=bquote(italic(gene)))
}

wgcnaCorPlot <- function(df, modCol, modMem, pheno, label, color, yaxis, title) {
  df %>% filter(moduleColors == modCol & {{modMem}} > 0 & {{pheno}} > 0) %>%
  ggplot(aes(x = {{modMem}}, y = {{pheno}}, label = {{label}})) +
    geom_point(aes(color = {{color}})) +
    geom_text_repel(aes(fontface="italic")) +
    theme_linedraw(base_size = 12) +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Module Membership", y = yaxis) +
    theme(legend.position = "none",
        plot.title=element_text(hjust=0.5,face="bold")) + 
    ggtitle(title)
  }

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

#######################
#### GEX Functions ####
#######################

sumtozero_contrast <- function(K){
  u <- 1/((K - 1)^(0.5))
  v <- (-1 + (K^(0.5))) * ((K - 1)^(-1.5))
  w <- (K - 2) * v + u
  C <- matrix(-v, K - 1, K - 1)
  diag(C) <- rep(w, K - 1)
  rbind(C, rep(-u, K - 1))
}

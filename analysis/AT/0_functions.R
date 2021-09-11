#####################
#### Environment ####
#####################

#### USING R version 4.1.1

tc_packages <- c("openxlsx", "ggplot2", "cowplot", "ggfortify", "ggcorrplot",
                   "rstatix", "dplyr", "tidyverse", "ggrepel", "reshape2")
lapply(tc_packages, require, character.only=TRUE)

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

############################
#### Plotting Functions ####
############################

phenoPlot <- function(df, pheno, ylab, title) { ggplot(df,
                       aes(x=interaction(tx,strain),
                           y={{pheno}},
                           fill=strain,
                           group=interaction(tx,strain))) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=tx,color = sex, stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=tc_labels) +
  scale_color_manual(values=c("black","red")) +
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

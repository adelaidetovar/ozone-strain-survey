source("0_functions.R")
source("1_data_cleaning.R")

output_dir <- file.path("../../plots")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

set.wd(output_dir)

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

#############################
#### 1 ppm Ozone Plots ####
#############################

one_ppm <- read.csv('../../data/oneppm_pilot_data.csv')
one_ppm <- one_ppm[,c(1:3,5:8)]

me84_neutp <- ggplot(me84,aes(x=Strain,
                              y=Per_neu,
                              fill=Strain)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=Sex,stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(17,24)) +
  scale_fill_manual(values=me93_colors_sub) +
  scale_x_discrete(labels = me93_labels_sub) +
  scale_y_continuous(breaks=c(seq(0,20,by=5)),limits=c(0,20),labels=c(seq(0,20,by=5))) +
  theme_linedraw(base_size=14) +
  xlab(NULL) +
  labs(y="% neutrophils in BAL") +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.margin=unit(c(1,1,1,1),"lines"),
        legend.position="none")

t.test(log10(me84$Per_neu[1:6]),log10(me84$Per_neu[7:12])) # 0.001915

ggsave(plot=me84_neutp,'me84_percentneuts.png',dpi=300,
       units = "in", height = 4, width = 4)

me84_protein <- ggplot(me84,aes(x=Strain,
                                y=Protein_concentration,
                                fill=Strain)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=Sex,stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(17,24)) +
  scale_fill_manual(values=me93_colors_sub, guide = "none") +
  scale_x_discrete(labels = me93_labels_sub) +
  scale_y_continuous(breaks=c(seq(0,2600,by=500)),limits=c(0,2600)) +
  theme_linedraw(base_size=14) +
  xlab(NULL) +
  ylab("Concentration (ug/ml)")+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.margin=unit(c(1,1,1,1),"lines"),
        legend.position="none")

t.test(log(me84$Protein_concentration[1:6]),log(me84$Protein_concentration[7:12])) # 0.0003132

ggsave(plot=me84_protein,'me84_totalprotein.png',dpi=300,
       units = "in", height = 4, width = 4)

me84_total <- ggplot(me84,aes(x=Strain,
                              y=Total_cells,
                              fill=Strain)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=Sex,stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(17,24)) +
  scale_fill_manual(values=me93_colors_sub, guide = "none") +
  scale_x_discrete(labels = me93_labels_sub) +
  scale_y_continuous(breaks=c(seq(0,4e5,by=1e5)),labels=c(0,1,2,3,4),limits=c(0,4e5)) +
  theme_linedraw(base_size=14) +
  xlab(NULL) +
  labs(y=parse(text=paste0('"# of cells"','~ (10^5)'))) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.margin=unit(c(1,1,1,1),"lines"),
        legend.position="none")

t.test(log(me84$Total_cells[1:6]),log(me84$Total_cells[7:12])) # 0.003974


ggsave(plot=me84_total,'me84_totalcells.png',dpi=300,
       units = "in", height = 4, width = 4)

me84_neut <- ggplot(me84,aes(x=Strain,
                             y=Total_neu,
                             fill=Strain)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=Sex,stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(17,24)) +
  scale_fill_manual(values=me93_colors_sub, guide = "none") +
  scale_x_discrete(labels = me93_labels_sub) +
  scale_y_continuous(breaks=c(seq(0,4.5e4,by=1e4)),labels=c(0,1,2,3,4),limits=c(0,4.5e4)) +
  theme_linedraw(base_size=14) +
  xlab(NULL) +
  labs(y=parse(text=paste0('"# of neutrophils"','~ (10^4)'))) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.margin=unit(c(1,1,1,1),"lines"),
        legend.position="none")

t.test(log(me84$Total_neu[1:6]),log(me84$Total_neu[7:12])) # 0.001038

ggsave(plot=me84_neut,'me84_totalneuts.png',dpi=300,
       units = "in", height = 4, width = 4)


#BAL For samir

me93_neuts_labels = c("FA.CC017"="FA",
                      "O3.CC017"=expression("O"[3]),
                      "FA.CC003"="FA",
                      "O3.CC003"=expression("O"[3]))

sub_neuts$no_neu <- as.numeric(as.character(sub_neuts$no_neu))

sub_neuts$strain <- as.factor(sub_neuts$strain)
sub_neuts$strain <- factor(sub_neuts$strain, levels = c("CC003","CC017"))
me93_colors_neuts <- c( "#009E73","#56B4E9")
me93_neuts_strains <- c("CC003","CC017")

stat.test <- sub_neuts %>%
  group_by(strain) %>%
  t_test(no_neu ~ tx)

stat.test <- stat.test %>%
  add_xy_position(x = "strain", dodge = 0.8)

stat.test <- stat.test %>%
  mutate(xmax = c(2.2, 4.2),
         xmin = c(0.8, 2.8))

stat.test <- stat.test %>%
  mutate(p.adj = c("*","**"))

stat.test <- stat.test %>%
  mutate(y.position = c(657474,109444))

stat.test2 <- sub_neuts %>%
  filter(tx == "O3") %>%
  t_test(no_neu ~ strain)

stat.test2 <- stat.test2 %>%
  add_xy_position(x = "tx", dodge = 0.8)

stat.test2 <- stat.test2 %>%
  mutate(xmin = c(1.8),
         xmax = c(4.2),
         y.position = c(7.25e5))

stat.test2 <- stat.test2 %>%
  mutate(strain = c("CC017"))

stat.test2 <- stat.test2 %>%
  mutate(p.adj = c("*"))

me93_sub_neut <- ggplot(sub_neuts,aes(x=interaction(tx,strain),
                                      y=no_neu,
                                      fill=strain)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(shape=tx,stroke=1),height=0.0,width=0.2,size=3) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=me93_colors_sub, guide = "none") +
  theme_linedraw(base_size=14) +
  xlab(NULL) +
  labs(y=parse(text=paste0('"# of neutrophils"','~ (10^5)'))) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.margin=unit(c(1,1,2.75,1),"lines"),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_discrete(labels = me93_neuts_labels) +
  scale_y_continuous(breaks=c(seq(0,6.75e5,by=1e5)),labels=c(0,1,2,3,4,5,6)) +
  coord_cartesian(ylim=c(-1e4,7.75e5), clip = "off", expand = FALSE) +
  annotate(geom = "text", x = seq(1,4,by=1), y = -(0.06)*8.5e5, label = me93_neuts_labels, size = 4) +
  annotate(geom = "text", x = 1.5 + 2 * (0:1),y = -(0.105)*9e5, label = me93_neuts_strains, size= 6, fontface=2) +
  stat_pvalue_manual(
    stat.test, label = "p.adj", tip.length = 0,
    bracket.nudge.y = -1, label.size = 7
  ) +
  stat_pvalue_manual(
    stat.test2, label = "p.adj", tip.length = 0,
    bracket.nudge.y = -1, label.size = 7
  )


t.test(log(sub_neuts$no_neu[1,3,10,11,16,17]),log(sub_neuts$no_neu[7:12])) # 0.001038

ggsave(plot=me93_sub_neut,'me93_subset_neuts.png',dpi=300,
       units = "in", height = 4, width = 5)

#########################
#### ME75 pilot data ####
#########################

me93_colors_exp <- c("#56B4E9",rep("#000000",3),"#009E73","#F0E442",
                     "#0072B2",rep("#000000",5),"#D55E00")

me75_data <- read.csv("~/Downloads/me75.csv")
me75_data <- me75_data %>%
  group_by(Strain) %>%
  summarize(Prot = mean(Protein_concentration), Neut = mean(Per_neu))

me75_plot <- me75_data %>% ggplot(aes(x = Prot, y = Neut, fill = Strain)) +
  geom_jitter(shape = 21, aes(stroke=1),height=0.0,width=0.2,size=3) +
  theme_linedraw(base_size=14) +
  scale_fill_manual(values = me93_colors_exp) +
  labs(x = "Protein concentration (ug/mL)", y = "Percent neutrophils in BAL") +
  guides(fill = "none") +
  geom_text_repel(aes(label = Strain))

ggsave(me75_plot, filename = "For_resub/me75_plot.png", units = "in", height = 5, width = 5)

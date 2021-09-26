source("0_functions.R")
source("1_data_cleaning.R")

output_dir <- file.path("../output/pheno_data_stats/")

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

set.wd(output_dir)
rm(output_dir)

###############
#### ANOVA ####
###############

# Total cells LRT for interaction

base <- lm(log(total_cells) ~ strain + tx + sex, data = bal)
full <- lm(log(total_cells) ~ strain + tx + sex + strain:tx, data = bal)
sink("total_cells_LRT.txt")
anova(base,full,test="LRT")
sink()

# Total cells ANOVA

sink("total_cells_anova.txt")
anova(lm(log(total_cells) ~ strain + tx + sex + strain:tx, data = bal))
sink()


# Neutrophil percent and number LRT for interaction

base <- lm(log(no_neu+1) ~ strain + tx + sex, data = bal)
full <- lm(log(no_neu+1) ~ strain + tx + sex + strain:tx, data = bal)
sink("neut_number_LRT.txt")
anova(base,full,test="LRT")
sink()

base <- lm(log10(per_neu+1) ~ strain + tx + sex, data = bal)
full <- lm(log10(per_neu+1) ~ strain + tx + sex + strain:tx, data = bal)
sink("neut_percent_LRT.txt")
anova(base,full,test="LRT")
sink()

# Neutrophil percent and number ANOVA

sink("neut_number_anova.txt")
anova(lm(log(no_neu+1) ~ strain + tx + sex + strain:tx, data = bal))
sink()

sink("neut_percent_anova.txt")
anova(lm(log10(per_neu+1) ~ strain + tx + sex + strain:tx, data = bal))
sink()

# Macrophage percent and number LRT for interaction

base <- lm(log(no_macs+1) ~ strain + tx + sex, data = bal)
full <- lm(log(no_macs+1) ~ strain + tx + sex + strain:tx, data = bal)
sink("mac_number_LRT.txt")
anova(base,full,test="LRT")
sink()

base <- lm(log10(per_macs+1) ~ strain + tx + sex, data = bal)
full <- lm(log10(per_macs+1) ~ strain + tx + sex + strain:tx, data = bal)
sink("mac_percent_LRT.txt")
anova(base,full,test="LRT")
sink()

# Macrophage percent and number ANOVA

sink("mac_number_anova.txt")
anova(lm(log(no_macs+1) ~ strain + tx + sex + strain:tx, data = bal))
sink()

sink("mac_percent_anova.txt")
anova(lm(log10(per_macs+1) ~ strain + tx + sex + strain:tx, data = bal))
sink()

# Albumin ANOVA

sink("albumin_anova.txt")
anova(lm(sqrt(albumin_conc) ~ strain + tx + sex + strain:tx, data = albumin))
sink()

# Albumin LRT for interaction

base <- lm(sqrt(albumin_conc) ~ strain + tx + sex, data = albumin)
full <- lm(sqrt(albumin_conc) ~ strain + tx + sex + strain:tx, data = albumin)
sink("albumin_LRT.txt")
anova(base,full,test="LRT")
sink()

# Protein ANOVA

sink("protein_anova.txt")
anova(lm(log10(protein_conc) ~ strain + tx + sex + strain:tx, data = protein))
sink()

# Protein LRT for interaction 
base <- lm(log10(protein_conc) ~ strain + tx + sex, data = protein)
full <- lm(log10(protein_conc) ~ strain + tx + sex + strain:tx, data = protein)
sink("protein_LRT.txt")
anova(base, full, test="LRT")
sink()

############################
#### T-Tests for Values ####
############################

bal_t <- bal %>%
  group_by(strain) %>%
  mutate(norm_total = log(total_cells)) %>%
  rstatix::pairwise_t_test(norm_total ~ tx)

neut_t <- bal %>%
  group_by(strain) %>%
  mutate(neut_norm = log(no_neu+1)) %>%
  rstatix::pairwise_t_test(neut_norm ~ tx)

neup_t <- bal %>%
  group_by(strain) %>%
  mutate(neup_norm = log10(per_neu+1)) %>%
  rstatix::pairwise_t_test(neup_norm ~ tx)

mact_t <- bal %>%
  group_by(strain) %>%
  mutate(mact_norm = log(no_macs+1)) %>%
  rstatix::pairwise_t_test(mact_norm ~ tx)

macp_t <- bal %>%
  group_by(strain) %>%
  mutate(macp_norm = log10(per_macs+1)) %>%
  rstatix::pairwise_t_test(macp_norm ~ tx)

prot_t <- protein %>%
  group_by(strain) %>%
  mutate(protein_norm = log10(protein_conc)) %>%
  rstatix::pairwise_t_test(protein_norm ~ tx)

albumin_t <- albumin %>%
  group_by(strain) %>%
  mutate(albumin_norm = sqrt(albumin_conc),na.rm=TRUE) %>%
  rstatix::pairwise_t_test(albumin_norm ~ tx)

cytokines_melt <- melt(cytokines, id.vars=c("strain","tx"),
                       measure.vars = c("eotaxin","il10","il6","ip10","kc","lix","mip1b"))

cytokines_t <- cytokines_melt %>%
  group_by(strain,variable) %>%
  mutate(value_norm = log(value + 1)) %>%
  rstatix::pairwise_t_test(value_norm ~ tx)

cytokines_lrt <- cytokines_melt %>%
  group_by(variable) %>%
  mutate(value_norm = log(value + 1))

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="eotaxin",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="eotaxin",])
anova(base,full,test="LRT") # 1.06e-7

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="il10",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="il10",])
anova(base,full,test="LRT") # 0.00233

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="il6",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="il6",])
anova(base,full,test="LRT") # 0.000282

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="ip10",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="ip10",])
anova(base,full,test="LRT") # 0.1147

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="kc",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="kc",])
anova(base,full,test="LRT") # 0.5926

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="lix",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="lix",])
anova(base,full,test="LRT") # 0.5541

base <- lm(value_norm ~ strain + tx, data = cytokines_lrt[cytokines_lrt$variable=="mip1b",])
full <- lm(value_norm ~ strain + tx + strain:tx, data = cytokines_lrt[cytokines_lrt$variable=="mip1b",])
anova(base,full,test="LRT") # 0.1293

## one-sample t-tests for those with drop-outs
# eotaxin, CC0059
t.test(eotaxin$value[c(34,35,39,40)], mu=eotaxin$value[33])

# il6, CC003
t.test(il6$value[c(1,3:5,7)], mu=il6$value[2])

# mip1b, CC003
t.test(mip1b$mip1b[c(1,3:5,7)], mu=mip1b$mip1b[2])

# mip1b, CC017
t.test(mip1b$mip1b[c(11,16)],mu=mip1b$mip1b[14])

# mip1b, CC025
t.test(mip1b$mip1b[c(17,18,21)],mu=mip1b$mip1b[24])

# mip1b, CC039
t.test(mip1b$mip1b[c(25,26,30,31)], mu=mip1b$mip1b[29])
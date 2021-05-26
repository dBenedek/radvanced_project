#########################################################################################
#         Analysis of the genes location on the chromosome in case of given genes
#########################################################################################

## Load libraries and the datatable ----

library(tidyverse)
library(dplyr)
library(ggplot2)

GO_Terms <- read.table("D:/University/4_SEMESTER/R_advanced/miRNA_project/GO_terms_output_2.csv", header = TRUE)
head(GO_Terms)

View(GO_Terms)

#########################################################################################
#########################################################################################

# 1. Filtering for interesting genes: ----

# cholesterol homeostasis
chol_hom <- "cholesterol homeostasis"
chol_hom_go_terms <- GO_Terms %>% filter(name_1006 %in% chol_hom)
head(chol_hom_go_terms)

# positive regulation of translation - moderately general
pos_reg_translation <- "positive regulation of translation"
prt_go_terms <- GO_Terms %>% filter(name_1006 %in% pos_reg_translation)
head(prt_go_terms)

# gene silencing by miRNA - very general
gene_silencing_miRMA <- "gene silencing by miRNA"
silenc_go_terms <- GO_Terms %>% filter(name_1006 %in% gene_silencing_miRMA)

# positive regulation of MAPK cascade - very rare
pos_reg_MAPKc <- "positive regulation of MAPK cascade"
mapk_go_terms <- GO_Terms %>% filter(name_1006 %in% pos_reg_MAPKc)

# negative regulation of cell migration - general
neg_reg_migr <- "negative regulation of cell migration"
migr_go_terms <- GO_Terms %>% filter(name_1006 %in% neg_reg_migr)

# cellular response to insulin stimulus - rare
cresp_ins_stim <- "cellular response to insulin stimulus"
ins_stim_go_terms <- GO_Terms %>% filter(name_1006 %in% cresp_ins_stim)
ins_stim_go_terms

#########################################################################################
filter_go_terms <- c("GO:0042632", "GO:0045727", "GO:0035195", "GO:0043410", "GO:0030336", "GO:0032869")
filtered_go_terms <- GO_Terms %>% filter(go_id %in% filter_go_terms)
head(filtered_go_terms)

chr_occurrences_df <- as.data.frame(table(filtered_go_terms$go_id, filtered_go_terms$chromosome_name))
chr_occurrences

selected_genes_plot <- ggplot(chr_occurrences, aes(x = Var2, y = Freq, fill = Var1))+
  geom_col()+
  theme_light()+
  labs(x = "Chromosomes",
       y = "Number of occurences of the function's genes on the crhomosomes",
       title = "Distribution of molecular functions on the chromosomes")
selected_genes_plot
ggsave("D:/University/4_SEMESTER/R_advanced/miRNA_project/selected_genes_loc.png", selected_genes_plot, units = 'in', width = 8)

#########################################################################################
#########################################################################################

# Checking the occurences

# Cholesterol homeostasis regulating genes:
# Counting the occurences of the genes on the chromosomes:
chol_hom_occ <- table(chol_hom_go_terms$go_id, chol_hom_go_terms$chromosome_name)
chol_hom_occ

#chol_hom_occurrences <- as.data.frame(table(chol_hom_go_terms$go_id, chol_hom_go_terms$chromosome_name))
chol_hom_occurrences <- as.data.frame(chol_hom_occ)
head(chol_hom_occurrences)

chol_hom_occurrences <- arrange(chol_hom_occurrences, Var1)
head(chol_hom_occurrences)

ggplot(chol_hom_occurrences, aes(x = Var2, y = Freq))+
  geom_col()+
  theme_light()+
  labs(x = "Chromosome",
      y = "Number of occurences of the function's genes on the crhomosomes",
      title = "Cholesterol homeostasis regulating genes distribution in the genome")

#########################################################################################

# Genes that gives gene silencing by miRNA:
silenc_go_terms

mireg_occurrences <- table(silenc_go_terms$go_id, silenc_go_terms$chromosome_name)
mireg_occurrences
mireg_occurrences_df <- as.data.frame(mireg_occurrences)
head(mireg_occurrences_df)

ggplot(mireg_occurrences_df, aes(x = Var2, y = Freq))+
  geom_col()+
  theme_light()+
  labs(x = "Chromosome",
       y = "Number of occurences of the function's genes on the crhomosomes",
       title = "Distribution of the genes that gives gene silencing by miRNA")


#########################################################################################

# Genes that gives the positive regulation of translation:
prt_go_terms

prt_occurrences <- table(prt_go_terms$go_id, prt_go_terms$chromosome_name)
prt_occurrences
prt_occurrences_df <- as.data.frame(prt_occurrences)
head(prt_occurrences_df)

ggplot(prt_occurrences_df, aes(x = Var2, y = Freq))+
  geom_col()+
  theme_light()+
  labs(x = "Chromosome",
       y = "Number of occurences of the function's genes on the crhomosomes",
       title = "Translation regulating genes distribution in the genome")

#########################################################################################
# Genes that gives the cellular response to insulin stimulus:
ins_stim_go_terms 

ins_stim_occurrences <- table(ins_stim_go_terms$go_id, ins_stim_go_terms$chromosome_name)
ins_stim_occurrences
ins_stim_occurrences_df <- as.data.frame(ins_stim_occurrences)
head(ins_stim_occurrences_df)

ggplot(ins_stim_occurrences_df, aes(x = Var2, y = Freq))+
  geom_col()+
  theme_light()+
  labs(x = "Chromosome",
       y = "Number of occurences of the function's genes on the crhomosomes",
       title = "Distribution of the genes that gives response to insulin stimuli in the genome")

#########################################################################################
# Genes that gives the cellular response to insulin stimulus:
mapk_go_terms # - all of them located on the chr9

mapk_occurrences <- table(mapk_go_terms$go_id, mapk_go_terms$chromosome_name)
mapk_occurrences
mapk_occurrences_df <- as.data.frame(mapk_occurrences)
head(mapk_occurrences_df)

ggplot(mapk_occurrences_df, aes(x = Var2, y = Freq))+
  geom_col()+
  theme_light()+
  labs(x = "Chromosome",
       y = "Number of occurences of the function's genes on the crhomosomes",
       title = "Distribution of the genes that gives positive regulation of MAPK cascade")

#########################################################################################
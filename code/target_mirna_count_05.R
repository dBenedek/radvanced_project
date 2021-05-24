##############################################################
# Check miRNA-target gene interaction numbers per chromosome #
##############################################################

# 2021-04-15
# Benedek Dankó


# Load libraries ----
library("tidyverse")
library("tidyr")
library("data.table")
library("scales")
library("readxl")
require("rtracklayer")
library("patchwork")


# Load data & process data ----
setwd("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project")

# Load chromosome info:
chr_info <- fread("data/chromInfo.txt") %>% 
  dplyr::select(-V3)

gencode_annot <- readGFF("data/gencode.v37.annotation.gtf.gz", version=2L) %>% 
  mutate(gene_id = str_extract(gene_id, "ENSG\\d+"))

mirna_target_data <- read_tsv("data/Exp-miBRS_track_information_hg38.tsv") %>% 
  dplyr::select(c(MIRNA, GENES, CHROM)) %>% 
  separate_rows(MIRNA, sep = ",") %>%
  separate_rows(GENES, sep = ",") %>%
  filter(!is.na(CHROM)) %>% 
  mutate(MIRNA=str_replace(MIRNA, "miR", "mir")) %>% 
  mutate(MIRNA=str_replace(MIRNA, "hsa-", "")) %>% 
  mutate(MIRNA=str_replace(MIRNA, "-", "")) %>%
  mutate(MIRNA=str_replace(MIRNA, "mir", "MIR")) %>% 
  mutate(MIRNA=str_replace(MIRNA, "-\\dp", "")) %>% 
  mutate(MIRNA=str_replace(MIRNA, "let", "MIRLET")) %>% 
  mutate(MIRNA=toupper(MIRNA)) %>% 
  mutate(MIRNA=str_replace(MIRNA, "-", "")) %>% 
  left_join(gencode_annot, by=c("MIRNA"="gene_name")) %>% 
  dplyr::select(MIRNA:seqid) %>% 
  dplyr::rename(target_chr=CHROM,
                mirna_chr=seqid) %>% 
  filter(!is.na(mirna_chr)) %>% 
  filter(target_chr!="chrM")

plot1_data_genecount <- mirna_target_data %>% 
  dplyr::select(c(GENES, target_chr)) %>% 
  distinct() %>% 
  group_by(target_chr) %>% 
  summarize(count=n()) %>% 
  left_join(chr_info, by=c("target_chr"="V1")) %>% 
  mutate(norm_count=count/(V2/10000000))  # 10M
  

plot1_data_mirnacount <- mirna_target_data %>% 
  dplyr::select(c(MIRNA, mirna_chr)) %>% 
  distinct() %>% 
  group_by(mirna_chr) %>% 
  summarize(count=n()) %>% 
  left_join(chr_info, by=c("mirna_chr"="V1")) %>% 
  mutate(norm_count=count/(V2/1000000000)) # 10M
  

# plot1_data_mirnacount <- gencode_annot %>% 
#   filter(gene_type=="miRNA",
#          gene_name %in% mirna_target_data$MIRNA) %>% 
#   dplyr::select(c(seqid, gene_name)) %>% 
#   dplyr::rename(mirna_chr=seqid) %>% 
#   distinct() %>% 
#   group_by(mirna_chr) %>% 
#   summarize(count=n()) %>% 
#   left_join(chr_info, by=c("mirna_chr"="V1")) %>% 
#   mutate(norm_count=count/V2*1000000000) 


# Plot ----
cor <- cor(plot1_data_genecount$norm_count[1:23], plot1_data_mirnacount$norm_count)
level_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                 "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                 "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                 "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

barplot <- ggplot(plot1_data_genecount, 
                  aes(x=factor(target_chr, level=level_order), 
                      y=norm_count))+
  geom_bar(stat="identity", fill="steelblue")+
  geom_point(data=plot1_data_mirnacount, 
             aes(x=factor(mirna_chr),
                 y=norm_count),
             color="red",
             size=4)+
  scale_y_continuous(name = "Normalized target gene count",
                     sec.axis = sec_axis(~./100, 
                                         name="Normalized miRNA count"))+
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=10, hjust = 0.5, vjust = 0.5, angle=90),
        axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        axis.title.y.right = element_text(size=14,
                                          angle=90,
                                          margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  xlab("Chromosome")


# Save plot ----
ggsave("docs/barplot_mirna_target_count.png",
       barplot, units = "in", width = 10, height = 6)

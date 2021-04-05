##############################################################
# Check miRNA-target gene interactions based on target sites #
# (e.g. 5' UTR, promoter)                                    #
##############################################################

# 2021-04-05
# Benedek Dankó


# Load libraries ----
library("tidyverse")
library("tidyr")
library("data.table")


# Load data & process data ----
setwd("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project")
raw_data <- read.table("data/Exp-miBRS_track_information_hg38.tsv",
                       header=T) %>% 
  as.data.frame() %>% 
  dplyr::select(c(MIRNA, GENES, REGION, CHROM, REGION)) %>% 
  filter(!is.na("MIRNA")) %>% 
  separate_rows(MIRNA, sep=",") %>% 
  distinct() %>% 
  filter(!is.na("MIRNA") & MIRNA!="NA") %>% 
  dplyr::rename(TARGETG_CHR=CHROM) 

# write_tsv(x = raw_data,
#           file = "data/mirna_targetgene_pairs_v2.tsv")


# Process data ----
data_summary <- raw_data %>% 
  group_by(REGION) %>% 
  summarize(count=n())


# Plot ----
barplot <- ggplot(data_summary, 
                  aes(x=reorder(REGION, -count), y=count))+
  geom_bar(stat="identity")+
  geom_hline(aes(yintercept=100), color="red", linetype="dashed", alpha=.7)+
  scale_y_continuous(trans="log10")+
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=10, hjust = 0.5, vjust = 0.5, angle=90),
        axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  xlab("Target Gene Region")+
  ylab("Interaction count")


# Save plot ----
ggsave("docs/barplot_mirna_targetgene_by_region.png",
       barplot, units = "in", width = 8, height = 7)

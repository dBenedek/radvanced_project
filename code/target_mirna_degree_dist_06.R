# Load libraries ----
library("tidyverse")
library("tidyr")
library("data.table")
library("scales")
library("readxl")
require("rtracklayer")
library("patchwork")
library("igraph")
library("scales")


# Functions(s) ----
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  l <- gsub("%*%", "", l)
  l <- gsub("^+0", "^0", l)
  # return this as an expression
  parse(text=l)
}


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
  filter(!is.na(MIRNA) | !is.na(GENES) | !is.na(CHROM)) %>% 
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
  filter(!is.na(MIRNA),
         !is.na(mirna_chr)) %>% 
  distinct()


# Process data ----
mirna_graph <- graph_from_edgelist(as.matrix(mirna_target_data[,1:2]), 
                                   directed = F)
random_graph <- sample_gnm(n=length(c(mirna_target_data$MIRNA,   # number of vertices
                                      mirna_target_data$GENES)), 
                           m = nrow(mirna_target_data))          # number of edges

# Get the degree distribution data:
random_degree_dist <- data.frame(k=as.numeric(degree(random_graph))) %>% 
  group_by(k) %>% 
  summarize(freq=n()) %>% 
  mutate(freq=freq/sum(freq)) %>% 
  filter(k!=0) %>% 
  as.data.frame() 
mirna_degree_dist <- data.frame(k=as.numeric(degree(mirna_graph))) %>% 
  group_by(k) %>% 
  summarize(freq=n()) %>% 
  mutate(freq=freq/sum(freq)) %>% 
  as.data.frame()


# Plot ----
# Degree distribution of the miRNA - target mRNA network:
colors <- c("miRNA - target network" = "darkblue",
            "ER random network" = "green")
plot_01 <- ggplot(data=mirna_degree_dist, aes(x=k, y=freq))+
  geom_line(alpha=0.4, stat="identity", aes(color="miRNA - target network"))+
  geom_point(size=3, aes(color="miRNA - target network"))+
  geom_point(data=random_degree_dist, aes(x=k, y=freq, color="ER random network"),
             size=3)+
  geom_line(data=random_degree_dist, aes(x=k, y=freq, color="ER random network"), 
            alpha=0.4, stat="identity")+
  scale_x_continuous(trans="log10", labels=fancy_scientific)+
  scale_y_continuous(trans="log10", labels=fancy_scientific)+
  scale_color_manual(values = colors, name="")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, angle=0),
        axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size = 12))+
  xlab("k")+
  ylab("P(k)")


# Save plot ----
ggsave("docs/degree_dist.png",
       plot_01, units = "in", width = 10, height = 6)

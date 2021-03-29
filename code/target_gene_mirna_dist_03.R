#######################################################
# Check intrachromosomal miRNA - taget gene distances #
#######################################################


# Load libraries ----
library("tidyverse")
library("data.table")
library("scales")


# Load data ----
data <- read_tsv("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/mirna_targetgene_pairs.tsv")

# Load Gencode annotation, format:
gencode <- "/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/gencode.v37.annotation.gtf.gz"
gencode_annot <- fread(gencode,
                       col.names = c("chr", "ensemble", "type",
                                     "start", "end", "V6", "strand",
                                     "V8", "info")) %>% 
  mutate(gene_id = str_extract(info, "ENSG\\d+")) %>% 
  mutate(gene_name = str_extract(info, 'gene_name\\s.\\w+')) %>% 
  mutate(gene_name = str_replace(gene_name, 'gene_name \"', ""))


# Process data ----
data_annot <- data %>% 
  left_join(gencode_annot, by=c("target_gene_id"="gene_id")) %>% 
  filter(!is.na(target_gene_id) & !is.na(mirna_gene_id)) %>% 
  dplyr::rename(target_gene_start=start, # target gene start location
                target_gene_end=end) %>% # target gene end location
  dplyr::select(c(mirna, target_gene, mirna_gene_id, target_gene_id,
                  mirna_chr, target_gene_chr, target_gene_start,
                  target_gene_end)) %>% 
  left_join(gencode_annot, by=c("mirna_gene_id"="gene_id")) %>% 
  dplyr::rename(mirna_start=start,  # miRNA gene start location
                mirna_end=end) %>%  # miRNA gene end location
  filter(mirna_chr==target_gene_chr) %>% 
  mutate(start_dist=target_gene_start-mirna_start,
         end_dist=target_gene_end-mirna_end)


# Plot data ----
# Plot the distribution of intrachromosomal miRNA - target gene distances
plot_dist <- ggplot(data_annot, aes(x=start_dist))+
  geom_histogram(color="darkblue", fill="lightblue", bins=90)+
  scale_x_continuous(breaks = c(-200000000, 
                                -150000000,
                                -100000000,
                                -50000000,
                                0,
                                50000000,
                                100000000,
                                150000000,
                                200000000),
                     labels=c("-200 Mb", "-150 Mb",
                              "-100 Mb", "-50 Mb",
                              "0 Mb", "50 Mb",
                              "100 Mb", "150 Mb",
                              "200 Mb"))+
  xlab("Intrachromosomal distance of miRNA - target gene start sites")+
  ylab("Count")+
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, angle=0),
        axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size = 12))


# Save plot ----
ggsave("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/docs/intrachr_mirna_tagetgene_hist.png",
       plot_dist, units = "in", width = 10, height = 6)

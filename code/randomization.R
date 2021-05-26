####################################
# Randomization # by Sophie Farkas #
####################################

# Libraries ----

        library(tidyverse)
        # library("ComplexHeatmap")

# Data cleaning ----

        setwd("~/R")
        raw_data <- read.table("Exp-miBRS_track_information_hg38.tsv",header = TRUE) # Loading the source table, colum names from the first row
        miRNA_data <- subset(raw_data,subset=!is.na(MIRNA)) # Deleting NA rows
        miRNA_data_min <- miRNA_data[,c(1,10,12)] # Deleting useless information
        miRNA_data_min[miRNA_data_min=="chr1"] <- "chr01" #renameing the chromosomes for later ordering
        miRNA_data_min[miRNA_data_min=="chr2"] <- "chr02"
        miRNA_data_min[miRNA_data_min=="chr3"] <- "chr03"
        miRNA_data_min[miRNA_data_min=="chr4"] <- "chr04"
        miRNA_data_min[miRNA_data_min=="chr5"] <- "chr05"
        miRNA_data_min[miRNA_data_min=="chr6"] <- "chr06"
        miRNA_data_min[miRNA_data_min=="chr7"] <- "chr07"
        miRNA_data_min[miRNA_data_min=="chr8"] <- "chr08"
        miRNA_data_min[miRNA_data_min=="chr9"] <- "chr09"
        miRNA_data_3 <- separate_rows(miRNA_data_min,MIRNA,sep=",") # Separating the miRNAs to rows for each regulated gene
        miRNA_data_4 <- subset(miRNA_data_3,subset=MIRNA!="NA") # Deleting NA rows

# Randomizing ----

        miRNA_data_randomized <- miRNA_data_4 # Loading the original data to randomize
        miRNA_data_randomized$CHROM <- sample(unique(miRNA_data_4$CHROM),
                                              size=nrow(miRNA_data_4),
                                              replace=TRUE) # Randomizing the target gene's location
        # miRNA_data_randomized$GENES <- sample(miRNA_data_4$GENES,
        #                                       size=nrow(miRNA_data_4),
        #                                       replace=FALSE) # Randomizing the target gene
# Checking the occurrence ----

# original
        (occurence <- (table(miRNA_data_4$MIRNA,miRNA_data_4$CHROM)))
        occurence_df <- as.data.frame(table(miRNA_data_4$MIRNA,miRNA_data_4$CHROM))
        occurence_df <- arrange(occurence_df,Var1)
        occurence_df %>% 
                group_by(Var1) %>% 
                summarize(Var1,Var2,sum=sum(Freq),max=max(Freq),percent=(Freq)/sum(Freq))
        occ_percent <- occurence_df %>% 
                group_by(Var1) %>% 
                summarize(Var1,Var2,sum=sum(Freq),max=max(Freq),percent=(Freq)/sum(Freq))
        target100 <- occ_percent %>% 
                group_by(Var1) %>% 
                filter(sum>100) %>% 
                summarise(Var1,Var2,percent)

# random
        (occurence_random <- (table(miRNA_data_randomized$MIRNA,miRNA_data_randomized$CHROM)))
        occurence_random_df <- as.data.frame(table(miRNA_data_randomized$MIRNA,miRNA_data_randomized$CHROM))
        occurence_random_df <- arrange(occurence_random_df,Var1)
        occurence_random_df %>% 
                group_by(Var1) %>% 
                summarize(Var1,Var2,sum=sum(Freq),max=max(Freq),percent=(Freq)/sum(Freq))
        occ_rand_percent <- occurence_random_df %>% 
                group_by(Var1) %>% 
                summarize(Var1,Var2,sum=sum(Freq),max=max(Freq),percent=(Freq)/sum(Freq))
        target100_rand <- occ_rand_percent %>% 
                group_by(Var1) %>% 
                filter(sum>100) %>% 
                summarise(Var1,Var2,percent)

# Visualization ----

# original 1. plot
        png("original1.png", width = 700, height = 700)
        heatmap(occurence[,],
                scale="row",
                Rowv=TRUE,Colv=NA,
                xlab="Chromosome",ylab="miRNA",
                main="The number of miRNA target genes in each chromosome",
                margins=c(4,8))
        dev.off()

# random 1. plot
        png("random1.png", width = 700, height = 700)
        heatmap(occurence_random[],
                scale="row",
                Rowv=TRUE,Colv=NA,
                xlab="Chromosome",ylab="miRNA",
                main="The number of miRNA target genes in each chromosome \n after randomization",
                margins=c(4,8))
        # Heatmap(occurence_random[],
        #         name="The",
        #         cluster_columns = FALSE,
        #         clustering_method_rows ="complete",
        #         show_row_names = FALSE,
        #         row_title="miRNA",row_title_gp=gpar(fontsize = 14),
        #         column_title="Chromosome",column_title_gp=gpar(fontsize = 14))
        dev.off()

# original 2. plot
        original2 <- ggplot(occurence_df[,],aes(Var2,Var1,fill=Freq))+
                labs(x="Chromosome",
                     y="miRNA",
                     title="The occurence of target genes for miRNAs")+
                geom_tile()+
                theme(axis.title.x = element_text(size = 14),
                      axis.text.x = element_text(size = 12),
                      axis.title.y = element_text(size = 14),
                      axis.text.y = element_blank(),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12),
                      plot.title = element_text(size = 20))+
                scale_fill_gradient(name="Target \nnumber",
                                    low="white", high="darkblue",
                                    limits=c(0,100))
        ggsave("C:/Users/Lenovo/Documents/R/original2.png",
               original2, units = "in", width = 14, height = 10)

# random 2. plot
        random2 <- ggplot(occurence_random_df[,],aes(Var2,Var1,fill=Freq))+
                labs(x="Chromosome",
                     y="miRNA",
                     title="After randomization, the occurence of target genes for miRNAs")+
                geom_tile()+
                theme(axis.title.x = element_text(size = 14),
                      axis.text.x = element_text(size = 12),
                      axis.title.y = element_text(size = 14),
                      axis.text.y = element_blank(),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12),
                      plot.title = element_text(size = 20))+
                scale_fill_gradient(name="Target \nnumber",
                                    low="white", high="darkblue",
                                    limits=c(0,100))
        ggsave("C:/Users/Lenovo/Documents/R/random2.png",
               random2, units = "in", width = 14, height = 10)

# original 3. plot
        original3 <- ggplot(target100,aes(Var2,Var1,fill=percent))+
                labs(x="Chromosome",
                     y="miRNA",
                     title="The occurence of target genes for miRNAs, \nwhich have more than 100 targets")+
                geom_tile()+
                theme(axis.title.x = element_text(size = 14),
                      axis.text.x = element_text(size = 12),
                      axis.title.y = element_text(size = 14),
                      axis.text.y = element_text(size = 12),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12),
                      plot.title = element_text(size = 20))+
                scale_fill_gradient(name="Target \nfrequency",
                                    low="white", high="darkblue",
                                    limits=c(0,0.2))
        ggsave("C:/Users/Lenovo/Documents/R/original3.png",
               original3, units = "in", width = 15, height = 10)

# random 3. plot
        random3 <- ggplot(target100_rand,aes(Var2,Var1,fill=percent))+
                labs(x="Chromosome",
                     y="miRNA",
                     title="After randomization, the occurence of target genes for miRNAs, \nwhich have more than 100 targets")+
                geom_tile()+
                theme(axis.title.x = element_text(size = 14),
                      axis.text.x = element_text(size = 12),
                      axis.title.y = element_text(size = 14),
                      axis.text.y = element_text(size = 12),
                      legend.title = element_text(size = 14),
                      legend.text = element_text(size = 12),
                      plot.title = element_text(size = 20))+
                scale_fill_gradient(name="Target \nfrequency",
                                    low="white", high="darkblue",
                                    limits=c(0,0.2))
        ggsave("C:/Users/Lenovo/Documents/R/random3.png",
               random3, units = "in", width = 15, height = 10)
        
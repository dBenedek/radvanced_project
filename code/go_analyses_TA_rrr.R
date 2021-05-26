###################################################################################
#               A GO Term analyses using the biomaRt package
###################################################################################

library("biomaRt")
library("GO.db")

setwd("D:/University/4_SEMESTER/R_advanced/miRNA_project")

# Load the data: -----

mirnas <- read.table("D:/University/4_SEMESTER/R_advanced/miRNA_project/mirna_targetgene_pairs.tsv", header = TRUE)
head(mirnas)

# data cleaning is needed to remove the NA values: ----
gene_ids <- mirnas$mirna_gene_id

print(gene_ids)

gene_ids <- gene_ids[!is.na(gene_ids)]
print(gene_ids)

# //way2//Removing the NA values:
gene_ids <- na.omit(gene_ids)
head(gene_ids)

# GO analyses: ----

listMarts()
listDatasets(ensembl) # "hsapiens_gene_ensembl" 
listDatasets(ensembl)$dataset # just ..._gene_ensembl datasets printed out

# Choosing the ensembl mart and homo sapiens dataset:
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# Searching for filters (for the query which will be built)
filters = listFilters(ensembl)
head(filters)
# Check all the filters to find those what can be important for now:
View(filters)

# Check all the filters to find those what can be important for now:
attributes = listAttributes(ensembl)
head(attributes)
View(attributes)
# Selected attributes:
# ensembl_gene_id - Gene stable ID
# chromosome_name - Chromosome/scaffold name
# external_gene_name - Gene name
# transcript_count - Transcript count
# go_id - GO term accession
# name_1006 - GO term name

# Another way for writing out the attributes:
grep("go", attributes$name, ignore.case = T, value = T)
# Apply the same for filters:
grep("ensembl", filters$name, ignore.case = TRUE, value = TRUE)
# ensembl_gene_id - I will use this for filtering (I have these ensemble gene ids)

# Probes:
gene <- "ENSG00000188290"
geneSet <- c('ENSG00000188290', 'ENSG00000176022', 'ENSG00000176022')

# Probe on one gene:
getBM(attributes = c("ensembl_gene_id", "name_1006"), filters = "ensembl_gene_id", values = gene, mart = ensembl)
# Probe on one gene with all the attributes:
getBM(attributes = c("ensembl_gene_id", "name_1006", "go_id", "external_gene_name", "transcript_count", "chromosome_name"), filters = "ensembl_gene_id", values = gene, mart = ensembl)

# Probe on a set of genes:
# just 2 attributes
getBM(attributes =  c("ensembl_gene_id", "name_1006"), filters = "ensembl_gene_id", values = geneSet, mart = ensembl)
# all the chosen attributes:
getBM(attributes = c("ensembl_gene_id", "name_1006", "go_id", "external_gene_name", "transcript_count", "chromosome_name"), filters = "ensembl_gene_id", values = geneSet, mart = ensembl)


# Query on all the gene ids: 
GO_Terms_2_atrs <- getBM(attributes = "go_id", filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)
GO_Terms_2_atrs <- as.data.frame(GO_Terms_2_atrs)
write.table(GO_Terms_2_atrs, file = "GO_terms_2_atrs_output.csv", sep = "\t")

# Extended filtering - filtering data for GO enrichment analyses:
GO_Terms <- getBM(attributes = c("ensembl_gene_id", "name_1006", "go_id", "external_gene_name", "transcript_count", "chromosome_name"), filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)
GO_Terms <- as.data.frame(GO_Terms)
View(GO_Terms)

write.table(GO_Terms, file = "GO_terms_output_2.csv", sep = "\t")
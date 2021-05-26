setwd("D:/University/4_SEMESTER/R_advanced/miRNA_project")

# topGo tutorial on gene set enrichment analysis

## Quick start guide: ----

# Needed libraries:
library(topGO)
library(ALL)
# ALL = Acute Lymphoblastic Leukemia
data(ALL)
data("geneList")
head(ALL)

# a list of genes and the respective p-values for differential expression
# Data preparation

# needs the gene group, the GO terms and mapping that associate each gene with one or more GO terms:
# where to find GO annotations info is stored in the ALL object
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

# with the geneList, the topDiffGenes object is also loaded
# this function assumes that the provided argument (geneList) is a named vector of p-values
# this function give back the 50 genes with raw p-values less than 0.01 out of a total 323 genes
sum(topDiffGenes(geneList))

# Building the topGOdata object:
# - what is containing: 
#     the gene IDs and their scores,
#     the GO annotations,
#     the GO hierarchical structure,
#     and other needed informations
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10, #prune the GO hierarchy from the terms which have less than 10 annotated genes
                    annot = annFUN.db, affyLib = affyLib) # annFUN.db is used to extract the gene-to-GO mappings from the affyLib object
# Gives a summary about the object:
sampleGOdata
# this object facilitates the access to identifers, annotations and to basic data statistics
# This object is needed to the enrichment analyses.

# 3.2 Performing the enrichment tests ----

# Classical enrichment analyses:

# Fisherc exact test based on gene counts
# Kolmogorov-Smirnov test computes enrichment based on scores
# We can use both these tests since each gene has a score (representing how differentially expressed a gene is
# and by the means of topDiffGenes functions the genes are categorized into differentially expressed
# or not differentially expressed genes. --> stored in sampleGOdata 

# Fisher test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher

# Kolmogorov-Smirnov test - classic and elim method:

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

resultKS
resultKS.elim

# Analyses of results: ----

# GenTable function for analysing the go terms and p-values
# here:
# the top 10 significant GO terms is identified by the elim method
# + comparing the ranks and p-values of these terms (top10) obtained by classic method
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

# GenTable will return a dataframe with the top nodes
#includes some statistics on the GO terms and the p-values corresponding to each of the topGOresult object specifed as arguments
allRes


# score function to return the GO term's p-values from toGOresult
# in the example: differences between the classic and elim methods of Kolmogorov-Smirnov test:
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
# Scatter plot for visualization:
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
# colMap function:
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize, col = gCol)

# The plot also shows differences between the methods: classical methods found some GO terms for significant,
# but these were not significant by the elim method.

# Find and identify the number of annotated genes for the GO terms:
sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
        elim = pValue.elim[sel.go],
        classic = pValue.classic[sel.go])
# Go praph to visualize the data:
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')


## 4. Loading genes and annotations data: ----

# 4.1 Getting started

# libraries: topGO, ALL, (data - ALL)
# When the topGO package is loaded three environments GOBPTerm, GOMFTerm and GOCCTerm are created and bound to the package environment.
# GOBPTerm function split the GOTERM environment into three different ontologies:
# BP - biological function, MF - molecular function, CC - cellular component
BPterms <- ls(GOBPTerm)
head(BPterms)
BPterms

# Usually one needs to remove probes/genes with low expression value as well as probes with very small variability across samples.
# "genefilter" package is for  filtering genes
library(genefilter)
selProbes <- genefilter(ALL, filterfun(pOverA(0.20, log2(100)), function(x) (IQR(x) > 0.25)))
eset <- ALL[selProbes, ]
# This filter selects 4101 probesets out of 12625 probesets available on the hgu95av2 array.

# 4.2 The topGOdata object ----

# Needed datas for it:
# - gene ids + gene-wise scores
#     -->scores can cam from: t-test statistic or the p-value for differential expression, correlation with a phenotype, or any other relevant score
# - mapping between gene identifers and GO terms
# - GO hierarchical structure - from GO.db package

# 4.3 Custom annotations ----
# How to build topGOdata object from custom GO annotations:
# gene-to-GO or GO-to-gene mappings should be provided
# example (examples/geneid2go.map) contains gene-to-GOs
#   For each gene IDs are listed the GO terms to which gene is specifically annotated
# readMappings function for parse this file
geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
str(head(geneID2GO))

# List of character vectors: list names - gene IDs
# character vectors contains the GO ids annotated to the specific gene
# readMappings function also able to read mappings from a text file.

# inverseList function creates the inverse of the mappings
# from gene-to-GOs to GO-to-genes (or vice-versa)
GO2geneID <- inverseList(geneID2GO)
View(GO2geneID)
str(head(GO2geneID))

# 4.4 Predefined list of interesting genes: ----
# test the enrichment of GO terms with regard to this list of interesting genes
# In this scenario, when only a list of interesting genes is provided,
# the user can use only tests statistics that are based on gene counts,
# like Fisher's exact test, Z score and alike

geneNames <- names(geneID2GO)
head(geneNames)
# select 10% of the genes, randomly (store in 'myInterestingGenes')
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
# geneList factor indicates that genes are interesting or not
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

# Building GO data:
# with Molecular Function ontology
# mapping given by geneID2GO list, it is used for annFUN.gene2GO function
# annFUN function: to compile a list of GO terms such that each element in the list is a character vector
# containing all the gene identifiers that are mapped to the respective GO term
GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata

#### On my data ####
goterm_data <- read.table("D:/University/4_SEMESTER/R_advanced/miRNA_project/GO_terms_output_2.csv", header = TRUE)
head(goterm_data)
got_gene_go_ids <- goterm_data$go_id
names(got_gene_go_ids) <- goterm_data$ensembl_gene_id
got_gene_go_ids
geneNames <- names(goterm_data$ensembl_gene_id)
geneNames
names(goterm_data$ensembl_gene_id)
# Here I got this "NULL" and I just don't succeed to repair it.
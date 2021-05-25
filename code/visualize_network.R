##############################################################################
#                       Analyze miRNA - target networks                     #
##############################################################################

# 2021-03-24
# Luca Csabai

install.packages("igraph")
if (!requireNamespace("RCy3", quietly = TRUE)) 
  BiocManager::install("RCy3")

# Load libraries ----
library("tidyverse")
library("data.table")
library("igraph")
library("RCy3")



# Read data ----
network_raw <- read_tsv("input/Exp-miBRS_track_information_hg38.tsv")

# Functions from Benedek----
create_pairs_1 <- function(x){
  mirnas <- strsplit(x[1], ",") %>% 
    `[[`(1)
  target_gene <- x[2]
  df <- data.frame(mirnas=mirnas,
                   target_gene=rep(target_gene, length(mirnas)))
  return(df)
}

create_pairs_2 <- function(x){
  genes <- strsplit(x[1], ",") %>% 
    `[[`(1)
  mirna <- x[2]
  df <- data.frame(target_gene=genes,
                   mirna=rep(mirna, length(genes)))
  return(df)
}

# Process data ----
# Select the miRNA(s) and target gene(s) interactions only exonic:
network_raw <- network_raw[network_raw$REGION == 'exonic', ] 

int_data <- network_raw %>% 
  dplyr::select(c(MIRNA, GENES)) %>% 
  filter(!is.na(MIRNA))

# Process data, to have single, pairwise target gene - miRNA interactions:
int_df <- apply(int_data, 1, create_pairs_1)  
int_df <- do.call("rbind", int_df) %>% 
  dplyr::select(target_gene, mirnas)
int_df <- apply(int_df, 1, create_pairs_2)  
int_df <- do.call("rbind", int_df)

# Remove NA values that were in a list 
int_df[ int_df == "NA" ] <- NA
int_df <- int_df %>% drop_na()

# Select onmy first 50 interactions for visualization purposes
int_df <- int_df[1:70,]

# Separate nodes
target_gene <- int_df %>%
  distinct(target_gene) %>%
  rename(label = target_gene)

source_mirna <- int_df %>%
  distinct(mirna) %>%
  rename(label = mirna)


node_df <- full_join(source_mirna, target_gene, by = "label")


# Change direction 
int_df <- int_df[,c(2, 1)]

# Create an igraph object
net <- graph.data.frame(int_df[1:70,], directed=T)
# Visualize with igraph
#V(net)$size <- 10

#l <- layout.kamada.kawai(net)
#l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

#plot(net, rescale=F, layout=l*2, vertex.label=NA)

calc_network_stats <- function(network){
  # Calculate network density
  # The proportion of present edges from all possible edges in the network.
  density_res <- edge_density(network, loops=F)
  
  # Calculate transitivity
  transitivity_res <- triad_census(network)
  
  # Calculate degree (both in and out)
  degree_res <- degree(network, mode="all")
  # Add it to a data frame
  degree_res2 <- data.frame(degree=sort(degree_res, decreasing=TRUE)) %>% tibble::rownames_to_column()
  nodes_2 <- left_join(node_df, degree_res2, by = c("label"="rowname"))
  
  # Calculate betweenness centrality
  # centrality based on a broker position connecting others
  betweenness_res <- betweenness(network, v = V(network), directed = FALSE, weights = NULL,
                                 nobigint = TRUE, normalized = FALSE)
  # Add it to a data frame
  betweenness_res2 <- data.frame(betweenness_centrality=sort(betweenness_res, decreasing=TRUE)) %>% tibble::rownames_to_column()
  nodes_3 <- left_join(nodes_2, betweenness_res2, by = c("label"="rowname"))
  
  # Calculate hub score
  hub_res <- hub_score(network, scale = TRUE, weights = NULL,
                       options = arpack_defaults)
  # Add it to a data frame
  hub_res2 <- data.frame(hub_score=sort(hub_res$vector, decreasing=TRUE)) %>% tibble::rownames_to_column()
  nodes_4 <- left_join(nodes_3, hub_res2, by = c("label"="rowname"))
  
  return (nodes_4)
}

nodes_stat_df <- calc_network_stats(net)

nodes_stat_df

# Set color to be blue for miRNAs
nodes_stat_df$color <- ifelse(grepl("hsa", nodes_stat_df$label), "lightblue", "orange")


# this ensures the starting random position is the same
# for the layouts that use a random starting position
set.seed(1492) 

l <- layout_with_fr(net)

l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)


# Plot node sizes according to node degree
plot(net, rescale=F, layout=l*1.0,
       edge.arrow.size=0.5, 
       vertex.label=NA,
       vertex.shape="circle", 
       vertex.size=nodes_stat_df$degree*3, 
       vertex.label.color="black", 
       edge.width=0.5,
       vertex.color= nodes_stat_df$color,
       main = 'Degrees')
  
# Plot node sizes according to hubness
plot(net, rescale=F, layout=l*1.0,
       edge.arrow.size=0.5, 
       vertex.label=NA,
       vertex.shape="circle", 
       vertex.size=nodes_stat_df$hub_score*10, 
       vertex.label.color="black", 
       edge.width=0.5,
       vertex.color= nodes_stat_df$color,
       main = 'Hubs')
dev.off()


##### Import into Cytoscape ##### ----

# Preprocess
net2 <- int_df %>% dplyr::rename("source"= "mirna", "target"="target_gene")
net2['interaction'] = 'inhibits'
nodes2 <- nodes_stat_df %>% dplyr::rename("id" = "label")

# Determine if betweenness centrality calcualtions have been carried out
if("betweenness_centrality" %in% colnames(nodes_stat_df)){
  bet = TRUE
} else {
  bet=FALSE 
}

# See if cytoscape is open
cyto_error <- FALSE
tryCatch( {msg <- cytoscapePing () } , error = function(e) {cyto_error <<- TRUE})

if (!cyto_error){
  continue = TRUE
  message('Successfully connected to Cytoscape - continuing with visualisation')
} else {
  continue = FALSE
  message('Could not connect to Cytoscape - skipping visualisation')
}

# Run visualisation if cytoscape open
if(continue){
  createNetworkFromDataFrames(nodes2, edges=net2[,1:2], title="exonic network")
  edges <- net2 %>%
    mutate(key=paste(source, "(interacts with)", target))
  print(edges)
  loadTableData(edges[,3:4], data.key.column = 'key', table = 'edge')
  
  # Colour by Betweenness centrality if data exists in network
  if (bet){
    setNodeColorMapping("betweenness_centrality",mapping.type="c", network="exonic network")
  }
  
  # Save cys file
  saveSession(filename = file.path("output/exonic_network.cys"))
}



#### Lema-Datura microbiome multilayer network ####

##### Packages #####

# For download the required R packages, use the following code

# mlBioNets
#library("devtools")
#devtools::install_github("Nertekkad/mlBioNets")

# phyloseq
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("phyloseq")

# minet
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("minet")

# SpiecEasi
#devtools::install_github("zdk123/SpiecEasi")

# igraph
#install.packages("igraph")

# muxViz
#devtools::install_github("manlius/muxViz")


###### Load packages and data #####

library("mlBioNets")
library("phyloseq")
library("minet")
library("SpiecEasi")
library("igraph")
library("muxViz")

# Data from Mayoral et al. (2022)
data("ps2")

# Otus and taxa tables
T_table <- as.data.frame(tax_table(ps2))
O_table <- as.data.frame(t(otu_table(ps2)))

# Let's collapse the tables at genus level
T_Collapsed <- T_collapse(T_table = T_table, O_table = O_table, 
                          names_level = "Genus")

# Separate insect and plant data
Insect <- which(sample_data(ps2)$Type =="Insect")
Insect <- sample_data(ps2)$ID[Insect]
Plant <- which(sample_data(ps2)$Type =="Plant")
Plant <- sample_data(ps2)$ID[Plant]
Insectmat <- T_Collapsed[Insect,]
Plantmat <- T_Collapsed[Plant,]
# Node colors at phylum level
unq <- unique(T_table[,"Phylum"])
colors <- rainbow(length(unq))

##### ARACNe algorithm #####

# Co-abundance networks

# Beetle
mim <- build.mim(Insectmat,estimator="spearman")
aracne_mat <- aracne(mim)
Insect_minet <- graph.adjacency(aracne_mat)
Insect_minet <- as.undirected(Insect_minet)

# Solanaceous
mim <- build.mim(Plantmat,estimator="spearman")
aracne_mat <- aracne(mim)
Plant_minet <- graph.adjacency(aracne_mat)
Plant_minet <- as.undirected(Plant_minet)

# Color assignment for the nodes of the layers
g.list <- list(Insect_minet, Plant_minet)
g.list <- v_colored_ml(g.list, T_table, g_tax = "Phylum",
                       p_tax = "Genus", g_colors = colors)

###### Plot multilayer network by phylum ######
abs_IP <- list(Insectmat, Plantmat)
lay <- layoutMultiplex(g.list, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(g.list, layer.layout=lay,
                 layer.colors=c("blue3", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_IP, g.list, 10),
                 node.colors=node_color_mat(g.list, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

###### Plot multilayer network by centrality ######
g.list<-ctr_ml(g.list, ctr_type = "degree")
lay <- layoutMultiplex(g.list, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(g.list, layer.layout=lay,
                 layer.colors=c("blue3", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_IP, g.list, 20),
                 node.colors=node_color_mat(g.list, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)


#### SparCC algorithm ####

# Coabundance networks

# Insect (beetle)
sparccNet<-sparcc(Insectmat)
sparccNet <- abs(sparccNet$Cor) >= 0.4
insect_sparCC<-adj2igraph(sparccNet)
vertex.attributes(insect_sparCC) <- list(name = colnames(Insectmat))

# Plant (solanaceous)
sparccNet<-sparcc(Plantmat)
sparccNet <- abs(sparccNet$Cor) >= 0.4
plant_sparCC<-adj2igraph(sparccNet)
vertex.attributes(plant_sparCC) <- list(name = colnames(Plantmat))

# Color assignment for the nodes of the layers
g.list2 <- list(insect_sparCC, plant_sparCC)
g.list2 <- v_colored_ml(g.list2, T_table, g_tax = "Phylum",
                       p_tax = "Genus", g_colors = colors)

###### Plot multilayer network by phylum ######
abs_IP<-list(Insectmat, Plantmat)
lay <- layoutMultiplex(g.list, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(g.list2, layer.layout=lay,
                 layer.colors=c("red3", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_IP, g.list2, 20),
                 node.colors=node_color_mat(g.list2, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)


###### Plot multilayer network by centrality ######
g.list2<-ctr_ml(g.list2, ctr_type = "degree")
abs_IP<-list(Insectmat, Plantmat)
lay <- layoutMultiplex(g.list2, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(g.list2, layer.layout=lay,
                 layer.colors=c("blue3", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_IP, g.list2, 20),
                 node.colors=node_color_mat(g.list2, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)


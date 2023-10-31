#===========================================================================
############################ Network analysis ##############################
############################   [Workshop]  #################################
#===========================================================================
# R version 4.0.2 (2020-06-22)
#===========================================================================
#========== Customizing the directory of our workshop =====================
#===========================================================================

# 0) Create a new file
# 1) Create a directory called NetworkAnalysis in the right-bottom window and 
#    get in it.
# 2) Set this directory as our working directory running the next command:

setwd("~/Desktop/NetworkAnalysis")  # Ctrl + enter
getwd()

# Windows:
# setwd("c:/docs/mydir") 
# Linux and mac:
# setwd("/usr/sol/mydir") 

# 3) Save the script in the current location with the name networks_workshop.R   
# 4) Inside NetworkAnalysis create a directory called data
# 5) Inside NetworkAnalysis create another directory called results

#===========================================================================
#================ Install R packages if not installed =====================
#===========================================================================

if (!require("igraph")) install.packages("igraph")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("tidyr")) install.packages("tidyr")

#===========================================================================
#========================== Attach R packages ==============================
#===========================================================================

library(igraph)
library(RColorBrewer)
library(tidyr)

#===========================================================================
#========================= First example ===================================
#===========================================================================

set.seed(25032022)
g <- sample_gnp(15, 0.2, directed = FALSE, loops = FALSE)
plot(g)

plot(g,
     vertex.color="darkorchid",  # nodes colors
     vertex.size=20,             # size of nodes
     edge.color="black")         # edges colors

V(g)
E(g)
degree(g)
distances(g)
is_connected(g)
count_components(g)
components(g)
all_shortest_paths(g, 1, to = 5)
average.path.length(g, directed=FALSE, unconnected=TRUE)
diameter(g)
# ED: the portion of the potential connections in a network that are actual connections
# Calculated: how many edges there are in a network divided by the total possible number of edges
edge_density(g)  
# BC: is a measure of centrality in a graph based on shortest paths
# A node with higher betweenness centrality would have more control over the network, because more information will pass through that node
betweenness(g) 

write.graph(g,"~/Desktop/fist_net.graphml", format="graphml")

set.seed(25032022)

png("fist_net.png", width = 300*10, height = 300*8, res = 300, units = "px")

plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color="gray80")            # color de las aristas

dev.off()

#===========================================================================
#=========================== Study Case ====================================
#===========================================================================
# Data (11 Kb)
#===========================================================================
#=============================== Get the data ==============================
#===========================================================================

download.file(
  url = "https://github.com/biocomplexnet/NetworkAnalysis/blob/main/data/coexp_h_5.tsv",
  destfile = "data/coexp_h_5.tsv"
)

download.file(
  url = "https://github.com/biocomplexnet/NetworkAnalysis/blob/main/data/coexp_t_5.tsv",
  destfile = "data/coexp_t_5.tsv"
)

download.file(
  url = "https://github.com/biocomplexnet/NetworkAnalysis/blob/main/data/treatment_1.tsv",
  destfile = "data/treatment_1.tsv"
)

download.file(
  url = "https://github.com/biocomplexnet/NetworkAnalysis/blob/main/data/treatment_2.tsv",
  destfile = "data/treatment_2.tsv"
)

#===========================================================================
#=============================== Manage the data ===========================
#===========================================================================

# For this part of the workshop we are going to use: treatment_1.tsv file,
# but you can use anyone of the databases that we downloaded.

df <- read.csv("data/treatment_1.tsv", sep = "\t",header = TRUE)
View(df)

# Manage the dataset:

B <- as.data.frame(table(df)) # Create an edge weight column named "Freq"
B1 <- subset(B,Freq>0) # Delete all the edges having weight equal to 0

#===========================================================================
#=======================  Create an igraph object  =========================
#===========================================================================

net <- graph.data.frame(B1, directed = FALSE)
View(net)

# Display vertices and edges:
V(net)
E(net)

# Assigning weight to each edge
E(net)$weight<-E(net)$Freq
net

#==================================================================#
#===================  Exploring the Network  =====================#
#==================================================================#

# size and order of the network:
gsize(net)   # Number of edges
gorder(net)    # Number of vertices

# Vertices:
V(net)

# Edges:
E(net)


# Adjacency matrix
order <- gorder(net)
net[c(1:order),c(1:order)]

#==================================================================#
#===================  Studying the centrality  =====================#
#==================================================================#

# Assigning degree to the vertices

net_degree <- degree(net, mode = c("All"))
View(net_degree)
V(net)$degree <- net_degree
V(net)$degree

which.max(net_degree)

# Agregar camino más corto 
all_shortest_paths(net, 1, to = 5)
average.path.length(net, directed=FALSE, unconnected=TRUE)

# Betweenness centrality
net_bw <-betweenness(net, directed = FALSE)
V(net)$betweenness<-net_bw
V(net)$betweenness
which.max(net_bw)

# Eigenvector centrality
net_eig <- evcent(net)$vector
V(net)$Eigen<-net_eig
V(net)$Eigen
which.max(net_eig)

DF<-as_long_data_frame(net)
DF  # Explicar este dataframe


#==================================================================#
#===================  Network global structure  =====================#
#==================================================================#

# Network density
edge_density(net) # Global density 

# Conectivity
is_connected(net)
count_components(net)
components(net)

# Shortest and average path

all_shortest_paths(net, 1, to = 5)
average.path.length(net, directed=FALSE, unconnected=TRUE)

# Distance matrix

distances(net)


#==================================================================#
#===================== Network Visualization ======================#
#==================================================================#


# Plotting a network with the vertices size based on degree centrality

set.seed(1001) 

pal<-brewer.pal(length(unique(V(net)$name)), "Set3") # Vertex color assigned per each class number
plot(net,
     edge.color = 'black',
     vertex.label.cex =0.1,
     vertex.color=pal[as.numeric(as.factor(vertex_attr(net, "name")))],
     vertex.size = sqrt(net_degree)/3, 
     edge.width=sqrt(E(net)$weight/10),
     layout = layout.fruchterman.reingold
)


# Plotting a network with the vertices size based on eigenvector centrality

set.seed(1001)
plot(net,
     edge.color = 'black',
     vertex.label.cex =0.1,
     vertex.color=pal[as.numeric(as.factor(vertex_attr(net, "name")))],
     vertex.size = sqrt(net_eig)*10,
     edge.width=sqrt(E(net)$weight/10),
     layout = layout.fruchterman.reingold
)

# Plotting a network with the vertices size based on betweenness centrality
set.seed(1001)
plot(net,
     edge.color = 'black',
     vertex.label.cex =0.1,
     vertex.color=pal[as.numeric(as.factor(vertex_attr(net, "name")))],
     vertex.size = sqrt(net_bw)/3, edge.width=sqrt(E(net)$weight/10),
     layout = layout.fruchterman.reingold
)

## Plotting a scatter plot to see the correlation

#3.1. Between degree and betweenness centrality

plot(V(net)$degree, V(net)$betweenness)
cor(V(net)$degree, V(net)$betweenness)

#3.2. Between degree and eigenvector centrality

plot(V(net)$degree, V(net)$Eigen)
cor(V(net)$degree, V(net)$Eigen)


#==================================================================#
#====================== Community Detection =======================#
#==================================================================#


# Louvain clustering

lc <- cluster_louvain(net) # Create a cluster based on the Louvain method
communities(lc) # You can check which vertices belongs to which clusters.

#2. Plotting and saving figure of Betweenness Centrality network with the community detection

set.seed(1001) # To duplicate the computer process and create exactly the same network repetitively you should set the seed.

#png("results/community_detection.png", width = 300*10, height = 300*8, res = 300, units = "px")

plot(lc, 
     net,
     edge.color = 'black',
     vertex.label.cex =0.1,
     vertex.color=pal[as.numeric(as.factor(vertex_attr(net, "Class")))],
     vertex.size = sqrt(net_degree)/3, edge.width=sqrt(E(net)$weight/10),
     layout = layout.fruchterman.reingold)

#dev.off()

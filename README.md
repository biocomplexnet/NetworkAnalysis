# RiaBio Wokshop 2023: Bioinformatics for biodata analysis using Machine Learning and Network approaches

## Co-expression networks
#### Maribel Hernandez Rosales, CINVESTAV, Mexico
#### Marisol Navarro Miranda, CINVESTAV, Mexico
#### César Alfonso Díaz Mijangos, CINVESTAV, Mexico

Lets look at [igraph funtions/tools](https://igraph.org/r/html/latest/).

``` r
setwd("~/Desktop/")

# install.packages('igraph')
library(igraph)         # Network analysis library

```

### We generate a random network of 15 undirected nodes. We set a seed for random numbers generation.

``` r
set.seed(25032022)
g <- sample_gnp(15, 0.2, directed = FALSE, loops = FALSE)
plot(g)

plot(g,
     vertex.color="darkorchid",  # nodes colors
     vertex.size=20,             # size of nodes
     edge.color="black")         # edges colors
```

# We are going to use igraph functions to obtain some metrics:

Definition of the metrics can be foun in [Network Science](http://networksciencebook.com/) from Albert-László Barabási  

# 1) Vertexes

``` r
V(g)
```

# 2) Edges

``` r
E(g)
```

# 3) Degree

``` r
degree(g)

```
# 4) Distance matrix, paths

``` r
distances(g)
```

# 5) Connected components

``` r
is_connected(g)
count_components(g)
components(g)
```

# 6) Shortest path

``` r
all_shortest_paths(g, 1, to = 5)
average.path.length(g, directed=FALSE, unconnected=TRUE)
```

# 7) Diameter

``` r
diameter(g)
```

# 8) Density

``` r
edge_density(g)
```

# 9) Betweenness centrality

``` r
betweenness(g)
```

# Write our network as graphml

``` r
write.graph(g,"/Users/solouli/Desktop/yeast_protein_interaction.graphml", format="graphml")
```

# Save our network as png

``` r
png("random_network.png", width = 300*10, height = 300*8,
    res = 300, units = "px")

set.seed(25032022)

plot(g,
     vertex.size=2*sqrt(V(g)$degree),
     vertex.color=V(g)$color,
     vertex.label=NA,
     layout=layout_nicely,
     edge.color="gray80")            # color de las aristas

dev.off()
```

# Case Study: Co-expression networks from "Deciphering the Tissue-Specific Regulatory Role of Intronless Genes Across Cancers"

A **gene co-expression network (GCN)** is an undirected graph, where each node corresponds to a gene, and a pair of nodes is connected with an edge if there is a significant co-expression relationship between them.

**Intronless genes (IGs) or single-exon genes lacking introns** constitute approximately 3% of the human genome. IGs also often display tissue-specific expression. These characteristics translate into IG-associated diseases, mainly neuropathies, developmental disorders, and cancer.

We are going to work with subsets from the work [Deciphering the Tissue-Specific Regulatory Role of Intronless Genes Across Cancers](hhttps://link.springer.com/chapter/10.1007/978-3-031-06220-9_18).

## References

Stuart, Joshua M; Segal, Eran; Koller, Daphne; Kim, Stuart K (2003). "A gene-coexpression network for global discovery of conserved genetic modules". Science. 302 (5643): 249–55.

Grzybowska EA. Human intronless genes: functional groups, associated diseases, evolution, and mRNA processing in absence of splicing. Biochem Biophys Res Commun. 2012 Jul 20;424(1):1-6. doi: 10.1016/j.bbrc.2012.06.092. Epub 2012 Jun 23. PMID: 22732409.

Aviña-Padilla, K. et al. (2022). Deciphering the Tissue-Specific Regulatory Role of Intronless Genes Across Cancers. In: Jin, L., Durand, D. (eds) Comparative Genomics. RECOMB-CG 2022. Lecture Notes in Computer Science(), vol 13234. Springer, Cham. https://doi.org/10.1007/978-3-031-06220-9_18

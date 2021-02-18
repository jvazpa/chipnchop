################################################################
#                                                              #
#                                                              #
#                  TRANSCRIPTIONAL NETWORKS                    #
#                                                              #
#                                                              #
#  Author: Joaquín Tamargo Azpilicueta (joatamazp@alum.us.es)  #
#          José Vázquez Pacheco                                #
#          Ana González Toro                                   #
#                                                              #
#                   University of Seville                      #
#                                                              #
################################################################

## For elaboration of networks, igraph library is installed

library(igraph)

## tf_targets file stores all the cistrome genes of each 20 transcription
## factor. The names of all of these transcription factors is written at
## targets.txt file included in tf_targets directory. Workspace is thus
## chosen to be tf_targets file.

## 1. Reading transcription factor targets

## 1.1. Reading transcription factors names file (at its first column):

tf_names<-read.table(file = "targets.txt", header = FALSE, as.is = TRUE)[[1]]

## 1.2. Reading genes affected by those transcription factors:

genes_in_network <- c()

for(i in 1:length(tf_names))
{
  genes_in_network<-c(genes_in_network,read.table(file = paste0(tf_names[i], ".txt"), header = F, as.is = TRUE)[[1]])
}

genes_in_network<-unique(genes_in_network)
length(genes_in_network)


## 2. Generation of adjacency matrix

## 2.1. Creating an empty square matrix with as much rows and columns as unique genes
## are present in the cistrome

adjacency_matrix <- matrix(0,nrow = length(genes_in_network), ncol = length(genes_in_network))
colnames(adjacency_matrix) <- genes_in_network
rownames(adjacency_matrix) <- genes_in_network

## 2.2. Filling in the matrices with cistrome information

for(i in 1:length(tf_names))
{
  current_tf_target<-read.table(file = paste0(tf_names[i],".txt"), header = F, as.is = TRUE)[[1]]
  adjacency_matrix[tf_names[i],current_tf_target] <- adjacency_matrix[tf_names[i],current_tf_target]+1
}

## 3. Creation of graphs from adjacency matrices

transcriptional_network <- graph_from_adjacency_matrix(adjmatrix = adjacency_matrix, mode = "directed")

tf_subgraph <- induced_subgraph(graph = transcriptional_network, vids = tf_names)

## Writing gml graphs for visualization

write_graph(graph = tf_subgraph, file = "tf_subgraph.gml",format = "gml")


## 4. Network characterization

## 4.1. Node degrees distribution

## Let's check if this network is scale-free. Firstly, we will represent a histogram
## of the node degree distribution:

hist(x = degree(tf_subgraph, mode="total"), xlab = "Node degree", ylab = "")

## 4.2. Statistical analysis

## 4.2.1. Kolmogorov-Smirnoff test

## We can statistically analyze if the
## distribution suits a negative potential (Kolmogorov-Smirnoff test):

## H0 (Null hypothesis): It follows a negative potential distribution
## H1 (Alternative hypothesis): It doesn't follow a negative potential distribution

## Nevertheless, analyzing this sort of datasets with so little nodes may 
## have important errors due to the nature of the test: it doesn't follow regular
## null/alternative hypothesis structure where null hypothesis is the one
## tested to be wrong. Eventually, this leads to loss of control on the errors
## derived from this approach.

## 4.2.2. Linear regression

## Studying the linear regression of the negative potential (base e):

node_degree_total<-degree(tf_subgraph,mode="total")

node_degree_freq<-table(node_degree_total)

plot(log(as.numeric(names(node_degree_freq))), log(node_degree_freq),pch=19,cex=1)

lm.res <- lm(log(node.degree.freq) ~ log(as.numeric(names(node.degree.freq))))

summary(lm.res)

## When r^2 is really low and p-value is high, it means that
## the null hypothesis (consisting of not sticking
## to a negative potential linear regression), is not discarded.
## 
## To sum up, low r^2 and high p-values may indicate that the network
## doesn't follow a negative potential.

## 4.2.3. Statistical analysis comparing with other randomly generated networks 

## 4.2.3.1. Check the degree of each node

## We will need to build networks that have nodes with same degree but different targets.

node.degree.out <- degree(graph = tf_subgraph, mode = "out")

## 4.2.3.2. For loop that creates and evaluate randomly generated networks

number.randomisations <- 5000
random.autoreg <- vector(mode = "numeric", length = number.randomisations)
for (j in 1:number.randomisations)
{
  ## For each iteration, create a new blank adjacency network
  
  random_adjacency_network <- matrix(0, nrow = length(node.degree.out), ncol = length(node.degree.out))
  print(j)
  
  i=1
  pvalue=0
  for (i in 1:length(node.degree.out))
  {
    ## For each transcription factor, randomly choose any of the transcription factors and
    ## add 1 to that position in the matrix. This assignation would mean "TF A regulates
    ## TF B expression".
    
    random_adjacency_network[i,sample(x = 1:21,size = node.degree.out[i],replace = FALSE)] <- 1
  }
  
  ## When autorregulation exist, diagonal "1"s will appear in the adjacency matrix.
  ## We will sum up all of the diagonal values and store them in a vector.
  
  random.autoreg[j] <- sum(diag(random_adjacency_network))
}

hist(random.autoreg)

## 4.2.3.3. Comparison with the original adjacency matrix

## Adjacency matrix is extracted from the subgraph obtained in appendix 3 and diagonal sum is done.

tf_adjacency_matrix <- as.matrix(as_adjacency_matrix(graph = tf_subgraph))

sum(diag(tf_adjacency_matrix))

## The comparison will be done between all of the values contained in random.autoreg (the vector 
## containing all the diagonal sums for each iteration) and the sum of the diagonal values of the
## original adjacency matrix. 

sum(random.autoreg > sum(diag(tf_adjacency_matrix)))

## p-value is calculated on the basis of that calculation, dividing by all the iterations carried out:

sum(random.autoreg > sum(diag(tf_adjacency_matrix)))/number.randomisations

## As this value is really low (0, indeed), it can be established that AUTORREGULATION IS A 
## NETWORK MOTIF.


## 5. Network motifs analysis

## igraph has an implemented function that allows the generation of all possible isocratic
## networks given a fixed number of nodes.

help("graph.isocreate")

## For 3 nodes, there are 16 possible isocratic networks:

plot.igraph(graph.isocreate(size = 3, directed = T, number = 0))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 1))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 2))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 3))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 4))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 5))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 6))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 7))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 8))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 9))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 10))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 11))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 12))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 13))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 14))
plot.igraph(graph.isocreate(size = 3, directed = T, number = 15))

## To determine which of these networks constitute a network motif, a function 
## implemented in igraph allows to recognize and count how many times each of 
## these subgraphs appear in a given network.

help("graph.motifs")

graph.motifs(graph = tf_subgraph,size = 3)

## graph.motif function returns a vector that include, for each isomorphic 
## network, the number of times they appear. This means the vector will
## be as long as the possible number of isomorphic networks. For 3 nodes
## networks, this would mean a length of 16.

length(graph.motifs(graph = tf_subgraph,size = 3))

## To determine whether a motif is significantly overrepresented in respect
## to others or not, it is mandatory to build an statistical method.
## Recurring to the previously described method of randomly generating
## networks and including the graph.motif function:

number.randomisations <- 5000
motifs.matrix <- matrix(0, ncol = length(graph.motifs(graph = tf_subgraph,size = 3)), nrow = number.randomisations)
random.autoreg <- vector(mode = "numeric", length = number.randomisations)
for (j in 1:number.randomisations)
{
  random_adjacency_network <- matrix(0, nrow = length(node.degree.out), ncol = length(node.degree.out))
  print(j)
  
  i=1
  pvalue=0
  
  for (i in 1:length(node.degree.out))
  {
   random_adjacency_network[i,sample(x = 1:21,size = node.degree.out[i],replace = FALSE)] <- 1
  }
  
  random.autoreg[j] <- sum(diag(random_adjacency_network))
  
  random.graph <- graph.adjacency(adjmatrix = random_adjacency_network, mode = "directed")
  
  motifs.matrix[j,] <- graph.motifs(graph = random.graph, size = 3)
}

## In the original network, the motif identification resulted in:

tfs_motifs <- graph.motifs(graph = tf_subgraph,size = 3)

## Statistical analysis must identify if a given motif is overrepresented
## in respect to others. Thus, we will count how many more times is the motif
## represented in the random networks than on the original. This value
## will be normalized using the number of iterations.

sum(motifs.matrix[,3] > tfs_motifs[3])/number.randomisations

## It results in a p-value ~ 0,2 > 0,05, what means that it is not a network
## motif in the original network.


## A loop might be useful for looking into each of the possible motifs:

pvalues_motifs <- vector(length = length(graph.motifs(graph = tf_subgraph,size = 3)),mode = "numeric")
for (k in 1:length(graph.motifs(graph = tf_subgraph,size = 3)))
{
    pvalues_motifs[k] <-(sum(motifs.matrix[,k] > tfs_motifs[k]))/number.randomisations
}

pvalues_motifs

## There is an implemented function at igraph that allows to show
## all the instances of a given network for the specified motif size.

help("graph.subisomorphic.lad")

graph.subisomorphic.lad(graph.isocreate(size = 3,directed = T,number = 14), tf_subgraph, all.maps = T)[["maps"]]

## Frenquently, there will be feedback loops with multiple exits that.




---
bibliography: references.bib
link-citations: true
output:
  pdf_document: default
  html_document: default
---

## Introduction

## Methods

Communicability quantify how well information transit between two nodes by considering all possible path in a network and penalizing longer ones. It is compute with the exponent of the adjacency matrix of the network.$$G=\sum_{k=0}^{\infty} \frac{\left(\mathbf{A}^k\right)}{k !}=e^{\mathbf{A}}$$

where $G$ is the communicability matrix, $A$ the adjacency matrix and k is used as a penalizing term. It is possible to compute the exponential of a matrix with the graph spectrum :

$$
G=\sum^n_{j=1}\varphi_j\varphi_j^T e^{\lambda_j}
$$

where $\varphi_j$ and $\lambda_j$ are respectively the $j^{th}$ eigenvectors and eigenvalues of the matrix $A$. The spectral form of $G$ can be decompose by the following way :

$$
G=\varphi_1\varphi_1^{T} e^{\lambda_1}+ 
\sum^n_{j=2}\varphi_j^{+}\varphi_j^{+T} e^{\lambda_j}+ 
\sum^n_{j=2}\varphi_j^{-}\varphi_j^{-T} e^{\lambda_j}+
\sum^n_{j=2}\varphi_j^{-}\varphi_j^{+T} e^{\lambda_j}
$${#eq:cluster1}

\
\

where $\varphi^+$ or $\varphi^-$ indicate respectively all the positives or negatives values of the $j^{th}$ eigenvector. The first is not include because all the values the eigenvector the same sign. @Estrada2008Communicability explain that "two nodes have the same sign in an eigenvector if they can be considered as being in the same partition of the network, while those pairs having different signs correspond to nodes in different partitions." So to make partition we are mostly interested by the sign of the sums in @eq:cluster1 : $$\sum^{intracluster}_{j=2}\varphi_j\varphi_j^{T} e^{\lambda_j} = \sum^n_{j=2}\varphi_j^{+}\varphi_j^{+T} e^{\lambda_j}+ \sum^n_{j=2}\varphi_j^{-}\varphi_j^{-T} e^{\lambda_j}$$

and

$$\sum^{intercluster}_{j=2}\varphi_j\varphi_j^{T} e^{\lambda_j} = \sum^n_{j=2}\varphi_j^{-}\varphi_j^{+T} e^{\lambda_j}$$

so the clustering matrix is obtain with

$$\Delta G = \sum^{intracluster}_{j=2}\varphi_j\varphi_j^{T} e^{\lambda_j} -
\left|\sum^{intercluster}_{j=2}\varphi_j\varphi_j^{T} e^{\lambda_j}\right|
$$

in short it is in fact

$$
\Delta G = 
G - \varphi_{1} \varphi_{1}^T e^{\lambda_1}
$$

#### Example

```{r}
library(tidyverse)
library(lattice)
library(igraph)
library(colorRamps)
A  = matrix(c(0,1,0,1,1,0,0,0,0,0,0,
              1,0,1,1,1,0,0,0,0,0,0,
              0,1,0,1,1,0,0,0,0,0,0,
              1,1,1,0,1,0,0,0,0,0,0,
              1,1,1,1,0,1,0,0,0,0,0,
              0,0,0,0,1,0,1,0,0,0,0,
              0,0,0,0,0,1,0,1,1,0,1,
              0,0,0,0,0,0,1,0,1,1,1,
              0,0,0,0,0,0,1,1,0,1,1,
              0,0,0,0,0,0,0,1,1,0,1,
              0,0,0,0,0,0,1,1,1,1,0), nrow =11, ncol =11)
grap = graph_from_adjacency_matrix(A, mode = "undirected")
plot(grap)
```

A graph with 11 nodes and 2 distinct group. First we need to compute the graph spectrum

```{r}
spectra = eigen(A)
levelplot(spectra$vectors, ylab = "eigenvectors", xlab ="j th position")
```

Now let's take the $2^{nd}$ dimension as an example.

```{r}
##
G_dim2 = spectra$vectors[,2]%*%t(spectra$vectors[,2])*exp(spectra$values[2])
intra = G_dim2[G_dim2]
levelplot(G_dim2, ylab = "node", xlab ="node",col.regions = rev(matlab.like(16)))
```

And that it ! The second dimension of the graph communicability identify 2 cluster (blue).

We can compute for the third dimention

```{r}
G_dim3 = spectra$vectors[,3]%*%t(spectra$vectors[,3])*exp(spectra$values[3])
levelplot(G_dim3, ylab = "node", xlab ="node",col.regions = rev(matlab.like(16)))
```

Which identify clusters between 5:6 and 6:7. The cluster of the third dimension are less "obvious" than those from the second dimension

Now if we want to use other communicability dimension we just have to add

```{r}
levelplot(G_dim2+G_dim3, ylab = "node", xlab ="node",col.regions = rev(matlab.like(16)))
```

We could continue like that till the last dimension (11th), but it was for the explanation. So now we can compute directly $\Delta G$ by adding all the 10 dimension

```{r}
delta_G = matrix(0, nrow =11, ncol =11)
for(dim in 2:11){
  delta_G = delta_G + spectra$vectors[,dim]%*%t(spectra$vectors[,dim])*exp(spectra$values[dim])
}
levelplot(delta_G, ylab = "node", xlab ="node",col.regions = rev(matlab.like(16)))
```

## Results

![Figure 1 : Global matrix of clustering communicability. Positive values indicate species in same cluster, negative value species in "opposite" cluster. Only host order and virus order names are display on x and y.](figures/host_virus_cluster.jpg){#fig:global}

![Figure 2](figures/host_host_cluster.jpg){#fig:host}

![Figure 3](figures/host_host_recap.jpg){#fig:host_recap}

![Figure 4](figures/virus_virus_cluster.jpg){#fig:virus}

![Figure 5](figures/virus_virus_recap.jpg){#fig:virus_recap}

![Figure 4](figures/host_virus.jpg){#fig:host_virus}

![Figure 5](figures/host_virus_recap.jpg){#fig:host_virus_recap}

## Conclusion

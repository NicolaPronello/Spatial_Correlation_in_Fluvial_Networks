# ------------------------------------------------------------------------------
# This script allows to analyze spatial correlation patterns of water quality 
# variables in fluvial networks. In particular, we consider a toy example with
# simulated data
#
# Paper: Analyzing Spatial Correlation Patterns of Water Quality Variables in Fluvial Networks
# Authors: N. Pronello, S. Castiglia, V. Frontuto, N. Golini, R. Ignaccolo, L. Ippoliti
# Journal: Journal of Agricultural, Biological, and Environmental Statistics 
# Year: 2026
# DOI: 10.1007/s13253-026-00726-9
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Load the necessary R libraries
# ------------------------------------------------------------------------------
library(igraph)
library(ggspatial)
library(gridExtra)
library(ggnetwork)
library(e1071)
library(mgcv)
library(ggpubr)
library(dendextend)
library(corpcor)
library(dplyr)

# ------------------------------------------------------------------------------
# Load the ad-hoc functions to implement this analysis
# ------------------------------------------------------------------------------

# Set the working directory to the directory containing Functions4covAnalysis.R”
source("../R/Functions4covAnalysis.R")


# ------------------------------------------------------------------------------
# Graphical representation of the selected domain
# ------------------------------------------------------------------------------

# Set the working directory to the directory containing the data
setwd("../data/")

Adj_matrix <- read.table("adjacency_matrix.csv") # Adjacency matrix among water bodies
coord_points <- read.table("coord_points.csv") # Matrix of the coordinates of water bodies centroids
D <- read.table("Domain.csv") # Dataframe to plot in ggplot() the fluvial network
ID_wb <- rownames(Adj_matrix) # Water bodies ID
L <- diag(apply(Adj_matrix,1,sum))-Adj_matrix # Laplacian matrix
# res <- eigen(L) Warning: eigen() may produce different results depending 
                  # on the software version and the operating system (macOS or Windows)  
res <- readRDS("eigen_decom.rds") # the ouput of eigen() applied to the Laplacian matrix
thre <- -.02 # Threshold used to select a subdomain
nv <- length(which(res$vectors[,536]>thre)) # Number of points in sub-domain
id_rid <- which(res$vectors[,536]>thre) # id points in sub-domain
clrs <- rep("lightgrey",dim(coord_points)[1])
clrs[id_rid] <- "violet"

# Network representation where vertices correspond to water bodies. 
# In violet the selected sites which constitute the sub-network utilized as new domain.
g1 <- ggplot() +
  geom_segment( data=D, aes( x=x.from, y=y.from, xend=x.to, yend=y.to ),lineend="round",col="lightgray")+
  geom_point(data=coord_points[-id_rid,], aes( X,Y),size=2,color = scales::alpha("violet", 0.8),shape=21,fill = "white",stroke=.5)+
  geom_point(data=coord_points[id_rid,], aes( X,Y),size=2,shape = 21, fill = scales::alpha("violet", 0.6), color = "violet", stroke = .5) +
  theme_blank()+theme_bw()+theme_void()
g1

# Routine to obtain a dataframe D1 useful to plot in ggplot() the sub-domain as a network
D1 <- c()
G <- Adj_matrix[id_rid,id_rid]
diag(G) <- 0
value <- apply(G,1,function(x){sum(x!=0)})
for(h in 1:dim(G)[1]){
  if(value[h] != 0){
    for(j in 1:value[h]){
      D1 <- rbind(D1, c(ID_wb[id_rid][h],ID_wb[id_rid][which(G[h,]!=0)[j]],G[h,G[h,]!=0][j],value[h]))
    }
  }
}

D1 <- as.data.frame(D1)
colnames(D1) <- c("from","to","factor","value")
class(D1$factor) <- "numeric"
class(D1$value) <- "numeric"
D1$x.from <- rep(0,dim(D1)[1])
D1$y.from <- rep(0,dim(D1)[1])
D1$x.to <- rep(0,dim(D1)[1])
D1$y.to <- rep(0,dim(D1)[1])
layout <- coord_points[id_rid,]
rownames(layout) <- ID_wb[id_rid]
for(i in 1:dim(D1)[1]){
  id.from <- which(D1$from[i]==rownames(layout))
  id.to <- which(D1$to[i]==rownames(layout))
  D1$x.from[i] <- layout[id.from,1]
  D1$y.from[i] <- layout[id.from,2]
  D1$x.to[i] <- layout[id.to,1]
  D1$y.to[i] <- layout[id.to,2]
}


# Network representation of the sub-domain, where vertices correspond to water bodies. 
g2 <- ggplot() +
  geom_segment( data=D1, aes( x=x.from, y=y.from, xend=x.to, yend=y.to ),lineend="round",col="lightgray")+
  geom_point(data=coord_points[id_rid,], aes( X,Y),size=2,shape = 21, fill = scales::alpha("violet", 0.6), color = "violet", stroke = .5) +
  theme_blank()+theme_bw()+theme_void()
g2

ggarrange(g1,g2,ncol=2)
 
# -------------------------------------------------------------------------------------
# Estimation of correlation matrices with just one observation in some sample locations 
# -------------------------------------------------------------------------------------
L <- diag(apply(Adj_matrix[id_rid,id_rid],1,sum))-Adj_matrix[id_rid,id_rid] # Laplacian matrix of the sub-domain
true_covs <- readRDS("true_covs.rds") # True covariance matrices 5x5 obtained as described in the paper
set.seed(12)
id_sample <- sample(1:nv,25,replace=F) # Choose 25 sites at random over the domain

# Network representation of the sub-domain, where vertices in violet correspond to selected water bodies.
ggplot() +
  geom_segment( data=D1, aes( x=x.from, y=y.from, xend=x.to, yend=y.to ),lineend="round",col="lightgrey")+
  geom_point(data=coord_points[id_rid,], aes( X,Y),size=2,color = scales::alpha("violet", 0.8),shape=21,fill = "white",stroke=.5)+
  geom_point(data=coord_points[id_rid[id_sample],], aes( X,Y),size=2,shape = 21, fill = scales::alpha("violet", 0.6), color = "violet", stroke = .5)+
  theme_blank()+theme_bw()+theme_void()+ggtitle("")


# Data simulation: data simulated where no mean trend is assumed
nvo <- length(id_sample) # Number of selected vertices
cov_star <- list()
C_estim <- list()
p <- dim(true_covs[[1]])[1]

for(i in 1:nvo){
     X <- rmvn(1,rep(0,p),true_covs[[id_sample[i]]])
     cov_star[[i]] <- X%*%t(X)
}

gamma_values <- seq(.1,20,by=1)
res_cv <- cv_fun_cov_holes(nv,gamma_values,L,cov_star,id_sample) # Cross-validation to select gamma 
gamma <- gamma_values[which.min(res_cv[,2])]
K <- solve(diag(nv)+gamma*L) # Kernel matrix for all sub-domain
K_r_rid <- matrix(0,nv,nvo) # Kernel matrix for nvo observed sites
for (j in 1:nv) {
  weights <- K[j,id_sample] 
  K_r_rid[j,] <- weights/sum(weights)
}

for(i in 1:nv){ # Interpolation on all nv vertices of the network
  R <- Reduce("+", Map(function(M, w) M * w, cov_star, K_r_rid[i, ]))
  C_estim[[i]] <- R
}

# -------------------------------------------------------------------------------------
# Clustering Betti numbers from estimated partial correlation matrices 
# Be careful: the simulation may not be consistent with the clustering structure.
# -------------------------------------------------------------------------------------
curves_betti <- matrix(0,ncol=p,nrow=nv)
for(i in 1:nv){ # Calculate Betti_0 curves 
  P <- solve(C_estim[[i]]) # Precision matrix
  cor <- (diag(sqrt(1/diag(P)))%*%P%*%diag(sqrt(1/diag(P))))*(-1) # Diagonal will be removed in the next function
  result <- betti_0_filtration_mst(1-round(cor,8))
  curves_betti[i,] <- result$eta
}

ks_matrix <- matrix(0, nrow = nv, ncol = nv)
for (i in 1:(nv - 1)) {
  for (j in (i + 1):nv) {
    d <- ks_distance(curves_betti[i, ], curves_betti[j, ])  # Calculate the KS-distance between two Betti_0 curves 
    ks_matrix[i, j] <- d
    ks_matrix[j, i] <- d 
  }
}

dend_oss <- hclust(as.dist(ks_matrix), method = "average") # Hierarchical clustering of all sites using the average linkage method
k_opt <- 4 # choose the number of clusters
cls <- cutree(dend_oss, k = k_opt) # assigns each water bodie to a cluster 
dend <- as.dendrogram(dend_oss)
pal <- c("dodgerblue","darkorange","#1C1C1C","#D40078")
d_oss <- color_branches(
  dend, 
  k = 4,
  col = pal[unique(cls[labels(dend)])]
)
par(mar = c(3,3, 1, 1))
par(cex.axis = 2)
d_oss <- set(d_oss, "branches_lwd", 2)
plot(d_oss, leaflab = "none")

# Spatial clustering map, where water bodies are colored based on their cluster assignment
ggplot() +
  geom_segment( data=D1, aes( x=x.from, y=y.from, xend=x.to, yend=y.to ),lineend="round",col="lightgray")+
  geom_point(data=coord_points[id_rid,], aes( X,Y),size=4,shape = 21, fill = pal[cls], color = pal[cls], stroke = .5) +
  theme_blank()+theme_bw()+theme_void()




# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')


## Real Data Application: Detecting disease related SNPs.

# Compute the z-scores and distance matrix
prepare_data <- function(){
  # read the ld info file
  lddata<-read.delim("./data/SNPs/myld.ld",sep="")

  #calculate the z-scores
  z_snp<-testgenom$Beta/testgenom$SE

  #calculating the distance matrix
  d_snp<-matrix(rep(1, m*m), m, m)
  x_ld<-rep(0,length(lddata$SNP_A))
  y_ld<-rep(0,length(lddata$SNP_B))

  for(i in 1:length(lddata$CHR_A)){
    cat(i)
    x_ld[i]<-which(testgenom$Pos==lddata$BP_A[i])
    y_ld[i]<-which(testgenom$Pos==lddata$BP_B[i])
    d_snp[x_ld[i],y_ld[i]]=1-lddata$R2[i]
    d_snp[y_ld[i],x_ld[i]]=d_snp[x_ld[i],y_ld[i]]
  }
}

data_stored = TRUE
if (!data_stored) {
  prepare_data()
}else{
  load("./data/SNPs/z_snp.RData")
  load("./data/SNPs/d_snp.RData")
}

# read the summary stats for the selected 5000 SNPs
testgenom<-read.delim("./data/SNPs/S5000.txt")
rs_snp<-read.delim("./data/SNPs/IDs.txt", header = FALSE)

m<-5000
mu0<-rep(0,m)
sd0<-rep(1,m)
fdr_level<-c(0.001,0.01,0.05,0.1)
nrep<-length(fdr_level)

pv<-2*pnorm(-abs(z_snp), 0, 1)

n_cluster.vec <- c(5, 10, 20)
n_cluster<-length(n_cluster.vec)

bon.nr<-rep(0,nrep)
bh.nr<-rep(0,nrep)
lasla.dd.nr<-rep(0,nrep)
adapt.nr_list<-replicate(n_cluster, rep(0,nrep), simplify = FALSE)

bon.de<-matrix(rep(0, m*nrep), nrep, m)
bh.de<-matrix(rep(0, m*nrep), nrep, m)
lasla.dd.de<-matrix(rep(0, m*nrep), nrep, m)
adapt.de_list<-replicate(n_cluster, matrix(rep(0, m*nrep), nrep, m), simplify = FALSE)


for (i in 1:nrep) {
  cat("Running with FDR level:", fdr_level[i], "\n")
  q<-fdr_level[i]
  
  bh.res<-bh.func(pv, q)

  bw=density(z_snp,bw="SJ-ste")$bw
  pis_lasla <- lasla_pis(z_snp, d_snp/1.4, pv, tau=bh.func(pv,0.8)$th, h=bw, eps=0)

  weight <- lasla_weights(z_snp, d_snp/1.4, pis_lasla, mu0, sd0, q, h=bw, eps=0)
  lasla.dd.res <- lasla_thres(pv, pis_lasla, weight, q)
  
  d <- as.dist(d_snp)
  hc <- hclust(d, method = "complete")
  for(k in seq_along(n_cluster.vec)) {
    n_cluster = n_cluster.vec[k]
    # Generate the clusters from the distance matrix
    clusters <- cutree(hc, k = n_cluster)
    
    # Generate the indicators
    s<-matrix(0,m,n_cluster)
    for(cluster_i in 1:n_cluster){
      s[, cluster_i] <- ifelse(clusters == cluster_i, 1, 0)
    }
    adapt.res <-  adapt_xgboost(s,pv,
                                  verbose = list(print = FALSE, 
                                                 fit = FALSE,
                                                 ms = FALSE),
                                  piargs = list("nrounds" = 50,
                                                "max_depth" = 1,
                                                "nthread" = 1,
                                                "verbose" = 0),
                                  muargs = list("nrounds" = 50,
                                                "max_depth" = 1,
                                                "nthread" = 1,
                                                "verbose" = 0),
                                  alphas = c(q),
                                  nfits = 5)
    adapt.de_list[[k]][i, which(adapt.res$qvals <= q)]<-1
    adapt.nr_list[[k]][i]<-length(which(adapt.res$qvals <= q))
  }
    
  bon.nr[i]<-length(which(pv<=q/m))
  bh.nr[i]<-bh.res$nr
  lasla.dd.nr[i]<-lasla.dd.res$nr
  
  bon.de[i,which(pv<=q/m)]<-1
  bh.de[i,]<-bh.res$de
  lasla.dd.de[i,]<-lasla.dd.res$de
}




#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

rejections <- data.frame()
nr = cbind(bon.nr, bh.nr,lasla.dd.nr)

adapt.names <- character(0)
for(k in 1:length(n_cluster.vec)) {
  method_name <- sprintf("adapt_cluster_%d", n_cluster.vec[k])
  adapt.names <- c(adapt.names, method_name)
  
  nr <- cbind(nr, adapt.nr_list[[k]])
  colnames(nr)[ncol(nr)] <- paste(method_name, ".nr", sep = "")
}

method_names <- c("Bonferroni", "BH", "LASLA", adapt.names)

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(fdr_level)),
                    FDR = fdr_level,
                    nr = nr[,i])
  rejections<- rbind(rejections, tmp)
}
save(rejections, file=sprintf("%s/SNPs.RData", data_dir))



#######################################
#          Plot networks             #
#######################################
if (FALSE){
  install.packages("network")
  install.packages('ergm')
}

library(network)
library(ergm)

# under FDR=0.05, analysis the rejections made by BH and LASLA
i = 3
bh.rej.idx <- which(bh.de[i,] == 1)
bh.rej.snp <- testgenom[bh.rej.idx,]
lasla.rej.idx <- which(lasla.dd.de[i,] == 1)
lasla.rej.snp <- testgenom[lasla.rej.idx,]

# Verify that rejections by BH is a subset of which by LASLA
if (all(bh.rej.idx %in% lasla.rej.idx)){
  cat("BH rejctions at FDR level ", format(fdr_level[i], nsmall=2), 
        "is a subset of LASLA's. \n")
}

# Find SNPs that are rejected by LASLA but not by BH
extra.rej.idx <- setdiff(lasla.rej.idx, bh.rej.idx)
extra.rej.snp <- testgenom[extra.rej.idx,]


# Full net rejected by LASLA
n <- length(lasla.rej.idx)
fullmat <- matrix(rep(0, n*n), n,n)
for (i in 1:n) {
  for (j in i:n) {
    fullmat[i,j] = 1-d_snp[lasla.rej.idx[i], lasla.rej.idx[j]]
    fullmat[j,i] = fullmat[i,j]
  }
}
fullmat <- ceiling(fullmat)
diag(fullmat) <- 0    ## to make sure there is no self-edge
fullnet <- as.network(x=fullmat,
                   directed = FALSE,
                   loops = FALSE,
                   matrix.type = 'adjacency'
)
network.vertex.names(fullnet) <- lasla.rej.idx

## create a attribute indicating whether the SNP is rejected by BH
in_bhnet <- rep(0,n)
for (i in 1:n) {
  in_bhnet[i] = is.element(lasla.rej.idx[i], bh.rej.idx)
}

set.vertex.attribute(fullnet, "in_bhnet", in_bhnet)
node_colors <- rep('', n)
for (i in 1:n) {
  if(get.node.attr(fullnet,'in_bhnet')[i] == 1){
    node_colors[i] <- 'lightblue'
  }
  else{
    node_colors[i] <- 'red'
  }
}

plot.network(fullnet,
             vertex.col = node_colors,
             vertex.cex = 5,
             displaylabels = T,
             label.pos = 5,
             displayisolates = F)




# Remove some of the isolated subnetworks detected by both methods 
remove_nodes <- c(1146,2366,3507,4347,1881,2048,1541,4995,1986,4316,3969)
remove_nodes <- as.integer(remove_nodes)
remove_idx <- rep(length(remove_nodes))
for (i in 1:length(remove_nodes)) {
  remove_idx[i] = which(lasla.rej.idx==remove_nodes[i])
}
submat <- fullmat[-remove_idx,-remove_idx]
sub.rej.idx <-lasla.rej.idx[-remove_idx]

subnet <- as.network(x=submat,
                      directed = FALSE,
                      loops = FALSE,
                      matrix.type = 'adjacency'
)
network.vertex.names(subnet) <- sub.rej.idx

## create a attribute indicating whether the snp is rejected by bh
in_bhnet <- rep(0,dim(submat)[1])
for (i in 1:dim(submat)[1]) {
  in_bhnet[i] = is.element(sub.rej.idx[i], bh.rej.idx)
}

set.vertex.attribute(subnet, "in_bhnet", in_bhnet)
node_colors <- rep('', dim(submat)[1])
for (i in 1:dim(submat)[1]) {
  if(get.node.attr(subnet,'in_bhnet')[i] == 1){
    node_colors[i] <- 'lightblue'
  }
  else{
    node_colors[i] <- 'red2'
  }
}

par(mar = c(0, 0, 0, 0))

plot.network(subnet,
             vertex.col = node_colors,
             vertex.cex = 1.3,
             # displaylabels = T,
             # label.pos = 5,
            displayisolates = F,
            pad=0,             
            edge.col = "#333333"
            )

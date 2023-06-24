# Set working directory to the LASLA folder
setwd("~/GitHub/LASLA")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')


## Real Data Application: Detecting disease related SNPs.
prepare_data <- function(){
  # read the ld info file
  lddata<-read.delim("./data/myld.ld",sep="")

  # read the summary stats for the selected 5000 SNPs
  testgenom<-read.delim("./data/S5000.txt")

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

if (FALSE) {
  prepare_data()
}


rs_snp<-read.delim("./data/IDs.txt", header = FALSE)

m<-5000
mu0<-rep(0,m)
sd0<-rep(1,m)
# fdr_level<-c(0.001,0.01,0.05,0.1)
fdr_level<-c(0.05)
nrep<-length(fdr_level)

de<-matrix(rep(0, m*nrep), nrep, m)

pv<-2*pnorm(-abs(z_snp), 0, 1)

bon.nr<-rep(0,nrep)
bh.nr<-rep(0,nrep)
lasla.dd.nr<-rep(0,nrep)

for (i in 1:nrep) {
  cat("iteration:", i, "\n")
  q<-fdr_level[i]
  
  bh.res<-bh.func(pv, q)

  bw=density(z_snp,bw="SJ-ste")$bw
  pis_lasla <- lasla_pis(z_snp, d_snp/1.4, pv, tau=bh.func(pv,0.8)$th, h=bw, eps=0)
  
  weight <- lasla_weights(z_snp, d_snp/1.4, pis_lasla, mu0, sd0, q, h=bw, eps=0)
  lasla.dd.res <- lasla_thres(pv, pis_lasla, weight, q)
    
  bon.nr[i]<-length(which(pv<=q/m))
  bh.nr[i]<-bh.res$nr
  lasla.dd.nr[i]<-lasla.dd.res$nr
  de[i,]<-lasla.dd.res$de
}

############## Plot the sub-network ####################
if (False){
  install.packages("network")
  library(network)
  install.packages('ergm')
  library(ergm)
}


# under fdr=0.05, analysis the rejections of bh and lasla
bh.rej.idx <- which(bh.res$de == 1)
bh.rej.snp <- testgenom[bh.rej.idx,]
lasla.rej.idx <- which(lasla.dd.res$de == 1)
lasla.rej.snp <- testgenom[lasla.rej.idx,]

## find snps that are rejected bh lasla but not by bh
extra.rej.idx <- setdiff(lasla.rej.idx, bh.rej.idx)
extra.rej.snp <- testgenom[extra.rej.idx,]


## full net rejected by lasla
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

## create a attribute indicating whether the snp is rejected by bh
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
             vertex.cex = 3,
             displaylabels = T,
             label.pos = 5,
             displayisolates = F)


## net by bh
n <- length(bh.rej.idx)
bhmat <- matrix(rep(0, n*n), n,n)
for (i in 1:n) {
  for (j in i:n) {
    bhmat[i,j] = 1-d_snp[bh.rej.idx[i], bh.rej.idx[j]]
    bhmat[j,i] = bhmat[i,j]
  }
}
bhmat <- ceiling(bhmat)
diag(bhmat) <- 0    ## to make sure there is no self-edge
bhnet <- as.network(x=bhmat,
                  directed = FALSE,
                  loops = FALSE,
                  matrix.type = 'adjacency'
)
network.vertex.names(bhnet) <- rs_snp[bh.rej.idx,]



## find some subnets
## delete unwanted nodes from the fullnet
# remove_nodes <- c(3015,4347,3507,1881,2048,1835,898,3765,2778,714,1105,837,4818,621,2520,368,476,920,1730,4112,1559,
#                   4873,2935,581,788,4940,1937,3522,2677,615,1622,2082,3932,1714,1615,291,3913,2711,1863,2344,
#                   207,4926,4927,2828,1076,3347,1541,4995,2760,4030,2543,2356,4594,2070,4367,741,4316,3969,1986,1146,2366)
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

plot.network(subnet,
             vertex.col = node_colors,
             vertex.cex = 1.5,
             # displaylabels = T,
             # label.pos = 5,
            displayisolates = F,
             )

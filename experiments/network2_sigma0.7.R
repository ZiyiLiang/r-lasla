# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(adaptMT)
library(ggdendro)

#######################################
#     Experiment Parameters           #
#######################################
nrep=100
q<-0.05
dis.vec <- seq(from=0, to=1, by=0.2)
np <- length(dis.vec)

m=1200    # number of observations
pis <- rep(0.1, m)


n_cluster.vec <- c(2, 5, 10, 20)
n_cluster<-length(n_cluster.vec)

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

adaptMT.fdr_list <- replicate(n_cluster, rep(0, np), simplify = FALSE)
adaptMT.etp_list <- replicate(n_cluster, rep(0, np), simplify = FALSE)

bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

lasla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)

adapt.fdp<-matrix(rep(0, np*nrep), nrep, np)
adapt.ntp<-matrix(rep(0, np*nrep), nrep, np)

adaptMT.fdp_list<-replicate(n_cluster, matrix(rep(0, np*nrep), nrep, np), simplify = FALSE)
adaptMT.ntp_list<-replicate(n_cluster, matrix(rep(0, np*nrep), nrep, np), simplify = FALSE)

pb <- progress_bar$new(total = nrep)   # show progress bar
for (i in 1:nrep){
  pb$tick()
  # fix the primary sequence, only auxiliary sequence changes
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0 <- rep(0,m)
  mu1 <- rep(3,m)
  sd0 <- rep(1,m)
  sd1 <- rep(1,m)
  x0<-rnorm(m, mean=mu0, sd=sd0)
  x1<-rnorm(m, mean=mu1, sd=sd1)
  x<-(1-theta)*x0+theta*x1
  pv<-2*pnorm(-abs(x), 0, 1)
  
  for (j in 1:np){
    set.seed(i*nrep+j) 
    dis<-dis.vec[j]
    # generate the distance matrix
    d_lasla<-matrix(rep(0,m*m),m,m)
    for (k in 1:m) {
      for (h in min((k+1),m):m) {
        d_lasla[k,h]=(theta[k]==theta[h])*abs(rnorm(1,dis,0.7))+(theta[k]!=theta[h])*abs(rnorm(1,1,0.7))
      }
    }
    for(k in 2:m) {
      for (h in 1:(k-1)) {
        d_lasla[k,h]=d_lasla[h,k]
      }
    }
    d_lasla <- d_lasla*2
    
    pis_lasla<-lasla_pis(x, d_lasla, pval=pv, bh.func(pv,0.8)$th)
    
    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)
    
    weight<-lasla_weights(x,d_lasla,pis_lasla,mu0,sd0, progress = FALSE)
    lasla.dd.res<-lasla_thres(pvs=pv, pis_lasla, ws=weight, q)
    lasla.dd.de<-lasla.dd.res$de
    lasla.dd.fdp[i, j]<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
    lasla.dd.ntp[i, j]<-sum(theta*lasla.dd.de)/sum(theta)
    
    d <- as.dist(d_lasla)
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
      adaptMT.res <-  adapt_xgboost(s,pv,
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
      adaptMT.rej <- which(adaptMT.res$qvals <= q)
      adaptMT.fdp_list[[k]][i,j] <- sum((1-theta)[adaptMT.rej])/max(length(adaptMT.rej),1)
      adaptMT.ntp_list[[k]][i,j] <- sum(theta[adaptMT.rej])/sum(theta)
    }
  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])
  
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[,i])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[,i])
  
  for(k in 1:length(n_cluster.vec)) {
    adaptMT.fdr_list[[k]][i]<-mean(adaptMT.fdp_list[[k]][,i])
    adaptMT.etp_list[[k]][i]<-mean(adaptMT.ntp_list[[k]][,i])
  }
}

nw_fdr2.mthd<-cbind(bh.fdr, lasla.dd.fdr)
nw_etp2.mthd<-cbind(bh.etp, lasla.dd.etp)

nw_fdr2.mthd<-cbind(bh.fdr, lasla.dd.fdr)
nw_etp2.mthd<-cbind(bh.etp, lasla.dd.etp)

adaptMT.names <- character(0)
for(k in 1:length(n_cluster.vec)) {
  method_name <- sprintf("adapt_cluster_%d", n_cluster.vec[k])
  adaptMT.names <- c(adaptMT.names, method_name)
  
  nw_fdr2.mthd <- cbind(nw_fdr2.mthd, adaptMT.fdr_list[[k]])
  colnames(nw_fdr2.mthd)[ncol(nw_fdr2.mthd)] <- paste(method_name, ".fdr", sep = "")
  
  nw_etp2.mthd <- cbind(nw_etp2.mthd, adaptMT.etp_list[[k]])
  colnames(nw_etp2.mthd)[ncol(nw_etp2.mthd)] <- paste(method_name, ".etp", sep = "")
}


#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(dis.vec, nw_fdr2.mthd, type="o", pch=1:4, lwd=2, main="Network-setting2 FDR Comparison", xlab=expression(mu[2]), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LASLA.DD", adaptMT.names), pch=1:3, col=1:6, lwd=2)

matplot(dis.vec, nw_etp2.mthd, type="o", pch=1:4, lwd=2, main="Network-setting2 Power Comparison", xlab=expression(mu[2]), ylab="Power")



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"
results <- data.frame()

method_names <- c("BH","LASLA.DD", adaptMT.names)

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(dis.vec)),
                    FDR = nw_fdr2.mthd[,i],
                    Power = nw_etp2.mthd[,i],
                    Distance = dis.vec)
  results <- rbind(results, tmp)
}
save(results, file=sprintf("%s/network2_sigma0.7.RData", data_dir))
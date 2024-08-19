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
mean.vec <- seq(from=2.5, to=3, by=0.1)
np <- length(mean.vec)

m=1200    # number of observations
pis <- rep(0.1, m)

n_cluster.vec <- c(5, 10, 20)
n_cluster<-length(n_cluster.vec)

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

adaptMT.fdr_list <- replicate(n_cluster, rep(0, np), simplify = FALSE)
adaptMT.etp_list <- replicate(n_cluster, rep(0, np), simplify = FALSE)

bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)

lasla.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
lasla.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)

adaptMT.fdp_list<-replicate(n_cluster, matrix(rep(0, np*nrep), np, nrep), simplify = FALSE)
adaptMT.ntp_list<-replicate(n_cluster, matrix(rep(0, np*nrep), np, nrep), simplify = FALSE)


for (i in 1:np){
  mean<-mean.vec[i]
  cat("\n", "Network_setting with signal strength:", mean)
  pb <- progress_bar$new(total = nrep)   # show progress bar
  
  for (j in 1:nrep)
  {
    pb$tick()
    set.seed(i*nrep+j)  
    theta<-rbinom(m, size=1, prob=pis)
    mu0 <- rep(0,m)        
    mu1 <- rep(mean,m)
    sd0 <- rep(1,m)
    sd1 <- rep(1,m)
    x0<-rnorm(m, mean=mu0, sd=sd0)
    x1<-rnorm(m, mean=mu1, sd=sd1)
    x<-(1-theta)*x0+theta*x1
    pv<-2*pnorm(-abs(x), 0, 1)
    
    
    # generate the distance matrix
    d_lasla<-matrix(rep(0,m*m),m,m)
    for (k in 1:m) {
      for (h in min((k+1),m):m) {
        d_lasla[k,h]<-(theta[k]==theta[h])*abs(rnorm(1,0,0.7))+(theta[k]!=theta[h])*abs(rnorm(1,1,0.7))
      }
    }
    for(k in 2:m) {
      for (h in 1:(k-1)) {
        d_lasla[k,h]<-d_lasla[h,k]
      }
    }
    d_lasla <- d_lasla*2
    pis_lasla<-lasla_pis(x, d_lasla, pval=pv, tau=bh.func(pv,0.8)$th)
    
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
  
  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[i,])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[i,])
  
  for(k in 1:length(n_cluster.vec)) {
    adaptMT.fdr_list[[k]][i]<-mean(adaptMT.fdp_list[[k]][i,])
    adaptMT.etp_list[[k]][i]<-mean(adaptMT.ntp_list[[k]][i,])
  }
}

nw_fdr1.mthd<-cbind(bh.fdr, lasla.dd.fdr)
nw_etp1.mthd<-cbind(bh.etp, lasla.dd.etp)

adaptMT.names <- character(0)
for(k in 1:length(n_cluster.vec)) {
  method_name <- sprintf("adapt_cluster_%d", n_cluster.vec[k])
  adaptMT.names <- c(adaptMT.names, method_name)
  
  nw_fdr1.mthd <- cbind(nw_fdr1.mthd, adaptMT.fdr_list[[k]])
  colnames(nw_fdr1.mthd)[ncol(nw_fdr1.mthd)] <- paste(method_name, ".fdr", sep = "")
  
  nw_etp1.mthd <- cbind(nw_etp1.mthd, adaptMT.etp_list[[k]])
  colnames(nw_etp1.mthd)[ncol(nw_etp1.mthd)] <- paste(method_name, ".etp", sep = "")
}


#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
matplot(mean.vec, nw_fdr1.mthd, type="o", pch=1:4, lwd=2, main="Network-setting1 FDR Comparison", xlab=expression(mu[1]), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LASLA.DD", adaptMT.names), pch=1:3, col=1:6, lwd=2)

matplot(mean.vec, nw_etp1.mthd, type="o", pch=1:4, lwd=2, main="Network-setting1 Power Comparison", xlab=expression(mu[1]), ylab="Power")



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

method_names <- c("BH","LASLA.DD", adaptMT.names)

# Setting 1
results <- data.frame()

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(mean.vec)),
                    FDR = nw_fdr1.mthd[,i],
                    Power = nw_etp1.mthd[,i],
                    Mean = mean.vec)
  results <- rbind(results, tmp)
}
save(results, file=sprintf("%s/network1_sigma0.7.RData", data_dir))

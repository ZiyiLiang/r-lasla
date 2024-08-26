# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(MASS) 


#######################################
#     Experiment Parameters           #
#######################################


nrep=100
q<-0.05
rho.vec <- seq(from=0, to=0.8, by=0.2)       # Correlation within each block
np <- length(rho.vec)

m=1000    # number of observations
pis <- rep(0.1, m)

# Specify the correlation type and parameters
cor_type <- "block"
block_size <- 100

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)

lasla.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
lasla.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)


for (i in 1:np){
  rho<-rho.vec[i]
  cat("\n", "Dependent network setting with base correlation:", rho)
  pb <- progress_bar$new(total = nrep)   # show progress bar
  
  for (j in 1:nrep)
  {
    pb$tick()
    set.seed(i*nrep+j)  
    theta<-rbinom(m, size=1, prob=pis)
    mu0 <- rep(0,m)        
    mu1 <- rep(3,m)
    sd0 <- rep(1,m)
    sd1 <- rep(1,m)
    
    cor_mat <- generate_correlation_matrix(m, cor_type = "block", block_size = block_size, rho = rho)
    x0<-mvrnorm(n = 1, mu = mu0, Sigma = cor_mat)
    x1<-mvrnorm(n = 1, mu = mu1, Sigma = cor_mat)
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
    pis_lasla<-lasla_pis(x, d_lasla, pval=pv, eps = 0, tau=bh.func(pv,0.8)$th)
    
    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)
    
    weight<-lasla_weights(x,d_lasla,pis_lasla,mu0,sd0, eps = 0, progress = FALSE)
    lasla.dd.res<-lasla_thres(pvs=pv, pis_lasla, ws=weight, q)
    lasla.dd.de<-lasla.dd.res$de
    lasla.dd.fdp[i, j]<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
    lasla.dd.ntp[i, j]<-sum(theta*lasla.dd.de)/sum(theta)
    
  }
  
  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[i,])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[i,])
}

dependent_fdr1.mthd<-cbind(bh.fdr, lasla.dd.fdr)
dependent_etp1.mthd<-cbind(bh.etp, lasla.dd.etp)



#######################################
#          Preview Results            #
#######################################

par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
matplot(rho.vec,dependent_fdr1.mthd, type="o", pch=1:4, lwd=2, main="Dependent-setting1 FDR Comparison", xlab=expression(rho), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LASLA.DD"), pch=1:3, col=1:6, lwd=2)

matplot(rho.vec, dependent_etp1.mthd, type="o", pch=1:4, lwd=2, main="Dependent-setting1 Power Comparison", xlab=expression(rho), ylab="Power")



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

save <- TRUE

if (save){
  method_names <- c("BH","LASLA.DD")
  
  # Setting 1
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i],length(rho.vec)),
                      FDR = dependent_fdr1.mthd[,i],
                      Power = dependent_etp1.mthd[,i],
                      Rho = rho.vec)
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/dependent1.RData", data_dir))
}
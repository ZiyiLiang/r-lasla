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
factor <- 1

m <- 1000    # number of observations
bw <- m**(-1/6) # larger bandwidth
pis <- rep(0.1, m)

# Specify the correlation type and parameters
cor_type <- "random"

bh.fdp<-rep(0, nrep)
bh.ntp<-rep(0, nrep)

lasla.dd.fdp<-rep(0, nrep)
lasla.dd.ntp<-rep(0, nrep)

pb <- progress_bar$new(total = nrep)   # show progress bar

for (i in 1:nrep)
{
  pb$tick()
  set.seed(i)  
  theta<-rbinom(m, size=1, prob=pis)
  mu0 <- rep(0,m)        
  mu1 <- rep(3,m)
  sd0 <- rep(1,m)
  sd1 <- rep(1,m)
  
  cor_mat <- generate_correlation_matrix(m, cor_type = cor_type, factor = factor)

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
  pis_lasla<-lasla_pis(x, d_lasla, pval=pv, eps = 0, h = bw, tau=bh.func(pv,0.8)$th)

  bh.res<-bh.func(pv, q)
  bh.de<-bh.res$de
  bh.fdp[i]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
  bh.ntp[i]<-sum(theta*bh.de)/sum(theta)

  weight<-lasla_weights(x,d_lasla,pis_lasla,mu0,sd0, eps = 0, h = bw, progress = FALSE)
  lasla.dd.res<-lasla_thres(pvs=pv, pis_lasla, ws=weight, q)
  lasla.dd.de<-lasla.dd.res$de
  lasla.dd.fdp[i]<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
  lasla.dd.ntp[i]<-sum(theta*lasla.dd.de)/sum(theta)

}

dependent_fdr3.mthd<-cbind(bh.fdp, lasla.dd.fdp)
dependent_etp3.mthd<-cbind(bh.ntp, lasla.dd.ntp)


#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

save <-FALSE

if (save){
  method_names <- c("BH","LASLA.DD")
  
  # Setting 2
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = method_names[i],
                      FDP = dependent_fdr3.mthd[,i],
                      Power = dependent_etp3.mthd[,i])
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/dependent3.RData", data_dir))
}
# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')



#######################################
#     Experiment Parameters           #
#######################################
m<-1000
mean<-3
# Note that the sigma refers to the variance of the Gaussian distribution,
# not the standard deviation.
sigma.vec<-seq(from=0.2, to=1, by=0.2)
pis<-rep(0.1, m)
q<-0.05
np<-length(sigma.vec)
nrep<-1000
pb <- progress_bar$new(total = nrep)   # show progress bar


#------ overall the FDR and power-------------
bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

law.or.fdr<-rep(0, np)
law.or.etp<-rep(0, np)

lasla.or.fdr<-rep(0, np)
lasla.or.etp<-rep(0, np)


#------FDP and power for all exp-------------
bh.fdp<-matrix(rep(0, nrep*np), nrep, np)
bh.ntp<-matrix(rep(0, nrep*np), nrep, np)

law.or.fdp<-matrix(rep(0, nrep*np), nrep, np)
law.or.ntp<-matrix(rep(0, nrep*np), nrep, np)

lasla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)


for (i in 1:nrep)
{
  pb$tick()
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  x0<-rnorm(m, mean=0, sd=1)
  x1_base <- rnorm(m, mean=mean, sd=0.1)
  isPositive<-rbinom(m, size=1, prob=0.5)
  x1_base <- isPositive*x1_base + (1-isPositive)*(-x1_base)
  
  for (j in 1:np)
  {
    set.seed(2000*i+j)
    sigma<-sigma.vec[j]
    isNoisy<-rbinom(m, size=1, prob=1/2)
    mu1<-rep(mean,m)
    sd1<-isNoisy*rep(1,m) + (1-isNoisy)*rep(sigma,m)
    x1<-isNoisy*(x1_base+rnorm(m,0,0.9)) + (1-isNoisy)*(x1_base+rnorm(m,0,sigma-0.1))
    x<-(1-theta)*x0+theta*x1
    pv<-2*pnorm(-abs(x), 0, 1)
    
    
    qq=q/(1-pii)
    bh.res<-bh.func(pv, qq)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)
    
    wor_laws <- pis/(1-pis)
    law.or.res<-lasla_thres(pvs=pv, pis=pis, ws=wor_laws, q)
    law.or.de<-law.or.res$de
    law.or.fdp[i, j]<-sum((1-theta)*law.or.de)/max(sum(law.or.de), 1)
    law.or.ntp[i, j]<-sum(theta*law.or.de)/sum(theta)
    
    wor_lasla <- lasla_oracle_weights.func(x, pis, dist_type="normal", q, mu1=mu1, sd1=sd1)
    lasla.or.res<-lasla_thres(pvs=pv, pis=pis, ws=wor_lasla, q)
    lasla.or.de<-lasla.or.res$de
    lasla.or.fdp[i, j]<-sum((1-theta)*lasla.or.de)/max(sum(lasla.or.de), 1)
    lasla.or.ntp[i, j]<-sum(theta*lasla.or.de)/sum(theta)
  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])
  
  lasla.or.fdr[i]<-mean(lasla.or.fdp[,i])
  lasla.or.etp[i]<-mean(lasla.or.ntp[,i])
  
  law.or.fdr[i]<-mean(law.or.fdp[,i])
  law.or.etp[i]<-mean(law.or.ntp[,i])
}

altshape_fdr.mthd<-cbind(bh.fdr, lasla.or.fdr, law.or.fdr)
altshape_etp.mthd<-cbind(bh.etp, lasla.or.etp, law.or.etp)



#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(sigma.vec, altshape_fdr.mthd, type="o", pch=1:4, lwd=2, main="FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("PV.OR",  "LASLA.OR", "LAWS.OR"), pch=1:3, col=1:3, lwd=2)

matplot(sigma.vec, altshape_etp.mthd, type="o", pch=1:3, lwd=2, main="Power Comparison", xlab=expression(sigma), ylab="Power")



#######################################
#             Save Results            #
#######################################
save = FALSE

if (save){
  data_dir <- "./results"
  
  method_names <- c("PV.OR", "LASLA.OR", "LAWS.OR")
  
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i],length(sigma.vec)),
                      FDR = altshape_fdr.mthd[,i],
                      Power = altshape_etp.mthd[,i],
                      Sigma = sigma.vec)
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/alt_dist.RData", data_dir))
}
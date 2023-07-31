# Set working directory to the LASLA folder
setwd("~/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')



#######################################
#     Experiment Parameters           #
#######################################
m<-1000
mean<-3
gamma.vec<-seq(from=0.5, to=1, by=0.1)  # controls the asymmetry level
pis<-rep(0.1, m)
q<-0.05
np<-length(gamma.vec)
nrep<-200
pb <- progress_bar$new(total = nrep)   # show progress bar

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)
bh.fdp<-matrix(rep(0, nrep*np), nrep, np)
bh.ntp<-matrix(rep(0, nrep*np), nrep, np)

law.or.fdr<-rep(0, np)
law.or.etp<-rep(0, np)
law.or.fdp<-matrix(rep(0, nrep*np), nrep, np)
law.or.ntp<-matrix(rep(0, nrep*np), nrep, np)

lasla.or.fdr<-rep(0, np)
lasla.or.etp<-rep(0, np)
lasla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)

for (i in 1:nrep)
{
  pb$tick()
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  x0<-rnorm(m, mean=0, sd=1)
  x1_base <- rnorm(m, mean=mean, sd=1)
 
   for (j in 1:np)
  { 
    # Make sure experiments are reproducible
    set.seed(1000*i+j)   
    gamma<-gamma.vec[j]
    gamma<-rep(gamma,m)
    isPositive<-rbinom(m, size=1, prob=gamma)
    mu0<-rep(0,m)
    mu1<-rep(mean,m)
    sd0<-rep(1,m)
    sd1<-rep(1,m)
    x1 <- isPositive*x1_base + (1-isPositive)*(-x1_base)
    x<-(1-theta)*x0+theta*x1
    pv<-2*pnorm(-abs(x), 0, 1)
    
    # generating the distance matrix
    d<-matrix(rep(0,m^2),m,m)
    for (k in 1:m) {
      d[k,]<-abs(1:m-k)
    }
    
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
    
    awor <- awor_asymmetric.func(x, pis, gamma, mu0, mu1, sd0, sd1, q)
    lasla.or.res<-lasla_thres(pvs=pv, pis=pis, ws=awor, q)
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

asymmetric_fdr.mthd<-cbind(bh.fdr, lasla.or.fdr, law.or.fdr)
asymmetric_etp.mthd<-cbind(bh.etp, lasla.or.etp, law.or.etp)



#######################################
#          Preview Results            #
#######################################
par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(gamma.vec, asymmetric_fdr.mthd, type="o", pch=1:4, lwd=2, main="FDR Comparison", xlab=expression(gamma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("PV.OR",  "LASLA.OR", "LAWS.OR"), pch=1:3, col=1:3, lwd=2)

matplot(gamma.vec, asymmetric_etp.mthd, type="o", pch=1:3, lwd=2, main="Power Comparison", xlab=expression(gamma), ylab="Power")



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

method_names <- c("PV.OR", "LASLA.OR", "LAWS.OR")

results <- data.frame()

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(gamma.vec)),
                    FDR = asymmetric_fdr.mthd[,i],
                    Power = asymmetric_etp.mthd[,i],
                    Gamma = gamma.vec)
  results <- rbind(results, tmp)
}
save(results, file=sprintf("%s/asymmetric.RData", data_dir))



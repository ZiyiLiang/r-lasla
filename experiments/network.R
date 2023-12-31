# Set working directory to the LASLA folder
setwd("~/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')


#######################################
#     Experiment Parameters           #
#######################################
nrep=100
q<-0.05
mean.vec <- seq(from=2.5, to=3, by=0.1)
np <- length(mean.vec)

m=1200    # number of observations
pis <- rep(0.1, m)

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)
lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)


bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)
lasla.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
lasla.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)


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
        d_lasla[k,h]<-(theta[k]==theta[h])*abs(rnorm(1,0,0.5))+(theta[k]!=theta[h])*abs(rnorm(1,1,0.5))
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

  }

  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[i,])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[i,])
}

nw_fdr1.mthd<-cbind(bh.fdr, lasla.dd.fdr)
nw_etp1.mthd<-cbind(bh.etp, lasla.dd.etp)


#######################################
#     Experiment Parameters           #
#######################################
nrep=100
q<-0.05
dis.vec <- seq(from=0, to=1, by=0.2)
np <- length(dis.vec)

m=1200    # number of observations
pis <- rep(0.1, m)

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)
lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)

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
    dis<-dis.vec[j]
    # generate the distance matrix
    d_lasla<-matrix(rep(0,m*m),m,m)
    for (k in 1:m) {
      for (h in min((k+1),m):m) {
        d_lasla[k,h]=(theta[k]==theta[h])*abs(rnorm(1,dis,0.5))+(theta[k]!=theta[h])*abs(rnorm(1,1,0.5))
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

  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])

  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[,i])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[,i])
}

nw_fdr2.mthd<-cbind(bh.fdr, lasla.dd.fdr)
nw_etp2.mthd<-cbind(bh.etp, lasla.dd.etp)



#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
matplot(mean.vec, nw_fdr1.mthd, type="o", pch=1:4, lwd=2, main="Network-setting1 FDR Comparison", xlab=expression(mu[1]), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LASLA.DD"), pch=1:3, col=1:6, lwd=2)

matplot(mean.vec, nw_etp1.mthd, type="o", pch=1:4, lwd=2, main="Network-setting1 Power Comparison", xlab=expression(mu[1]), ylab="Power")

matplot(dis.vec, nw_fdr2.mthd, type="o", pch=1:4, lwd=2, main="Network-setting2 FDR Comparison", xlab=expression(mu[2]), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LASLA.DD"), pch=1:3, col=1:6, lwd=2)

matplot(dis.vec, nw_etp2.mthd, type="o", pch=1:4, lwd=2, main="Network-setting2 Power Comparison", xlab=expression(mu[2]), ylab="Power")



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

method_names <- c("BH", "LASLA.DD")

# Setting 1
results <- data.frame()

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(mean.vec)),
                    FDR = nw_fdr1.mthd[,i],
                    Power = nw_etp1.mthd[,i],
                    Mean = mean.vec)
  results <- rbind(results, tmp)
}
save(results, file=sprintf("%s/network1.RData", data_dir))


# Setting 2
results <- data.frame()

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(dis.vec)),
                    FDR = nw_fdr2.mthd[,i],
                    Power = nw_etp2.mthd[,i],
                    Distance = dis.vec)
  results <- rbind(results, tmp)
}
save(results, file=sprintf("%s/network2.RData", data_dir))
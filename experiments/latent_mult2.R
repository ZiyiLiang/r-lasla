# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(adaptMT)

#######################################
#           Setting 2                 #
#######################################
set.seed(0)

m<-1200
noise_l.vec<-seq(from=2.5, to=3, by=0.5)
pis<-rep(0.1, m)
q<-0.05
np<-length(noise_l.vec)
nrep<-50

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)
or.fdr<-rep(0, np)
or.etp<-rep(0, np)

lasla.or.fdr<-rep(0, np)
lasla.or.etp<-rep(0, np)
lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

avg.fdr<-rep(0, np)
avg.etp<-rep(0, np)
adapt.fdr<-rep(0, np)
adapt.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

or.fdp<-matrix(rep(0, np*nrep), nrep, np)
or.ntp<-matrix(rep(0, np*nrep), nrep, np)

lasla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)
avg.fdp<-matrix(rep(0, np*nrep), nrep, np)
avg.ntp<-matrix(rep(0, np*nrep), nrep, np)

adapt.fdp<-matrix(rep(0, np*nrep), nrep, np)
adapt.ntp<-matrix(rep(0, np*nrep), nrep, np)

pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = nrep,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE, 
                       width = 100)      
for (i in 1:nrep){
  pb$tick()
  
  theta<-rbinom(m, size=1, prob=pis)
  dis_theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(3,m)
  sd0<-rep(1,m) 
  sd1<-rep(sqrt(1+1^2),m)
  l0<-rep(0,m)
  l1<-rnorm(m, mean=3, sd=1)
  dis_l0<-rnorm(m, mean=0, sd=1)
  dis_l1<-rnorm(m, mean=3, sd=1)
  latent<-(1-theta)*l0+theta*l1
  dis_latent<-(1-dis_theta)*dis_l0+dis_theta*dis_l1
  x<-rnorm(m,mean=latent, sd=1)        # Primary data
  for (j in 1:np)
  {
    noise<-noise_l.vec[j]
    s1<-rnorm(m, mean=latent, sd=noise)  # Two informative auxiliary sequences
    s2<-rnorm(m, mean=latent, sd=noise)
    s3<-rnorm(m, mean=dis_latent, sd=noise)  # Two irrelevant auxiliary sequences
    s4<-rnorm(m, mean=dis_latent, sd=noise)
    pv<-2*pnorm(-abs(x), 0, 1)
    
    # distance matrix for lasla
    S<-cbind(s1,s2,s3,s4)
    R<-cov(S)
    d_lasla<-matrix(rep(0,m^2),m,m)
    for (k in 1:m) {
      d_lasla[k,]=mahalanobis(S,S[k,],R)
    }
    d_lasla[lower.tri(d_lasla)] <- t(d_lasla)[lower.tri(d_lasla)]
    d_lasla <- normalize_distance(x, d_lasla)

    # distance matrix for averaging method
    d_avg<-matrix(rep(0,m^2),m,m)
    s_avg<-(s1+s2+s3+s4)/4
    for (k in 1:m) {
      d_avg[k,]<-abs(s_avg-s_avg[k])
    }
    d_avg <- normalize_distance(x, d_avg)

    pis_lasla<-lasla_pis(x, d_lasla, pval=pv, tau=bh.func(pv,0.8)$th)
    pis_avg<-lasla_pis(x, d_avg, pval=pv, tau=bh.func(pv,0.8)$th)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    lpis<-latent_pis_mult2.func(pis, s1,s2,s3,s4, mu1,noise)
    awor <- latent_awor_mult2.func(x,s1,s2,s3,s4,pis,mu1,noise,mu0,sd0)
    lasla.or.res<-lasla_thres(pvs=pv, lpis, ws=awor, q)
    lasla.or.de<-lasla.or.res$de
    lasla.or.fdp[i, j]<-sum((1-theta)*lasla.or.de)/max(sum(lasla.or.de), 1)
    lasla.or.ntp[i, j]<-sum(theta*lasla.or.de)/sum(theta)

    weight<-lasla_weights(x,d_lasla,pis_lasla,mu0,sd0, progress=FALSE)
    lasla.dd.res<-lasla_thres(pvs=pv, pis_lasla, ws=weight, q)
    lasla.dd.de<-lasla.dd.res$de
    lasla.dd.fdp[i, j]<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
    lasla.dd.ntp[i, j]<-sum(theta*lasla.dd.de)/sum(theta)

    weight1<-lasla_weights(x,d_avg,pis_avg,mu0,sd0,  progress=FALSE)
    avg.res<-lasla_thres(pvs=pv, pis_avg, ws=weight1, q)
    avg.de<-avg.res$de
    avg.fdp[i, j]<-sum((1-theta)*avg.de)/max(sum(avg.de), 1)
    avg.ntp[i, j]<-sum(theta*avg.de)/sum(theta)

    adapt.res <-  adapt_xgboost(S ,pv,
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
    adapt.rej <- which(adapt.res$qvals <= q)
    adapt.fdp[i,j] <- sum((1-theta)[adapt.rej])/max(length(adapt.rej),1)
    adapt.ntp[i,j] <- sum(theta[adapt.rej])/sum(theta)
  }
}


for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])
  lasla.or.fdr[i]<-mean(lasla.or.fdp[,i])
  lasla.or.etp[i]<-mean(lasla.or.ntp[,i])
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[,i])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[,i])
  avg.fdr[i]<-mean(avg.fdp[,i])
  avg.etp[i]<-mean(avg.ntp[,i])
  adapt.fdr[i]<-mean(adapt.fdp[,i])
  adapt.etp[i]<-mean(adapt.ntp[,i])
}
fdr_mult2.mthd<-cbind(bh.fdr, lasla.or.fdr, lasla.dd.fdr, avg.fdr, adapt.fdr)
etp_mult2.mthd<-cbind(bh.etp, lasla.or.etp, lasla.dd.etp, avg.etp, adapt.etp)


#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(noise_l.vec, fdr_mult2.mthd, type="o", pch=1:7, lwd=2, main="Mult-setting2 FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH", "LASLA.OR","LASLA.DD","AVG","ADAPT"), pch=1:6, col=1:6, lwd=2)

matplot(noise_l.vec, etp_mult2.mthd, type="o", pch=1:7, lwd=2, main="Mult-setting2 Power Comparison", xlab=expression(sigma), ylab="Power")



#######################################
#             Save Results            #
#######################################
save <- FALSE

if (save){
  data_dir <- "./results"
  
  method_names <- c("BH", "LASLA.OR","LASLA.DD","AVG", "ADAPT")
  
  # Setting 2
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i],length(noise_l.vec)),
                      FDR = fdr_mult2.mthd[,i],
                      Power = etp_mult2.mthd[,i],
                      Sigma = noise_l.vec)
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/latent_mult2.RData", data_dir))
}
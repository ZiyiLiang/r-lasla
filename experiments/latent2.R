# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')


#######################################
#     Experiment Parameters           #
#######################################
m<-1200
mu_l.vec<-seq(from=3, to=4, by=0.20)
noise<-1
pis<-rep(0.1, m)
q<-0.05
np<-length(mu_l.vec)
nrep<-100


bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

lasla.or.fdr<-rep(0, np)
lasla.or.etp<-rep(0, np)
lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

wbh.fdr<-rep(0, np)
wbh.etp<-rep(0, np)

awbh.fdr<-rep(0, np)
awbh.etp<-rep(0, np)

sabha.fdr<-rep(0, np)
sabha.etp<-rep(0, np)

adaptMT.fdr<-rep(0, np)
adaptMT.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, np*nrep), np, nrep)
bh.ntp<-matrix(rep(0, np*nrep), np, nrep)

lasla.or.fdp<-matrix(rep(0, np*nrep), np, nrep)
lasla.or.ntp<-matrix(rep(0, np*nrep), np, nrep)
lasla.dd.fdp<-matrix(rep(0, np*nrep), np, nrep)
lasla.dd.ntp<-matrix(rep(0, np*nrep), np, nrep)

wbh.fdp<-matrix(rep(0, np*nrep), np, nrep)
wbh.ntp<-matrix(rep(0, np*nrep), np, nrep)
awbh.fdp<-matrix(rep(0, np*nrep), np, nrep)
awbh.ntp<-matrix(rep(0, np*nrep), np, nrep)

sabha.fdp<-matrix(rep(0, np*nrep), np, nrep)
sabha.ntp<-matrix(rep(0, np*nrep), np, nrep)

adaptMT.fdp<-matrix(rep(0, np*nrep), np, nrep)
adaptMT.ntp<-matrix(rep(0, np*nrep), np, nrep)

for (i in 1:np)
{
  cat("\n", "Latent setting with signal strength:", mu_l.vec[i])
  pb <- progress_bar$new(total = nrep)   # show progress bar
  
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(mu_l.vec[i],m)
  sd0<-rep(1,m)  
  sd1<-rep(sqrt(1+1^2),m)

  for (j in 1:nrep)
  {
    pb$tick()
    l0<-rep(0,m)
    l1<-rnorm(m, mean=mu_l.vec[i], sd=1)
    latent<-(1-theta)*l0+theta*l1
    x<-rnorm(m,mean=latent, sd=1)       # Primary data
    s<-rnorm(m, mean=latent, sd=noise)  # Auxillary data
    pv<-2*pnorm(-abs(x), 0, 1)

    # Compute distance matrix
    d<-matrix(rep(0,m^2),m,m)
    for (k in 1:m) {
      d[k,]<-abs(s-s[k])
    }

    pis_lasla <- lasla_pis(x, d=d, pval=pv, tau=bh.func(pv,0.8)$th)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    lpis<-latent_pis.func(pis, s, mu1,noise)
    awor <- latent_awor.func(x,s,pis,mu1,noise,mu0,sd0)
    lasla.or.res<-lasla_thres(pvs=pv, lpis, ws=awor, q)
    lasla.or.de<-lasla.or.res$de
    lasla.or.fdp[i, j]<-sum((1-theta)*lasla.or.de)/max(sum(lasla.or.de), 1)
    lasla.or.ntp[i, j]<-sum(theta*lasla.or.de)/sum(theta)

    weight<-lasla_weights(x,d,pis_lasla,mu0,sd0, progress = FALSE)
    lasla.dd.res<-lasla_thres(pvs=pv, pis=pis_lasla, ws=weight, q)
    lasla.dd.de<-lasla.dd.res$de
    lasla.dd.fdp[i, j]<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
    lasla.dd.ntp[i, j]<-sum(theta*lasla.dd.de)/sum(theta)

    adaptMT.res <-  adapt_xgboost(matrix(s, m, 1) ,pv,
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
    adaptMT.fdp[i,j] <- sum((1-theta)[adaptMT.rej])/max(length(adaptMT.rej),1)
    adaptMT.ntp[i,j] <- sum(theta[adaptMT.rej])/sum(theta)
    
    sweight<-weight*m/sum(weight)
    wbh.res<-bh.func(pv/sweight, q)
    wbh.de<-wbh.res$de
    wbh.fdp[i, j]<-sum((1-theta)*wbh.de)/max(sum(wbh.de), 1)
    wbh.ntp[i, j]<-sum(theta*wbh.de)/sum(theta)
  
    # qq=q/(1-pii)
    ll <- bh.func(pv,0.8)$th
    phat <- (max(sweight) + sum(sweight * (pv > ll))) / (m * (1 - ll)) 
    awbh.res<-bh.func(pv/sweight, q/phat)
    awbh.de<-awbh.res$de
    awbh.fdp[i, j]<-sum((1-theta)*awbh.de)/max(sum(awbh.de), 1)
    awbh.ntp[i, j]<-sum(theta*awbh.de)/sum(theta)
    
    ws_sabha<-1/(1-pis_lasla)
    sabha.res<-bh.func(pv/ws_sabha, q)
    sabha.de<-sabha.res$de
    sabha.fdp[i, j]<-sum((1-theta)*sabha.de)/max(sum(sabha.de), 1)
    sabha.ntp[i, j]<-sum(theta*sabha.de)/sum(theta)
  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])

  lasla.or.fdr[i]<-mean(lasla.or.fdp[i,])
  lasla.or.etp[i]<-mean(lasla.or.ntp[i,])
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[i,])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[i,])

  wbh.fdr[i]<-mean(wbh.fdp[i,])
  wbh.etp[i]<-mean(wbh.ntp[i,])
  
  awbh.fdr[i]<-mean(awbh.fdp[i,])
  awbh.etp[i]<-mean(awbh.ntp[i,])
  
  adaptMT.fdr[i]<-mean(adaptMT.fdp[i,])
  adaptMT.etp[i]<-mean(adaptMT.ntp[i,])
  
  sabha.fdr[i]<-mean(sabha.fdp[i,])
  sabha.etp[i]<-mean(sabha.ntp[i,])
}
fdr_tf2.mthd<-cbind(bh.fdr, lasla.or.fdr, lasla.dd.fdr, adaptMT.fdr, sabha.fdr, wbh.fdr, awbh.fdr)
etp_tf2.mthd<-cbind(bh.etp, lasla.or.etp, lasla.dd.etp, adaptMT.etp, sabha.etp, wbh.etp, awbh.etp)



#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(mu_l.vec, fdr_tf2.mthd, type="o", pch=1:4, lwd=2, main="TF-setting2 FDR Comparison", xlab=expression(mu), ylab="FDR", ylim=c(0.0, 0.095))
legend("top", c("BH", "LASLA.OR","LASLA.DD", "AdaptMT", "SABHA", "WBH", "adaptive-WBH"), pch=1:4, col=1:8, lwd=2)

matplot(mu_l.vec, etp_tf2.mthd, type="o", pch=1:4, lwd=2, main="TF-setting2 Power Comparison", xlab=expression(mu), ylab="Power")



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

save = TRUE

if (save){
  method_names <- c("BH", "LASLA.OR","LASLA.DD", "ADAPT", "SABHA", "WBH", "AWBH")
  
  # Setting 1
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i],length(mu_l.vec)),
                      FDR = fdr_tf2.mthd[,i],
                      Power = etp_tf2.mthd[,i],
                      Mean = mu_l.vec)
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/latent2.RData", data_dir))
}

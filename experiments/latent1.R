# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

#######################################
#     Experiment Parameters           #
#######################################
m<-1200
#noise_l.vec<-seq(from=0.5, to=2, by=0.25)
noise_l.vec<-seq(from=2, to=2, by=0.25)
pis<-rep(0.1, m)
q<-0.05
np<-length(noise_l.vec)
nrep <- 100

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

lasla.or.fdr<-rep(0, np)
lasla.or.etp<-rep(0, np)
lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

adaptMT.fdr<-rep(0, np)
adaptMT.etp<-rep(0, np)

wbh.fdr<-rep(0, np)
wbh.etp<-rep(0, np)

sabha.fdr<-rep(0, np)
sabha.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

lasla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)

adaptMT.fdp<-matrix(rep(0, np*nrep), nrep, np)
adaptMT.ntp<-matrix(rep(0, np*nrep), nrep, np)

wbh.fdp<-matrix(rep(0, np*nrep), nrep, np)
wbh.ntp<-matrix(rep(0, np*nrep), nrep, np)

sabha.fdp<-matrix(rep(0, np*nrep), nrep, np)
sabha.ntp<-matrix(rep(0, np*nrep), nrep, np)

pb <- progress_bar$new(total = nrep)   # show progress bar
for (i in 1:nrep)
{
  pb$tick()
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(2.5,m)
  sd0<-rep(1,m)  
  sd1<-rep(sqrt(1+1^2),m)
  l0<-rep(0,m)
  l1<-rnorm(m, mean=2.5, sd=1)
  latent<-(1-theta)*l0+theta*l1
  x<-rnorm(m,mean=latent, sd=1)   # primary statistics
  for (j in 1:np)
  {
    set.seed(i*nrep+j)
    noise<-noise_l.vec[j]
    s<-rnorm(m, mean=latent, sd=noise)  # auxiliary statistics
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

    ws_sabha<-1/(1-pis_lasla)
    sabha.res<-bh.func(pv/ws_sabha, q)
    sabha.de<-sabha.res$de
    sabha.fdp[i, j]<-sum((1-theta)*sabha.de)/max(sum(sabha.de), 1)
    sabha.ntp[i, j]<-sum(theta*sabha.de)/sum(theta)

  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])

  lasla.or.fdr[i]<-mean(lasla.or.fdp[,i])
  lasla.or.etp[i]<-mean(lasla.or.ntp[,i])
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[,i])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[,i])
  
  adaptMT.fdr[i]<-mean(adaptMT.fdp[,i])
  adaptMT.etp[i]<-mean(adaptMT.ntp[,i])

  wbh.fdr[i]<-mean(wbh.fdp[,i])
  wbh.etp[i]<-mean(wbh.ntp[,i])

  sabha.fdr[i]<-mean(sabha.fdp[,i])
  sabha.etp[i]<-mean(sabha.ntp[,i])
}
fdr_tf1.mthd<-cbind(bh.fdr, lasla.or.fdr, lasla.dd.fdr, adaptMT.fdr, sabha.fdr, wbh.fdr)
etp_tf1.mthd<-cbind(bh.etp, lasla.or.etp, lasla.dd.etp, adaptMT.etp, sabha.etp, wbh.etp)



#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(noise_l.vec, fdr_tf1.mthd, type="o", pch=1:4, lwd=2, main="TF-setting1 FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0, 0.09))

matplot(noise_l.vec, etp_tf1.mthd, type="o", pch=1:4, lwd=2, main="TF-setting1 Power Comparison", xlab=expression(sigma), ylab="Power")
legend("topright", c("BH", "LASLA.OR","LASLA.DD", "AdaptMT", "SABHA", "WBH"), pch=1:4, col=1:7, lwd=2)



#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

method_names <- c("BH", "LASLA.OR","LASLA.DD", "ADAPTMT","SABHA", "WBH")

results <- data.frame()

for (i in 1:length(method_names)) {
  tmp <- data.frame(Method = rep(method_names[i],length(noise_l.vec)),
                    FDR = fdr_tf1.mthd[,i],
                    Power = etp_tf1.mthd[,i],
                    Sigma = noise_l.vec)
  results <- rbind(results, tmp)
}
save(results, file=sprintf("%s/latent1.RData", data_dir))
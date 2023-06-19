source('G:/My Drive/LATLA/LATLA-code/latla_funcs.R')
## TF setting: latent variable


## Simulation 1: noise level of the auxillary data change from 0.5 to 2 by 0.25

m<-1200
noise_l.vec<-seq(from=0.5, to=2, by=0.25)
pis<-rep(0.1, m)
q<-0.05
np<-length(noise_l.vec)
nrep<-3

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

bc.fdr<-rep(0, np)
bc.etp<-rep(0, np)

latent_or.fdr<-rep(0, np)
latent_or.etp<-rep(0, np)

latla.or.fdr<-rep(0, np)
latla.or.etp<-rep(0, np)
latla.dd.fdr<-rep(0, np)
latla.dd.etp<-rep(0, np)


wbh.fdr<-rep(0, np)
wbh.etp<-rep(0, np)
latla.dd2.fdr<-rep(0, np)
latla.dd2.etp<-rep(0, np)


sabha.fdr<-rep(0, np)
sabha.etp<-rep(0, np)


bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

bc.fdp<-matrix(rep(0, np*nrep), nrep, np)
bc.ntp<-matrix(rep(0, np*nrep), nrep, np)


latent_or.fdp<-matrix(rep(0, np*nrep), nrep, np)
latent_or.ntp<-matrix(rep(0, np*nrep), nrep, np)

latla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)


wbh.fdp<-matrix(rep(0, np*nrep), nrep, np)
wbh.ntp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd2.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd2.ntp<-matrix(rep(0, np*nrep), nrep, np)

sabha.fdp<-matrix(rep(0, np*nrep), nrep, np)
sabha.ntp<-matrix(rep(0, np*nrep), nrep, np)

for (i in 1:nrep)
{
  cat("\n", "TF iteration i= ", i, "\n", "iteration j=")

  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(2.5,m)
  sd0<-rep(1,m)  #the auxillary source sd=0.1
  sd1<-rep(sqrt(1+1^2),m)
  #generate primary and auxillary data
  l0<-rep(0,m)
  l1<-rnorm(m, mean=2.5, sd=1)
  latent<-(1-theta)*l0+theta*l1
  x<-rnorm(m,mean=latent, sd=1)
  for (j in 1:np)
  {
    cat(j)
    noise<-noise_l.vec[j]
    s<-rnorm(m, mean=latent, sd=noise)  # variation from source when collecting auxillary data
    pv<-2*pnorm(-abs(x), 0, 1)

    # generating the distance matrix
    d<-matrix(rep(0,m^2),m,m)
    for (k in 1:m) {
      d[k,]<-abs(s-s[k])
    }

    pis_latla <- latla_pis(x, d=d, pval=pv, tau=bh.func(pv,0.8)$th)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    lpis<-latent_pis.func(pis, s, mu1,noise)
    awor <- latent_awor.func(x,s,pis,mu1,noise,mu0,sd0)
    latla.or.res<-latla_thres(pvs=pv, lpis, ws=awor, q)
    latla.or.de<-latla.or.res$de
    latla.or.fdp[i, j]<-sum((1-theta)*latla.or.de)/max(sum(latla.or.de), 1)
    latla.or.ntp[i, j]<-sum(theta*latla.or.de)/sum(theta)

    # weight<-latla_weights(x,d,pis_latla,mu0,sd0)
    # latla.dd.res<-latla_thres(pvs=pv, pis=pis_latla, ws=weight, q)
    # latla.dd.de<-latla.dd.res$de
    # latla.dd.fdp[i, j]<-sum((1-theta)*latla.dd.de)/max(sum(latla.dd.de), 1)
    # latla.dd.ntp[i, j]<-sum(theta*latla.dd.de)/sum(theta)
# 
#     sweight<-weight*m/sum(weight)
#     wbh.res<-bh.func(pv/sweight, q)
#     wbh.de<-wbh.res$de
#     wbh.fdp[i, j]<-sum((1-theta)*wbh.de)/max(sum(wbh.de), 1)
#     wbh.ntp[i, j]<-sum(theta*wbh.de)/sum(theta)
# 
#     ws_sabha<-1/(1-pis_latla)
#     sabha.res<-bh.func(pv/ws_sabha, q)
#     sabha.de<-sabha.res$de
#     sabha.fdp[i, j]<-sum((1-theta)*sabha.de)/max(sum(sabha.de), 1)
#     sabha.ntp[i, j]<-sum(theta*sabha.de)/sum(theta)

  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])

  latla.or.fdr[i]<-mean(latla.or.fdp[,i])
  latla.or.etp[i]<-mean(latla.or.ntp[,i])
  latla.dd.fdr[i]<-mean(latla.dd.fdp[,i])
  latla.dd.etp[i]<-mean(latla.dd.ntp[,i])

  wbh.fdr[i]<-mean(wbh.fdp[,i])
  wbh.etp[i]<-mean(wbh.ntp[,i])

  sabha.fdr[i]<-mean(sabha.fdp[,i])
  sabha.etp[i]<-mean(sabha.ntp[,i])
}
fdr_tf1.mthd<-cbind(bh.fdr, latla.or.fdr, latla.dd.fdr, sabha.fdr, wbh.fdr)
etp_tf1.mthd<-cbind(bh.etp, latla.or.etp, latla.dd.etp, sabha.etp, wbh.etp)



par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(noise_l.vec, fdr_tf1.mthd, type="o", pch=1:4, lwd=2, main="TF-setting1 FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0, 0.09))

matplot(noise_l.vec, etp_tf1.mthd, type="o", pch=1:4, lwd=2, main="TF-setting1 Power Comparison", xlab=expression(sigma), ylab="Power")
legend("topright", c("BH", "LATLA.OR","LATLA.DD", "SABHA", "WBH"), pch=1:4, col=1:6, lwd=2)

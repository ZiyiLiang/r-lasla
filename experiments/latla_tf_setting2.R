source('C:/Users/liang/OneDrive/Desktop/LATLA/LATLA-code/latla_funcs.R')
## TF setting: latent variable


## Simulation 2: signal strength of the latent variable

m<-1200
mu_l.vec<-seq(from=3, to=4, by=0.20)
noise<-1
pis<-rep(0.1, m)
q<-0.05
np<-length(mu_l.vec)
nrep<-100


bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

latla.or.fdr<-rep(0, np)
latla.or.etp<-rep(0, np)
latla.dd.fdr<-rep(0, np)
latla.dd.etp<-rep(0, np)

wbh.fdr<-rep(0, np)
wbh.etp<-rep(0, np)
sabha.fdr<-rep(0, np)
sabha.etp<-rep(0, np)


bh.fdp<-matrix(rep(0, np*nrep), np, nrep)
bh.ntp<-matrix(rep(0, np*nrep), np, nrep)

latla.or.fdp<-matrix(rep(0, np*nrep), np, nrep)
latla.or.ntp<-matrix(rep(0, np*nrep), np, nrep)
latla.dd.fdp<-matrix(rep(0, np*nrep), np, nrep)
latla.dd.ntp<-matrix(rep(0, np*nrep), np, nrep)

wbh.fdp<-matrix(rep(0, np*nrep), np, nrep)
wbh.ntp<-matrix(rep(0, np*nrep), np, nrep)
sabha.fdp<-matrix(rep(0, np*nrep), np, nrep)
sabha.ntp<-matrix(rep(0, np*nrep), np, nrep)


for (i in 1:np)
{
  cat("\n", "TF iteration i= ", i, "\n", "iteration j=")

  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(mu_l.vec[i],m)
  sd0<-rep(1,m)  #the auxillary source sd=0.1
  sd1<-rep(sqrt(1+1^2),m)
  #generate primary and auxillary data

  for (j in 1:nrep)
  {
    cat(j)
    l0<-rep(0,m)
    l1<-rnorm(m, mean=mu_l.vec[i], sd=1)
    latent<-(1-theta)*l0+theta*l1
    x<-rnorm(m,mean=latent, sd=1)
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

    # using the asymmetricity adapted weights
    lpis<-latent_pis.func(pis, s, mu1,noise)
    awor <- latent_awor.func(x,s,pis,mu1,noise,mu0,sd0)
    latla.or.res<-latla_thres(pvs=pv, lpis, ws=awor, q)
    latla.or.de<-latla.or.res$de
    latla.or.fdp[i, j]<-sum((1-theta)*latla.or.de)/max(sum(latla.or.de), 1)
    latla.or.ntp[i, j]<-sum(theta*latla.or.de)/sum(theta)

    weight<-latla_weights(x,d,pis_latla,mu0,sd0)
    latla.dd.res<-latla_thres(pvs=pv, pis=pis_latla, ws=weight, q)
    latla.dd.de<-latla.dd.res$de
    latla.dd.fdp[i, j]<-sum((1-theta)*latla.dd.de)/max(sum(latla.dd.de), 1)
    latla.dd.ntp[i, j]<-sum(theta*latla.dd.de)/sum(theta)

    sweight<-weight*m/sum(weight)
    wbh.res<-bh.func(pv/sweight, q)
    wbh.de<-wbh.res$de
    wbh.fdp[i, j]<-sum((1-theta)*wbh.de)/max(sum(wbh.de), 1)
    wbh.ntp[i, j]<-sum(theta*wbh.de)/sum(theta)

    ws_sabha<-1/(1-pis_latla)
    sabha.res<-bh.func(pv/ws_sabha, q)
    sabha.de<-sabha.res$de
    sabha.fdp[i, j]<-sum((1-theta)*sabha.de)/max(sum(sabha.de), 1)
    sabha.ntp[i, j]<-sum(theta*sabha.de)/sum(theta)
  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])

  latla.or.fdr[i]<-mean(latla.or.fdp[i,])
  latla.or.etp[i]<-mean(latla.or.ntp[i,])
  latla.dd.fdr[i]<-mean(latla.dd.fdp[i,])
  latla.dd.etp[i]<-mean(latla.dd.ntp[i,])

  wbh.fdr[i]<-mean(wbh.fdp[i,])
  wbh.etp[i]<-mean(wbh.ntp[i,])
  sabha.fdr[i]<-mean(sabha.fdp[i,])
  sabha.etp[i]<-mean(sabha.ntp[i,])
}
fdr_tf2.mthd<-cbind(bh.fdr, latla.or.fdr, latla.dd.fdr, sabha.fdr, wbh.fdr)
etp_tf2.mthd<-cbind(bh.etp, latla.or.etp, latla.dd.etp, sabha.etp, wbh.etp)



par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(mu_l.vec, fdr_tf2.mthd, type="o", pch=1:4, lwd=2, main="TF-setting2 FDR Comparison", xlab=expression(mu), ylab="FDR", ylim=c(0.0, 0.095))
legend("top", c("BH", "LATLA.OR","LATLA.DD", "SABHA", "WBH"), pch=1:4, col=1:6, lwd=2)

matplot(mu_l.vec, etp_tf2.mthd, type="o", pch=1:4, lwd=2, main="TF-setting2 Power Comparison", xlab=expression(mu), ylab="Power")

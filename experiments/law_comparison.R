source('G:/My Drive/LATLA/LATLA-code/latla_funcs.R')



m<-1000
mean.vec<-seq(from=3, to=3.5, by=0.1)
pis<-rep(0.1, m)
# pis[201:300]<-0.6
# pis[601:700]<-0.6
q<-0.05
np<-length(mean.vec)
nrep<-100

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

#------ overall the FDR and power-------------
law.or.fdr<-rep(0, np)
law.or.etp<-rep(0, np)
law.dd.fdr<-rep(0, np) 
law.dd.etp<-rep(0, np)

latla.or.fdr<-rep(0, np)
latla.or.etp<-rep(0, np)
latla.dd.fdr<-rep(0, np)
latla.dd.etp<-rep(0, np)


#------FDP and power for all exp-------------
bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)

law.or.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.or.ntp<-matrix(rep(0, nrep*np), np, nrep)
law.dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
law.dd.ntp<-matrix(rep(0, nrep*np), np, nrep)

latla.or.fdp<-matrix(rep(0, np*nrep), np, nrep)
latla.or.ntp<-matrix(rep(0, np*nrep), np, nrep)
latla.dd.fdp<-matrix(rep(0, np*nrep), np, nrep)
latla.dd.ntp<-matrix(rep(0, np*nrep), np, nrep)


for (i in 1:np)
{
  cat("\n", "iteration i= ", i, "\n", "repetition j=")
  mean<-mean.vec[i]
  
  for (j in 1:nrep)
  {
    set.seed(1000*i+j)
    cat(j)
    theta<-rbinom(m, size=1, prob=pis)
    pii<-sum(theta)/m
    mu0<-rep(0,m)
    mu1<-rep(mean,m)
    sd0<-rep(1,m)
    sd1<-rep(1,m)
    x0<-rnorm(m, mean=0, sd=1)
    x1<-rnorm(m, mean=mean, sd=1)
    x<-(1-theta)*x0+theta*x1
    pv<-2*pnorm(-abs(x), 0, 1)
    
    # generating the distance matrix
    d<-matrix(rep(0,m^2),m,m)
    for (k in 1:m) {
      d[k,]<-abs(1:m-k)
    }
    
    h = density(1:m)$bw
    pis_latla <- latla_pis(x, d=d, pval=pv, tau=bh.func(pv,0.8)$th, h=h, eps=0)
    # pis_latla <- latla_pis(x, d=d, pval=pv, tau=0.6, h=h, eps=0)
    
    qq=q/(1-pii)
    bh.res<-bh.func(pv, qq)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)
    
    wor_laws <- pis/(1-pis)
    law.or.res<-latla_thres(pvs=pv, pis=pis, ws=wor_laws, q)
    law.or.de<-law.or.res$de
    law.or.fdp[i, j]<-sum((1-theta)*law.or.de)/max(sum(law.or.de), 1)
    law.or.ntp[i, j]<-sum(theta*law.or.de)/sum(theta)
    
    pis_law <- pis_latla
    nu<-10e-5
    pis_law[which(pis_law<nu)]<-nu # stabilization
    pis_law[which(pis_law>1-nu)]<-1-nu # stabilization
    ws_laws<-pis_law/(1-pis_law)
    
    law.dd.res<-latla_thres(pvs=pv, pis=pis_law, ws=ws_laws, q)
    law.dd.de<-law.dd.res$de
    law.dd.fdp[i, j]<-sum((1-theta)*law.dd.de)/max(sum(law.dd.de), 1)
    law.dd.ntp[i, j]<-sum(theta*law.dd.de)/sum(theta)
    
    awor <- awor.func(x, pis, mu0, mu1, sd0, sd1, q)
    latla.or.res<-latla_thres(pvs=pv, pis=pis, ws=awor, q)
    latla.or.de<-latla.or.res$de
    latla.or.fdp[i, j]<-sum((1-theta)*latla.or.de)/max(sum(latla.or.de), 1)
    latla.or.ntp[i, j]<-sum(theta*latla.or.de)/sum(theta)

    weight<-latla_weights(x,d,pis_latla,mu0,sd0, h=h, eps=0)
    latla.dd.res<-latla_thres(pvs=pv, pis=pis_latla, ws=weight, q)
    latla.dd.de<-latla.dd.res$de
    latla.dd.fdp[i, j]<-sum((1-theta)*latla.dd.de)/max(sum(latla.dd.de), 1)
    latla.dd.ntp[i, j]<-sum(theta*latla.dd.de)/sum(theta)
  }
  
  bh.fdr[i]<-mean(bh.fdp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  law.or.fdr[i]<-mean(law.or.fdp[i,])
  law.or.etp[i]<-mean(law.or.ntp[i,])
  law.dd.fdr[i]<-mean(law.dd.fdp[i,])
  law.dd.etp[i]<-mean(law.dd.ntp[i,])
  latla.or.fdr[i]<-mean(latla.or.fdp[i,])
  latla.or.etp[i]<-mean(latla.or.ntp[i,])
  latla.dd.fdr[i]<-mean(latla.dd.fdp[i,])
  latla.dd.etp[i]<-mean(latla.dd.ntp[i,])
  
}

laws_new_fdr.mthd<-cbind(bh.fdr, latla.or.fdr, latla.dd.fdr, law.or.fdr, law.dd.fdr)
laws_new_etp.mthd<-cbind(bh.etp, latla.or.etp, latla.dd.etp, law.or.etp, law.dd.etp)


par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(2, 2, 2, 1)+0.1)

matplot(mean.vec, laws_new_fdr.mthd, type="o", pch=1:4, lwd=2, main="LAWS setting FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH",  "LATLA.OR", "LATLA.DD", "LAW.OR", "LAW.DD"), pch=1:3, col=1:6, lwd=2)

matplot(mean.vec, laws_new_etp.mthd, type="o", pch=1:4, lwd=2, main="LAWS setting Power Comparison", xlab=expression(sigma), ylab="ETP")

save(laws_new_fdr.mthd, file='./results/laws_comp_new_fdr.Rdata')
save(laws_new_etp.mthd, file='./results/laws_comp_new_etp.Rdata')
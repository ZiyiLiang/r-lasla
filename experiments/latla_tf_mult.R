source('C:/Users/liang/OneDrive/Desktop/LATLA/LATLA-code/latla_funcs.R')
## Simulation 1: 4 auxillary sequences, all obey common latent variable, compare with averaging auxillary data
##              only picking one sequence.


m<-1200
noise_l.vec<-seq(from=0.5, to=2, by=0.25)
pis<-rep(0.1, m)
q<-0.05
np<-length(noise_l.vec)
nrep<-50

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

latla.or.fdr<-rep(0, np)
latla.or.etp<-rep(0, np)
latla.dd.fdr<-rep(0, np)
latla.dd.etp<-rep(0, np)

avg.fdr<-rep(0, np)
avg.etp<-rep(0, np)


bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

latla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)

avg.fdp<-matrix(rep(0, np*nrep), nrep, np)
avg.ntp<-matrix(rep(0, np*nrep), nrep, np)



for (i in 1:nrep)
{
  cat("\n", "TF iteration i= ", i, "\n", "iteration j=")

  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(3,m)
  sd0<-rep(1,m)  #the auxillary source sd=0.1
  sd1<-rep(sqrt(1+1^2),m)
  #generate primary and auxillary data
  l0<-rep(0,m)
  l1<-rnorm(m, mean=3, sd=1)
  latent<-(1-theta)*l0+theta*l1
  x<-rnorm(m,mean=latent, sd=1)
  for (j in 1:np)
  {
    cat(j)
    noise<-noise_l.vec[j]
    s1<-rnorm(m, mean=latent, sd=noise)  # variation from source when collecting auxillary data
    s2<-rnorm(m, mean=latent, sd=noise)
    s3<-rnorm(m, mean=latent, sd=noise)
    s4<-rnorm(m, mean=latent, sd=noise)
    pv<-2*pnorm(-abs(x), 0, 1)

    # distance matrix for latla
    d_latla<-matrix(rep(0,m^2),m,m)
    S<-cbind(s1,s2,s3,s4)
    R<-cov(S)
    for (k in 1:m) {
      for (h in k:m) {
        d_latla[k,h]=mahalanobis(S[k,],S[h,],R)/5
      }
    }
    for(k in 2:m) {
      for (h in 1:(k-1)) {
        d_latla[k,h]=d_latla[h,k]
      }
    }

    # distance matrix for averaging method
    d_avg<-matrix(rep(0,m^2),m,m)
    s_avg<-(s1+s2+s3+s4)/4
    for (k in 1:m) {
      d_avg[k,]<-abs(s_avg-s_avg[k])
    }

    pis_latla<-latla_pis(x, d_latla, pval=pv, tau=bh.func(pv,0.8)$th)
    pis_avg<-latla_pis(x, d_avg, pval=pv, tau=bh.func(pv,0.8)$th)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    # using the asymmetricity adapted weights
    lpis<-latent_pis_mult.func(pis, s1,s2,s3,s4, mu1,noise)
    awor <- latent_awor_mult1.func(x,s1,s2,s3,s4,pis,mu1,noise,mu0,sd0)
    latla.or.res<-latla_thres(pvs=pv, lpis, ws=awor, q)
    latla.or.de<-latla.or.res$de
    latla.or.fdp[i, j]<-sum((1-theta)*latla.or.de)/max(sum(latla.or.de), 1)
    latla.or.ntp[i, j]<-sum(theta*latla.or.de)/sum(theta)

    weight<-latla_weights(x,d_latla,pis_latla,mu0,sd0)
    latla.dd.res<-latla_thres(pvs=pv, pis_latla, ws=weight, q)
    latla.dd.de<-latla.dd.res$de
    latla.dd.fdp[i, j]<-sum((1-theta)*latla.dd.de)/max(sum(latla.dd.de), 1)
    latla.dd.ntp[i, j]<-sum(theta*latla.dd.de)/sum(theta)

    weight1<-latla_weights(x,d_avg,pis_avg,mu0,sd0)
    avg.res<-latla_thres(pvs=pv, pis_avg, ws=weight1, q)
    avg.de<-avg.res$de
    avg.fdp[i, j]<-sum((1-theta)*avg.de)/max(sum(avg.de), 1)
    avg.ntp[i, j]<-sum(theta*avg.de)/sum(theta)
  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])
  latla.or.fdr[i]<-mean(latla.or.fdp[,i])
  latla.or.etp[i]<-mean(latla.or.ntp[,i])
  latla.dd.fdr[i]<-mean(latla.dd.fdp[,i])
  latla.dd.etp[i]<-mean(latla.dd.ntp[,i])
  avg.fdr[i]<-mean(avg.fdp[,i])
  avg.etp[i]<-mean(avg.ntp[,i])

}
fdr_mult1.mthd<-cbind(bh.fdr, latla.or.fdr, latla.dd.fdr, avg.fdr)
etp_mult1.mthd<-cbind(bh.etp, latla.or.etp, latla.dd.etp, avg.etp)


## Simulation 2: 4 auxillary sequences, two obey common latent variable with primary seq, two irrelevant auxillary
##               data. Compare with averaging method and randomly picking one.

m<-1200
noise_l.vec<-seq(from=0.5, to=2, by=0.25)
pis<-rep(0.1, m)
q<-0.05
np<-length(noise_l.vec)
nrep<-50

bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)
or.fdr<-rep(0, np)
or.etp<-rep(0, np)

latla.or.fdr<-rep(0, np)
latla.or.etp<-rep(0, np)
latla.dd.fdr<-rep(0, np)
latla.dd.etp<-rep(0, np)

avg.fdr<-rep(0, np)
avg.etp<-rep(0, np)


bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

or.fdp<-matrix(rep(0, np*nrep), nrep, np)
or.ntp<-matrix(rep(0, np*nrep), nrep, np)

latla.or.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.or.ntp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)
avg.fdp<-matrix(rep(0, np*nrep), nrep, np)
avg.ntp<-matrix(rep(0, np*nrep), nrep, np)



for (i in 1:nrep)
{
  cat("\n", "TF iteration i= ", i, "\n", "iteration j=")

  theta<-rbinom(m, size=1, prob=pis)
  dis_theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  mu0<-rep(0,m)
  mu1<-rep(3,m)
  sd0<-rep(1,m)  #the auxillary source sd=0.1
  sd1<-rep(sqrt(1+1^2),m)
  #generate primary and auxillary data
  l0<-rep(0,m)
  l1<-rnorm(m, mean=3, sd=1)
  dis_l0<-rnorm(m, mean=0, sd=1)
  dis_l1<-rnorm(m, mean=3, sd=1)
  latent<-(1-theta)*l0+theta*l1
  dis_latent<-(1-dis_theta)*dis_l0+dis_theta*dis_l1
  x<-rnorm(m,mean=latent, sd=1)
  for (j in 1:np)
  {
    cat(j)
    noise<-noise_l.vec[j]
    s1<-rnorm(m, mean=latent, sd=noise)  # variation from source when collecting auxillary data
    s2<-rnorm(m, mean=latent, sd=noise)
    s3<-rnorm(m, mean=dis_latent, sd=noise)
    s4<-rnorm(m, mean=dis_latent, sd=noise)
    pv<-2*pnorm(-abs(x), 0, 1)

    # distance matrix for latla
    d_latla<-matrix(rep(0,m^2),m,m)
    S<-cbind(s1,s2,s3,s4)
    R<-cov(S)
    for (k in 1:m) {
      for (h in k:m) {
        d_latla[k,h]=mahalanobis(S[k,],S[h,],R)/5
      }
    }
    for(k in 2:m) {
      for (h in 1:(k-1)) {
        d_latla[k,h]=d_latla[h,k]
      }
    }

    # distance matrix for averaging method
    d_avg<-matrix(rep(0,m^2),m,m)
    s_avg<-(s1+s2+s3+s4)/4
    for (k in 1:m) {
      d_avg[k,]<-abs(s_avg-s_avg[k])
    }

    pis_latla<-latla_pis(x, d_latla, pval=pv, tau=bh.func(pv,0.8)$th)
    pis_avg<-latla_pis(x, d_avg, pval=pv, tau=bh.func(pv,0.8)$th)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    lpis<-latent_pis_mult2.func(pis, s1,s2,s3,s4, mu1,noise)
    awor <- latent_awor_mult2.func(x,s1,s2,s3,s4,pis,mu1,noise,mu0,sd0)
    latla.or.res<-latla_thres(pvs=pv, lpis, ws=awor, q)
    latla.or.de<-latla.or.res$de
    latla.or.fdp[i, j]<-sum((1-theta)*latla.or.de)/max(sum(latla.or.de), 1)
    latla.or.ntp[i, j]<-sum(theta*latla.or.de)/sum(theta)

    weight<-latla_weights(x,d_latla,pis_latla,mu0,sd0)
    latla.dd.res<-latla_thres(pvs=pv, pis_latla, ws=weight, q)
    latla.dd.de<-latla.dd.res$de
    latla.dd.fdp[i, j]<-sum((1-theta)*latla.dd.de)/max(sum(latla.dd.de), 1)
    latla.dd.ntp[i, j]<-sum(theta*latla.dd.de)/sum(theta)

    weight1<-latla_weights(x,d_avg,pis_avg,mu0,sd0)
    avg.res<-latla_thres(pvs=pv, pis_avg, ws=weight1, q)
    avg.de<-avg.res$de
    avg.fdp[i, j]<-sum((1-theta)*avg.de)/max(sum(avg.de), 1)
    avg.ntp[i, j]<-sum(theta*avg.de)/sum(theta)

  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])
  latla.or.fdr[i]<-mean(latla.or.fdp[,i])
  latla.or.etp[i]<-mean(latla.or.ntp[,i])
  latla.dd.fdr[i]<-mean(latla.dd.fdp[,i])
  latla.dd.etp[i]<-mean(latla.dd.ntp[,i])
  avg.fdr[i]<-mean(avg.fdp[,i])
  avg.etp[i]<-mean(avg.ntp[,i])
}
fdr_mult2.mthd<-cbind(bh.fdr, latla.or.fdr, latla.dd.fdr, avg.fdr)
etp_mult2.mthd<-cbind(bh.etp, latla.or.etp, latla.dd.etp, avg.etp)


par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(noise_l.vec, fdr_mult1.mthd, type="o", pch=1:7, lwd=2, main="Mult-setting1 FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH", "LATLA.OR", "LATLA.DD","AVG"), pch=1:6, col=1:6, lwd=2)

matplot(noise_l.vec, etp_mult1.mthd, type="o", pch=1:7, lwd=2, main="Mult-setting1 Power Comparison", xlab=expression(sigma), ylab="Power")

matplot(noise_l.vec, fdr_mult2.mthd, type="o", pch=1:7, lwd=2, main="Mult-setting2 FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH", "LATLA.OR","LATLA.DD","AVG"), pch=1:6, col=1:6, lwd=2)

matplot(noise_l.vec, etp_mult2.mthd, type="o", pch=1:7, lwd=2, main="Mult-setting2 Power Comparison", xlab=expression(sigma), ylab="Power")


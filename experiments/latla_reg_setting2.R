#----------------------------------------------------------------------------------------
# Name: regression_setting2.R
# Description:
# Compare the performance of latla with other methods under the regression setting
#----------------------------------------------------------------------------------------

source("C:/Users/liang/OneDrive/Desktop/LATLA/LATLA-code/latla_funcs.R")

## Regression simulation 2: change the signal strength of the auxiliary sequences

## simple linear regression model
m <- 1000
p <- 800
pis<-rep(0.1,m)
mu_l.vec<-seq(from=0.25, to=0.35, by=0.025)
q <- 0.05
np <- length(mu_l.vec)
nrep <- 100


bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

latla.dd.fdr<-rep(0, np)
latla.dd.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

latla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
latla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)

for (i in 1:nrep)
{
  cat("\n", "TF iteration i= ", i, "\n", "iteration j=")

  theta<-rbinom(p, size=1, prob=pis)
  pii<-sum(theta)/p
  ## null distribution of the test statistics
  ## t distribution in ols, asymptotical standard normal in inverse regression
  mu0<-rep(0,p)
  sd0<-rep(1,p)

  for (j in 1:np){
    cat(j)
    mu <- mu_l.vec[j]
    # generate the real data
    x<-matrix(rnorm(m*p),m,p)
    beta<-matrix(rep(0,p),1,p)
    # assign real parameter value for non-null coefficients
    sign <- rbinom(sum(theta), size=1, prob=0.2)  # controls the symmetry
    beta[which(theta==1)] = abs(rnorm(sum(theta),mu,0.1))*(-1)^sign
    y<-x%*%t(beta)+rnorm(m)

    # calculate the t score and p value of the coefficients
    model=lm(y~x)
    tmp_se <- sqrt(diag(vcov(model)))[-1]
    tmp_beta <- model$coefficients[-1]
    t <- tmp_beta/tmp_se
    pv <- 2*pnorm(-abs(t))

    # construct the first auxiliary sequence
    beta1 <- beta + rnorm(p,0,0.05)
    y1 <- x%*%t(beta1)+rnorm(m)
    model1=lm(y1~x)
    tmp_se <- sqrt(diag(vcov(model1)))[-1]
    tmp_beta <- model1$coefficients[-1]
    t1 <- tmp_beta/tmp_se

    # construct the second auxiliary sequence
    beta2 <- beta + rnorm(p,0,0.05)
    y2 <- x%*%t(beta2)+rnorm(m)
    model2=lm(y2~x)
    tmp_se <- sqrt(diag(vcov(model2)))[-1]
    tmp_beta <- model2$coefficients[-1]
    t2 <- tmp_beta/tmp_se

    # construct the second auxiliary sequence
    beta3 <- beta + rnorm(p,0,0.05)
    y3 <- x%*%t(beta3)+rnorm(m)
    model3=lm(y3~x)
    tmp_se <- sqrt(diag(vcov(model3)))[-1]
    tmp_beta <- model3$coefficients[-1]
    t3 <- tmp_beta/tmp_se

    #distance matrix for latla
    d_latla<-matrix(rep(0,p^2),p,p)
    S<-cbind(t1,t2,t3)
    R<-cov(S)
    for (k in 1:p) {
      for (h in k:p) {
        d_latla[k,h]=mahalanobis(S[k,],S[h,],R)/4
      }
    }
    for(k in 2:p) {
      for (h in 1:(k-1)) {
        d_latla[k,h]=d_latla[h,k]
      }
    }

    pis_latla<-latla_pis(t, d_latla, pval=pv, tau=bh.func(pv,0.8)$th, eps=0)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    weight<-latla_weights(t,d_latla,pis_latla,mu0,sd0, eps=0)
    latla.dd.res<-latla_thres(pvs=pv, pis_latla, ws=weight, q)
    latla.dd.de<-latla.dd.res$de
    latla.dd.fdp[i, j]<-sum((1-theta)*latla.dd.de)/max(sum(latla.dd.de), 1)
    latla.dd.ntp[i, j]<-sum(theta*latla.dd.de)/sum(theta)

  }
}

for (i in 1:np) {
  bh.fdr[i]<-mean(bh.fdp[,i])
  bh.etp[i]<-mean(bh.ntp[,i])
  latla.dd.fdr[i]<-mean(latla.dd.fdp[,i])
  latla.dd.etp[i]<-mean(latla.dd.ntp[,i])

}
fdr_reg2.mthd<-cbind(bh.fdr, latla.dd.fdr)
etp_reg2.mthd<-cbind(bh.etp, latla.dd.etp)

par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(mu_l.vec, fdr_reg2.mthd, type="o", pch=1:7, lwd=2, main="Regression-setting2 FDR Comparison", xlab=expression(mu), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LATLA.DD"), pch=1:6, col=1:6, lwd=2)

matplot(mu_l.vec, etp_reg2.mthd, type="o", pch=1:7, lwd=2, main="Regression-setting2 Power Comparison", xlab=expression(mu), ylab="Power")


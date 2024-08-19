# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(adaptMT)

## Regression simulation 1: change the noise level of the auxiliary sequences
## simple linear regression model
m <- 1000
p <- 800
pis<-rep(0.1,m)
mean <- 0.25
noise_l.vec <- seq(from=0.1, to=0.2, by=0.02)
q <- 0.05
np <- length(noise_l.vec)
nrep <- 100


bh.fdr<-rep(0, np)
bh.etp<-rep(0, np)

lasla.dd.fdr<-rep(0, np)
lasla.dd.etp<-rep(0, np)

adapt.fdr<-rep(0, np)
adapt.etp<-rep(0, np)

bh.fdp<-matrix(rep(0, np*nrep), nrep, np)
bh.ntp<-matrix(rep(0, np*nrep), nrep, np)

lasla.dd.fdp<-matrix(rep(0, np*nrep), nrep, np)
lasla.dd.ntp<-matrix(rep(0, np*nrep), nrep, np)

adapt.fdp<-matrix(rep(0, np*nrep), nrep, np)
adapt.ntp<-matrix(rep(0, np*nrep), nrep, np)

pb <- progress_bar$new(total = nrep)   # show progress bar
for (i in 1:nrep)
{
  pb$tick()
  theta<-rbinom(p, size=1, prob=pis)
  mu0<-rep(0,p)
  sd0<-rep(1,p)

  x<-matrix(rnorm(m*p),m,p)   # Primary data
  beta<-matrix(rep(0,p),1,p)
  # assign real parameter value for non-null coefficients
  sign <- rbinom(sum(theta), size=1, prob=0.2)  # controls the symmetry
  beta[which(theta==1)] = abs(rnorm(sum(theta),mean=mean,sd=0.1))*(-1)^sign
  y<-x%*%t(beta)+rnorm(m)

  # Calculate the t score and p value of the coefficients
  model=lm(y~x)
  tmp_se <- sqrt(diag(vcov(model)))[-1]
  tmp_beta <- model$coefficients[-1]
  t <- tmp_beta/tmp_se
  pv <- 2*pnorm(-abs(t))

  for (j in 1:np){
    noise <- noise_l.vec[j]

    # Construct the first auxiliary sequence
    beta1 <- beta + rnorm(p,0,noise)
    y1 <- x%*%t(beta1)+rnorm(m)
    model1=lm(y1~x)
    tmp_se <- sqrt(diag(vcov(model1)))[-1]
    tmp_beta <- model1$coefficients[-1]
    t1 <- tmp_beta/tmp_se

    # Construct the second auxiliary sequence
    beta2 <- beta + rnorm(p,0,noise)
    y2 <- x%*%t(beta2)+rnorm(m)
    model2=lm(y2~x)
    tmp_se <- sqrt(diag(vcov(model2)))[-1]
    tmp_beta <- model2$coefficients[-1]
    t2 <- tmp_beta/tmp_se

    # Construct the third auxiliary sequence
    beta3 <- beta + rnorm(p,0,noise)
    y3 <- x%*%t(beta3)+rnorm(m)
    model3=lm(y3~x)
    tmp_se <- sqrt(diag(vcov(model3)))[-1]
    tmp_beta <- model3$coefficients[-1]
    t3 <- tmp_beta/tmp_se

    # Distance matrix
    d_lasla<-matrix(rep(0,p^2),p,p)
    S<-cbind(t1,t2, t3)
    R<-cov(S)

    for (k in 1:p) {
      d_lasla[k,]=mahalanobis(S,S[k,],R)
    }
    d_lasla[lower.tri(d_lasla)] <- t(d_lasla)[lower.tri(d_lasla)]
    d_lasla <- normalize_distance(t, d_lasla)
    

    pis_lasla<-lasla_pis(t, d_lasla, pval=pv, tau=bh.func(pv,0.8)$th)

    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum(theta*bh.de)/sum(theta)

    weight<-lasla_weights(t,d_lasla,pis_lasla,mu0,sd0,progress = FALSE)
    lasla.dd.res<-lasla_thres(pvs=pv, pis_lasla, ws=weight, q)
    lasla.dd.de<-lasla.dd.res$de
    lasla.dd.fdp[i, j]<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
    lasla.dd.ntp[i, j]<-sum(theta*lasla.dd.de)/sum(theta)
    
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
  lasla.dd.fdr[i]<-mean(lasla.dd.fdp[,i])
  lasla.dd.etp[i]<-mean(lasla.dd.ntp[,i])
  adapt.fdr[i]<-mean(adapt.fdp[,i])
  adapt.etp[i]<-mean(adapt.ntp[,i])
}
fdr_reg1.mthd<-cbind(bh.fdr, lasla.dd.fdr, adapt.fdr)
etp_reg1.mthd<-cbind(bh.etp, lasla.dd.etp, adapt.etp)



#######################################
#          Preview Results            #
#######################################
par(mfrow=c(2, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

matplot(noise_l.vec, fdr_reg1.mthd, type="o", pch=1:7, lwd=2, main="Regression-setting1 FDR Comparison", xlab=expression(sigma), ylab="FDR", ylim=c(0.01, 0.10))
legend("top", c("BH","LASLA.DD", "ADAPT"), pch=1:6, col=1:6, lwd=2)

matplot(noise_l.vec, etp_reg1.mthd, type="o", pch=1:7, lwd=2, main="Regression-setting1 Power Comparison", xlab=expression(sigma), ylab="Power")



#######################################
#             Save Results            #
#######################################
save = TRUE

if (save){
  data_dir <- "./results"
  
  method_names <- c("BH", "LASLA.DD", "ADAPT")
  
  # Setting 1
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i],length(noise_l.vec)),
                      FDR = fdr_reg1.mthd[,i],
                      Power = etp_reg1.mthd[,i],
                      Sigma = noise_l.vec)
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/regression1.RData", data_dir))
}
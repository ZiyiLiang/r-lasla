# Set working directory to the LASLA folder
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')


#######################################
#     Experiment Parameters           #
#######################################
m <- 1000
pis <- rep(0.1, m)
q <- 0.05
nrep <- 100
pb <- progress_bar$new(total = nrep)   # show progress bar



#------FDP and power for all exp-------------
bh.fdp<-rep(0, nrep)
bh.ntp<-rep(0, nrep)

lasla.or.fdp<-rep(0, nrep)
lasla.or.ntp<-rep(0, nrep)

wbh.fdp<-rep(0, nrep)
wbh.ntp<-rep(0, nrep)

awbh.fdp<-rep(0, nrep)
awbh.ntp<-rep(0, nrep)

for (i in 1:nrep)
{
  set.seed(i)
  
  pb$tick()
  pis <- c(runif(2*m/10,0.8,0.9), rep(0.01, 8*m/10))
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  gamma <- rep(1, m)
  mu1 <- rep(2, m)
  sd1 <- rep(1, m)
  
  # Generate heterogeneous alternative distribution
  isPositive<-rbinom(m, size=1, prob=gamma)
  x0<-rnorm(m, mean=0, sd=1)
  x1 <- rnorm(m, mean=mu1, sd=sd1)
  x1<-isPositive*x1 - (1-isPositive)*x1
  x<-(1-theta)*x0+theta*x1
  pv<-2*pnorm(-abs(x), 0, 1)
  
  qq=q/(1-pii)
  bh.res<-bh.func(pv, q)
  bh.de<-bh.res$de
  bh.fdp[i]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
  bh.ntp[i]<-sum(theta*bh.de)/sum(theta)
  
  wor_lasla <- lasla_oracle_weights.func(x, pis, gamma, dist_type="normal", q, mu1=mu1, sd1=sd1)
  lasla.or.res<-lasla_thres(pvs=pv, pis=pis, ws=wor_lasla, q)
  lasla.or.de<-lasla.or.res$de
  lasla.or.fdp[i]<-sum((1-theta)*lasla.or.de)/max(sum(lasla.or.de), 1)
  lasla.or.ntp[i]<-sum(theta*lasla.or.de)/sum(theta)
  
  sweight <- wor_lasla*m/sum(wor_lasla)
  wbh.res<-bh.func(pv/sweight, q)
  wbh.de<-wbh.res$de
  wbh.fdp[i]<-sum((1-theta)*wbh.de)/max(sum(wbh.de), 1)
  wbh.ntp[i]<-sum(theta*wbh.de)/sum(theta)
  
  ll <- bh.func(pv,0.8)$th
  phat <- (max(sweight) + sum(sweight * (pv > ll))) / (m * (1 - ll)) 
  awbh.res<-bh.func(pv/sweight, q/phat)
  awbh.de<-awbh.res$de
  awbh.fdp[i]<-sum((1-theta)*awbh.de)/max(sum(awbh.de), 1)
  awbh.ntp[i]<-sum(theta*awbh.de)/sum(theta)  
}


bh.fdr<-mean(bh.fdp)
bh.etp<-mean(bh.ntp)

wbh.fdr<-mean(wbh.fdp)
wbh.etp<-mean(wbh.ntp)

awbh.fdr<-mean(awbh.fdp)
awbh.etp<-mean(awbh.ntp)

lasla.or.fdr<-mean(lasla.or.fdp)
lasla.or.etp<-mean(lasla.or.ntp)

print(cbind(bh.fdr, lasla.or.fdr, wbh.fdr, awbh.fdr))
print(cbind(bh.etp, lasla.or.etp, wbh.etp, awbh.etp))

hetero_sparsity_fdp.mthd<-cbind(lasla.or.fdp, awbh.fdp, wbh.fdp)
hetero_sparsity_ntp.mthd<-cbind(lasla.or.ntp, awbh.ntp, wbh.ntp)



#######################################
#         Preview  Results            #
#######################################

par(mfrow=c(2, 2))  # 1 row, 2 columns

boxplot(hetero_sparsity_fdp.mthd,
        main = "Comparison of FDP Across Methods",
        ylab = "FDP",
        col = c("lightblue", "lightgreen", "lightpink"),
        names = c("LASLA", "adaptive-WBH", "WBH"))

boxplot(hetero_sparsity_ntp.mthd,
        main = "Comparison of NTP Across Methods",
        ylab = "NTP",
        col = c("lightblue", "lightgreen", "lightpink"),
        names = c("LASLA", "adaptive-WBH", "WBH"))


#######################################
#             Save Results            #
#######################################
save <- TRUE

if(save){
  data_dir <- "./results"
  
  method_names <- c("LASLA.OR", "AWBH", "WBH")
  
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i], nrep),
                      rep = seq(1, nrep, 1),
                      FDP = hetero_sparsity_fdp.mthd[,i],
                      Power = hetero_sparsity_ntp.mthd[,i])
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/heterogeneous_sparsity.RData", data_dir))
}
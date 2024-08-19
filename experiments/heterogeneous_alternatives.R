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

law.or.fdp<-rep(0, nrep)
law.or.ntp<-rep(0, nrep)

lasla.or.fdp<-rep(0, nrep)
lasla.or.ntp<-rep(0, nrep)


for (i in 1:nrep)
{
  set.seed(i)
  
  pb$tick()
  pis <- rep(0.1, m)
  theta<-rbinom(m, size=1, prob=pis)
  pii<-sum(theta)/m
  gamma <- runif(m,0,1)
  mu1 <- runif(m, 2.5, 3.5)
  sd1 <- runif(m,0.1,1)

  # Generate heterogeneous alternative distribution
  isPositive<-rbinom(m, size=1, prob=gamma)
  x0<-rnorm(m, mean=0, sd=1)
  x1 <- rnorm(m, mean=mu1, sd=sd1)
  x1<-isPositive*x1 - (1-isPositive)*x1
  x<-(1-theta)*x0+theta*x1
  pv<-2*pnorm(-abs(x), 0, 1)
  
  
  qq=q/(1-pii)
  bh.res<-bh.func(pv, qq)
  bh.de<-bh.res$de
  bh.fdp[i]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
  bh.ntp[i]<-sum(theta*bh.de)/sum(theta)
  
  wor_laws <- pis/(1-pis)
  law.or.res<-lasla_thres(pvs=pv, pis=pis, ws=wor_laws, q)
  law.or.de<-law.or.res$de
  law.or.fdp[i]<-sum((1-theta)*law.or.de)/max(sum(law.or.de), 1)
  law.or.ntp[i]<-sum(theta*law.or.de)/sum(theta)
  
  wor_lasla <- lasla_oracle_weights.func(x, pis, gamma, dist_type="normal", q, mu1=mu1, sd1=sd1)
  lasla.or.res<-lasla_thres(pvs=pv, pis=pis, ws=wor_lasla, q)
  lasla.or.de<-lasla.or.res$de
  lasla.or.fdp[i]<-sum((1-theta)*lasla.or.de)/max(sum(lasla.or.de), 1)
  lasla.or.ntp[i]<-sum(theta*lasla.or.de)/sum(theta)
}


bh.fdr<-mean(bh.fdp)
bh.etp<-mean(bh.ntp)

lasla.or.fdr<-mean(lasla.or.fdp)
lasla.or.etp<-mean(lasla.or.ntp)

law.or.fdr<-mean(law.or.fdp)
law.or.etp<-mean(law.or.ntp)

print(cbind(bh.fdr, lasla.or.fdr, law.or.fdr))
print(cbind(bh.etp, lasla.or.etp, law.or.etp))

altshape_fdp.mthd<-cbind(bh.fdp, lasla.or.fdp, law.or.fdp)
altshape_ntp.mthd<-cbind(bh.ntp, lasla.or.ntp, law.or.ntp)




#######################################
#             Save Results            #
#######################################
save <- TRUE

if(save){
  data_dir <- "./results"
  
  method_names <- c("PV.OR", "LASLA.OR", "LAWS.OR")
  
  results <- data.frame()
  
  for (i in 1:length(method_names)) {
    tmp <- data.frame(Method = rep(method_names[i], nrep),
                      rep = seq(1, nrep, 1),
                      FDP = altshape_fdp.mthd[,i],
                      Power = altshape_ntp.mthd[,i])
    results <- rbind(results, tmp)
  }
  save(results, file=sprintf("%s/heterogeneous_distribution.RData", data_dir))
}
#############################################################
# File contains oracle functions or other utility functions #
#############################################################

bh.func<-function(pv, q)
{
  # the input
  # pv: the p-values
  # q: the FDR level
  # the output
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # de: the decision rule
  
  m=length(pv)
  st.pv<-sort(pv)
  pvi<-st.pv/1:m
  de<-rep(0, m)
  if (sum(pvi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    de[which(pv<=pk)]<-1
  }
  y<-list(nr=k, th=pk, de=de)
  return (y)
}


latent_pis.func <-function(pis, s, mu1, sigma)
{
  ## calculates the conditional sparsity level of the latent variable setting
  ## Arguments
  # pis: sparsity level
  # s: the auxiliary sequence
  # sigma: noise of the auxiliary sequence
  ## Values
  # lpis: the conditional sparsity level
  lpis <- rep(0,m)
  num <- (1-pis)*dnorm(s,0,sigma)
  den <- pis*dnorm(s,mu1,sqrt(1+sigma^2))
  lpis <- num/(num+den)
  return(1-lpis)
}

######## other oracle related functions used in simulation ##########
latent_awor.func <- function(x, s, pis, mu1, sigma, mu0, sd0, q=0.05){
  ## Arguments
  # x: vector of normal variables
  # pis: conditional probabilities
  # mu1: signal strength
  # sigma: noise level of the auxiliary sequence
  # mu0: mean vector of null distribution
  # sd0: sd vector of null distribution
  # q: fdr level
  # Values
  # awor: the asymmetricity adapted oracle weights
  m <- length(x)
  tor <- rep(0,m)
  # calculating the oracle statistics
  lpis <- rep(0,m)
  num_s <- (1-pis)*dnorm(s,0,sigma)
  den_s <- pis*dnorm(s,mu1,sqrt(1+sigma^2))
  lpis <- 1-(num_s/(num_s+den_s))
  # calculating the oracle statistics
  num <- (1-lpis)*dnorm(x,mu0,sd0)
  den <- pis*exp(-((1+sigma^2)*x*2+(-2*s-2*mu1*sigma^2)*x+2*s^2-2*mu1*s+mu1^2*sigma^2+mu1^2)/(4*sigma^2+2))/(2*pi*sqrt(1+2*sigma^2))
  tor <- num/(num+den/(num_s+den_s))
  
  # calculating the moving average
  st.tor <- sort(tor)
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i
  }
  # calculating the oracle threshold
  if (sum(mva<=q)==0){
    th <- -Inf
  }
  else{
    th <- max(which(mva<=q))
    th <- st.tor[th]
  }
  
  
  # calculating the weights
  awor <- rep(0,m)
  nt <- seq(min(x)-2, 0, 0.05) # negative searching window
  pt <- seq(0, max(x)+2, 0.05) # positive searching window
  t_star <- rep(0,m)
  for (i in 1:m){
    if(x[i]>=0){
      temp_num <- (1-lpis[i])*dnorm(pt,mu0[i],sd0[i])
      temp_den <- pis[i]*exp(-((1+sigma^2)*pt*2+(-2*s[i]-2*mu1[i]*sigma^2)*pt+2*s[i]^2-2*mu1[i]*s[i]+mu1[i]^2*sigma^2+mu1[i]^2)/(4*sigma^2+2))/(2*pi*sqrt(1+2*sigma^2))
      t_temp <- temp_num/(temp_num+temp_den/((1-pis[i])*dnorm(s[i],0,sigma)+pis[i]*dnorm(s[i],mu1[i],sqrt(1+sigma^2))))
      if(length(which(t_temp<=th))==0){
        t_star[i] <- max(x)+2
      }
      else{
        t_star[i] = pt[min(which(t_temp<=th))]
      }
      awor[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      temp_num <- (1-lpis[i])*dnorm(nt,mu0[i],sd0[i])
      temp_den <- pis[i]*exp(-((1+sigma^2)*nt*2+(-2*s[i]-2*mu1[i]*sigma^2)*nt+2*s[i]^2-2*mu1[i]*s[i]+mu1[i]^2*sigma^2+mu1[i]^2)/(4*sigma^2+2))/(2*pi*sqrt(1+2*sigma^2))
      t_temp <- temp_num/(temp_num+temp_den/((1-pis[i])*dnorm(s[i],0,sigma)+pis[i]*dnorm(s[i],mu1[i],sqrt(1+sigma^2))))
      if(length(which(t_temp<=th))==0){
        t_star[i] <- min(x)-2
      }
      else{
        t_star[i] = nt[max(which(t_temp<=th))]
      }
      awor[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  #awor <- awor*m/sum(awor)
  nu<-10e-5
  awor[which(awor<nu)]<-nu # stabilization
  return(awor)
}


#####assumed 0,1 as mu0 and sd0, edit if needed.
latent_awor_mult1.func <- function(x, s1,s2,s3,s4, pis, mu1, sigma, mu0, sd0, q=0.05){
  ## calculates the oracle weights for multiple seq latent setting 1
  ## Arguments
  # x: vector of normal variables
  # s1-s4: auxillary sequences
  # pis: conditional probabilities
  # mu1: signal strength
  # sigma: noise level of the auxiliary sequence
  # mu0: mean vector of null distribution
  # sd0: sd vector of null distribution
  # q: fdr level
  # Values
  # awor: the asymmetricity adapted oracle weights
  m <- length(x)
  tor <- rep(0,m)
  # calculating the oracle statistics
  num<-rep(0,m)
  den<-rep(0,m)
  num<-10000*(1-pis)*dnorm(x,0,1)*dnorm(s1,0,sigma)*dnorm(s2,0,sigma)*dnorm(s3,0,sigma)*dnorm(s4,0,sigma)
  for (i in 1:length(pis)) {
    int<-function(y){10000*dnorm(x[i],y,1)*dnorm(s1[i],y,sigma)*dnorm(s2[i],y,sigma)*dnorm(s3[i],y,sigma)*dnorm(s4[i],y,sigma)*pis[i]*dnorm(y,mu1[i],1)}
    den[i]<-integrate(int,-Inf,Inf)$value
  }
  tor<-num/(num+den)
  # calculating the moving average
  st.tor <- sort(tor)
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]
  
  # calculating the weights
  awor <- rep(0,m)
  nt <- seq(min(x), 0, 0.05) # negative searching window
  pt <- seq(0, max(x), 0.05) # postive searching window
  t_star <- rep(0,m)
  for (i in 1:m){
    if(x[i]>=0){
      num<-rep(0,length(pt))
      den<-rep(0,length(pt))
      num<-10000*(1-pis[i])*dnorm(pt,0,1)*dnorm(s1[i],0,sigma)*dnorm(s2[i],0,sigma)*dnorm(s3[i],0,sigma)*dnorm(s4[i],0,sigma)
      for (j in 1:length(pt)) {
        int<-function(y){10000*dnorm(pt[j],y,1)*dnorm(s1[i],y,sigma)*dnorm(s2[i],y,sigma)*dnorm(s3[i],y,sigma)*dnorm(s4[i],y,sigma)*pis[i]*dnorm(y,mu1[i],1)}
        den[j]<-integrate(int,-Inf,Inf)$value
      }
      t_temp<-num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- max(x)
      }
      else{
        t_star[i] = pt[min(which(t_temp<=th))]
      }
      awor[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      num<-rep(0,length(nt))
      den<-rep(0,length(nt))
      num<-10000*(1-pis[i])*dnorm(nt,0,1)*dnorm(s1[i],0,sigma)*dnorm(s2[i],0,sigma)*dnorm(s3[i],0,sigma)*dnorm(s4[i],0,sigma)
      for (j in 1:length(nt)) {
        int<-function(y){10000*dnorm(nt[j],y,1)*dnorm(s1[i],y,sigma)*dnorm(s2[i],y,sigma)*dnorm(s3[i],y,sigma)*dnorm(s4[i],y,sigma)*pis[i]*dnorm(y,mu1[i],1)}
        den[j]<-integrate(int,-Inf,Inf)$value
      }
      t_temp<-num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- min(x)
      }
      else{
        t_star[i] = nt[max(which(t_temp<=th))]
      }
      awor[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  awor <- awor*m/sum(awor)
  nu<-10e-5
  awor[which(awor<nu)]<-nu # stabilization
  return(awor)
}

latent_awor_mult2.func <- function(x, s1,s2,s3,s4, pis, mu1, sigma, mu0, sd0, q=0.05){
  ## calculates the oracle weights for multiple latent setting 2
  ## Arguments
  # x: vector of normal variables
  # s1-s4: auxiliary sequences
  # pis: conditional probabilities
  # mu1: signal strength
  # sigma: noise level of the auxiliary sequence
  # mu0: mean vector of null distribution
  # sd0: sd vector of null distribution
  # q: fdr level
  # Values
  # awor: the asymmetricity adapted oracle weights
  m <- length(x)
  tor <- rep(0,m)
  # calculating the oracle statistics
  num<-rep(0,m)
  den<-rep(0,m)
  num<-10000*(1-pis)*dnorm(x,0,1)*dnorm(s1,0,sigma)*dnorm(s2,0,sigma)
  for (i in 1:length(pis)) {
    int<-function(y){10000*dnorm(x[i],y,1)*dnorm(s1[i],y,sigma)*dnorm(s2[i],y,sigma)*pis[i]*dnorm(y,mu1[i],1)}
    den[i]<-integrate(int,-Inf,Inf)$value
  }
  tor<-num/(num+den)
  # calculating the moving average
  st.tor <- sort(tor)
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]
  
  # calculating the weights
  awor <- rep(0,m)
  nt <- seq(min(x), 0, 0.05) # negative searching window
  pt <- seq(0, max(x), 0.05) # postive searching window
  t_star <- rep(0,m)
  for (i in 1:m){
    if(x[i]>=0){
      num<-rep(0,length(pt))
      den<-rep(0,length(pt))
      num<-10000*(1-pis[i])*dnorm(pt,0,1)*dnorm(s1[i],0,sigma)*dnorm(s2[i],0,sigma)
      for (j in 1:length(pt)) {
        int<-function(y){10000*dnorm(pt[j],y,1)*dnorm(s1[i],y,sigma)*dnorm(s2[i],y,sigma)*pis[i]*dnorm(y,mu1[i],1)}
        den[j]<-integrate(int,-Inf,Inf)$value
      }
      t_temp<-num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- max(x)
      }
      else{
        t_star[i] = pt[min(which(t_temp<=th))]
      }
      awor[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      num<-rep(0,length(nt))
      den<-rep(0,length(nt))
      num<-10000*(1-pis[i])*dnorm(nt,0,1)*dnorm(s1[i],0,sigma)*dnorm(s2[i],0,sigma)
      for (j in 1:length(nt)) {
        int<-function(y){10000*dnorm(nt[j],y,1)*dnorm(s1[i],y,sigma)*dnorm(s2[i],y,sigma)*pis[i]*dnorm(y,mu1[i],1)}
        den[j]<-integrate(int,-Inf,Inf)$value
      }
      t_temp<-num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- min(x)
      }
      else{
        t_star[i] = nt[max(which(t_temp<=th))]
      }
      awor[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  awor <- awor*m/sum(awor)
  nu<-10e-5
  awor[which(awor<nu)]<-nu # stabilization
  return(awor)
}

latent_pis_mult.func <-function(pis, s1,s2,s3,s4, mu1, sigma)
{
  ## calculates the conditional sparsity level of the latent variable setting
  ## Arguments
  # pis: sparsity level
  # s: the auxiliary sequence
  # sigma: noise of the auxiliary sequence
  ## Values
  # lpis: the conditional sparsity level
  lpis<-rep(0,m)
  num<-(1-pis)*dnorm(s1,0,sigma)*dnorm(s2,0,sigma)*dnorm(s3,0,sigma)*dnorm(s4,0,sigma)
  den<-rep(0,m)
  for (i in 1:length(pis)) {
    int<-function(x){dnorm(s1[i],x,sigma)*dnorm(s2[i],x,sigma)*dnorm(s3[i],x,sigma)*dnorm(s4[i],x,sigma)*pis[i]*dnorm(x,mu1[i],1)}
    den[i]<-integrate(int,-Inf,Inf)$value
  }
  lpis<-num/(num+den)
  return(1-lpis)
}

latent_pis_mult2.func <-function(pis, s1,s2,s3,s4, mu1, sigma)
{
  ## calculates the conditional sparsity level of the latent variable setting
  ## Arguments
  # pis: sparsity level
  # s: the auxiliary sequence
  # sigma: noise of the auxiliary sequence
  ## Values
  # lpis: the conditional sparsity level
  lpis<-rep(0,m)
  num<-(1-pis)*dnorm(s1,0,sigma)*dnorm(s2,0,sigma)
  den<-rep(0,m)
  for (i in 1:length(pis)) {
    int<-function(x){dnorm(s1[i],x,sigma)*dnorm(s2[i],x,sigma)*pis[i]*dnorm(x,mu1[i],1)}
    den[i]<-integrate(int,-Inf,Inf)$value
  }
  lpis<-num/(num+den)
  return(1-lpis)
}

awor.func <- function(x, pis, mu0, mu1, sd0, sd1, q=0.05){
  ## awor_1D.func improves the weight by adapting asymmetricity
  ## Arguments
  # x: vector of normal variables
  # pis: conditional probabilities
  # mu0: null mean vector
  # mu1: alternative mean vector
  # sd0: null sd vector
  # sd1: alternative sd vector
  # q: fdr level
  # Values
  # awor: the asymmetricity adapted oracle weights
  m <- length(x)
  tor <- rep(0,m)
  # calculating the oracle statistics
  tor <- (1-pis)*dnorm(x,mu0,sd0)/((1-pis)*dnorm(x,mu0,sd0)+pis*dnorm(x,mu1,sd1))
  
  # calculating the moving average
  st.tor <- sort(tor)
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i 
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]
  
  # calculating the weights
  awor <- rep(0,m)
  nt <- seq(min(x)-2, 0, 0.05) # negative searching window
  pt <- seq(0, max(x)+2, 0.05) # postive searching window
  t_star <- rep(0,m)
  for (i in 1:m){
    if(x[i]>=0){
      t_temp <- (1-pis[i])*dnorm(pt,mu0[i],sd0[i])/
        ((1-pis[i])*dnorm(pt,mu0[i],sd0[i])+pis[i]*dnorm(pt,mu1[i],sd1[i]))
      if(length(which(t_temp<=th))==0){
        t_star[i] <- max(x)+2
      }
      else{
        t_star[i] = pt[min(which(t_temp<=th))]
      }
      awor[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      t_temp <- (1-pis[i])*dnorm(nt,mu0[i],sd0[i])/
        ((1-pis[i])*dnorm(nt,mu0[i],sd0[i])+pis[i]*dnorm(nt,mu1[i],sd1[i]))
      if(length(which(t_temp<=th))==0){
        t_star[i] <- min(x)-2
      }
      else{
        t_star[i] = nt[max(which(t_temp<=th))]
      }
      awor[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  #awor <- awor*m/sum(awor)
  nu<-10e-5
  awor[which(awor<nu)]<-nu # stabilization
  return(awor)
}


awor_asymmetric.func <- function(x, pis, gamma, mu1, sd1, q=0.05){
  ## oracle lasla weight function for asymmetric setting with gamma 
  # controlling the level of asymmetry.
  # x: vector of normal variables
  # pis: conditional probabilities
  # gamma: level of asymmetry
  # mu1: alternative mean vector
  # sd1: alternative sd vector
  # q: fdr level
  # Values
  # awor: the asymmetricity adapted oracle weights
  m <- length(x)
  
  # Null distribution is assumed to be standard normal
  mu0 <- rep(0, m)
  sd0 <- rep(1, m)
  
  tor <- rep(0,m)
  # calculating the oracle statistics
  num <- (1-pis)*dnorm(x,mu0,sd0)
  den <- pis*(gamma*dnorm(x,mu1,sd1)+(1-gamma)*dnorm(x,-mu1,sd1))
  tor <- num/(num + den)
  
  # calculating the moving average
  st.tor <- sort(tor)
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i 
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]
  
  # calculating the weights
  awor <- rep(0,m)
  nt <- seq(min(x)-2, 0, 0.05) # negative searching window
  pt <- seq(0, max(x)+2, 0.05) # postive searching window
  t_star <- rep(0,m)
  for (i in 1:m){
    if(x[i]>=0){
      num <- (1-pis[i])*dnorm(pt,mu0[i],sd0[i])
      den <- pis[i]*(gamma[i]*dnorm(pt,mu1[i],sd1[i])+(1-gamma[i])*dnorm(pt,-mu1[i],sd1[i]))
      t_temp <- num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- max(x)+2
      }
      else{
        t_star[i] = pt[min(which(t_temp<=th))]
      }
      awor[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      num <- (1-pis[i])*dnorm(nt,mu0[i],sd0[i])
      den <- pis[i]*(gamma[i]*dnorm(nt,mu1[i],sd1[i])+(1-gamma[i])*dnorm(nt,-mu1[i],sd1[i]))
      t_temp <- num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- min(x)-2
      }
      else{
        t_star[i] = nt[max(which(t_temp<=th))]
      }
      awor[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  #awor <- awor*m/sum(awor)
  nu<-10e-5
  awor[which(awor<nu)]<-nu # stabilization
  return(awor)
}


lasla_oracle_weights.func <- function(x, pis, gamma, dist_type, q = 0.05, ...){
  # This function computes the oracle lasla weights 
  #
  # Args:
  # x: vector of normal variables
  # pis: conditional probabilities
  # dist_type: type of the alternative distribution
  # q: fdr level
  #
  # Output:
  # weights: the oracle weights
  
  m <- length(x)
  
  # Null distribution is assumed to be standard normal
  mu0 <- rep(0, m)
  sd0 <- rep(1, m)
  
  # Initialize variables for distribution parameters
  mu1 <- sd1 <- loc <- scale <- shape <- NULL
  
  # Parse additional arguments
  args <- list(...)
  
  # Handle different distribution types
  if (dist_type == "normal") {
    mu1 <- args$mu1
    sd1 <- args$sd1
    if (is.null(mu1) || is.null(sd1)) {
      stop("For normal distribution, 'mu1' and 'sd1' must be provided.")
    }
  } else if (dist_type == "skewed_normal") {
    loc <- args$loc
    scale <- args$scale
    shape <- args$shape
    if (is.null(loc) || is.null(scale) || is.null(shape)) {
      stop("For skewed normal distribution, 'loc', 'scale', and 'shape' must be provided.")
    }
  } else {
    stop("Unsupported distribution type. Use 'normal' or 'skewed_normal'.")
  }
  
  
  # Compute the oracle threshold
  tor <- rep(0,m)
  # calculating the oracle statistics
  num <- (1-pis)*dnorm(x,mu0,sd0)
  
  if (dist_type == "normal") {
    den <- pis*(gamma*dnorm(x,mu1,sd1)+(1-gamma)*dnorm(x,-mu1,sd1))
  }
  else if (dist_type == "skewed_normal") {
    den <- pis*(gamma*dsn(x,loc, scale, shape) + (1-gamma)*dsn(x,-loc, scale, shape))
  }
  tor <- num/(num + den)
  
  # calculating the moving average
  st.tor <- sort(tor)
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i 
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]
  
  # calculating the weights
  weights <- rep(0,m)
  nt <- seq(min(x)-2, 0, 0.05) # negative searching window
  pt <- seq(0, max(x)+2, 0.05) # postive searching window
  t_star <- rep(0,m)
  for (i in 1:m){
    if(x[i]>=0){
      num <- (1-pis[i])*dnorm(pt,mu0[i],sd0[i])
      if (dist_type == "normal") {
        den <- pis[i]*(gamma[i]*dnorm(pt,mu1[i],sd1[i])+(1-gamma[i])*dnorm(pt,-mu1[i],sd1[i]))
      }
      else if (dist_type == "skewed_normal") {
        den <- pis[i]*(gamma[i]*dsn(pt,loc[i], scale[i], shape[i]) + (1-gamma[i])*dsn(pt,-loc[i], scale[i], shape[i]))
      }
      
      t_temp <- num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- max(x)+2
      }
      else{
        t_star[i] = pt[min(which(t_temp<=th))]
      }
      weights[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      num <- (1-pis[i])*dnorm(nt,mu0[i],sd0[i])
      if (dist_type == "normal") {
        den <- pis[i]*(gamma[i]*dnorm(nt,mu1[i],sd1[i])+(1-gamma[i])*dnorm(nt,-mu1[i],sd1[i]))
      }
      else if (dist_type == "skewed_normal") {
        den <- pis[i]*(gamma[i]*dsn(nt,loc[i], scale[i], shape[i]) + (1-gamma[i])*dsn(nt,-loc[i], scale[i], shape[i]))
      }
      
      t_temp <- num/(num+den)
      if(length(which(t_temp<=th))==0){
        t_star[i] <- min(x)-2
      }
      else{
        t_star[i] = nt[max(which(t_temp<=th))]
      }
      weights[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  nu<-10e-5
  weights[which(weights<nu)]<-nu # stabilization
  return(weights)
}

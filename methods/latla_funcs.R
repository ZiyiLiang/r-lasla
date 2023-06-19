############ functions for data-driven latla ##################
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


latla_thres<-function(pvs, pis, ws, q)
{
  ## implementing neat
  ## Arguments
  # pvs: p-values
  # pis: conditional probabilities
  # ws: weights proposed by neat
  # q: FDR level
  ## Values
  # de: the decision
  # th: the threshold for weighted p-values

  m<-length(pvs)
  pws<-pvs/ws
  st.pws<-sort(pws)
  fdps<-sum(ws*(1-pis))*st.pws/(1:m)
  de<-rep(0, m)
  if(sum(fdps<=q)==0)
  {
    k<-0
    pwk<-1
  }
  else
  {
    k<-max(which(fdps<=q))
    pwk<-st.pws[k]
    de[which(pws<=pwk)]<-1
  }
  y<-list(nr=k, th=pwk, de=de)
  return (y)
}


############## Data-driven Weights #####################################
latla_weights <- function(x, d, pis, mu0, sd0, q=0.05, h='auto', eps=0.1){
  ## data-driven function that estimates the weights proposed by ANEAT
  ## Arguments:
  # x: vector of normal variables
  # d: the distance matrix
  # pis: the estimated sparsity level
  # mu0: mean vector of null
  # sd0: sd vector of null
  # q: fdr level
  # h: bandwith for the kernel estimation, default value is auto-selected
  #    by 'ste' solve the equation in
  # eps: determines the number of neighbors to be used
  ## Outputs:
  # wneat: the oracle-assited weights

  m <- length(x)
  # check default argument of nb (number of neighbors to be used)
  nb <- floor(m^(1-eps))
  # check default arugument of bandwidth
  if (h=='auto'){
    ds=density(d[1,],bw="SJ-ste")
    h=ds$bw
  }
  zds=density(x,from = 2*min(x),to=2*max(x)+5,n=1000)
  hx=zds$bw
  
  tor.est <- rep(0, m)
  for (i in 1:m) {
    vhi <- dnorm(d[i,], 0, h)
    vhi[i] = 0
    kht <- dnorm(x-x[i], 0, hx)
    # estimated conditional density evaluated at ti
    fti <- sum(vhi*kht)/sum(vhi)
    tor.est[i] <- (1-pis[i])*dnorm(x[i],mu0[i],sd0[i])/fti
  }
  st.tor <- sort(tor.est)

  # calculating the moving average
  mva <- rep(0,m)
  for (i in 1:m) {
    mva[i] <- sum(st.tor[1:i])/i
  }
  # calculating the oracle threshold
  th <- max(which(mva<=q))
  th <- st.tor[th]

  # calculating the weights
  weights <- rep(0,m)
  pt <- seq(0, max(x)+2, 0.05)
  nt <- seq(min(x)-2, 0, 0.05)
  t_star <- rep(0,m)
  n_p <- length(pt)
  n_n <- length(nt)
  for(i in 1:m){
    vhi <- dnorm(d[i,], 0, h)
    vhi[i]=0
    # only consider the test points in the neighborhood
    if (nb < m){
      st.nb<- sort(vhi)[m-nb]
      vhi[which(vhi<st.nb)]<-0
    }
    # adapt the asymmetricity
    if(x[i]>=0){
      # ft <- rep(0,n_p)
      # tor <- rep(0,n_p)
      pt_mat <- matrix(rep(pt,m),m,n_p,byrow = TRUE)
      x_mat <- matrix(rep(x,n_p),m,n_p)
      kht <- dnorm(x_mat - pt_mat, 0, hx)
      ft <- c(vhi %*% kht)/sum(vhi)
      # for (j in 1:n_p) {
      #   kht <- dnorm(x-pt[j],0,h)
      #   ft[j] <- sum(vhi*kht)/sum(vhi)
      # }
      tor<-(1-pis[i])*dnorm(pt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- max(x)+1
      }
      else{
        t_star[i] = pt[min(which(tor<=th))]
      }
      weights[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      # ft <- rep(0,n_n)
      # tor <- rep(0,n_n)
      nt_mat <- matrix(rep(nt,m),m,n_n,byrow = TRUE)
      x_mat <- matrix(rep(x,n_n),m,n_n)
      kht <- dnorm(x_mat - nt_mat, 0, hx)
      ft <- c(vhi %*% kht)/sum(vhi)
      # for (j in 1:n_n) {
      #   kht <- dnorm(x-nt[j],0,h)
      #   ft[j] <- sum(vhi*kht)/sum(vhi)
      # }
      tor<-(1-pis[i])*dnorm(nt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- min(x)-1
      }
      else{
        t_star[i] = nt[max(which(tor<=th))]
      }
      weights[i] <- pnorm(t_star[i], mu0[i], sd0[i])
    }
  }
  nu<-10e-5
  weights[which(weights<nu)]<-nu # stabilization
  return(weights)
}

latla_pis <- function(x, d, pval, tau=0.9, h='auto', eps=0.1)
{
  ## latla_pis calculates the conditional proportions by defining the null set using p values
  ## Arguments:
  # x: vector of normal variables
  # d: distance matrix
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # h: bandwidth for kernel estimation, default is automatically chosen by method 'ste'.
  # eps: determines the number of neighbors to be used
  ## Outputs:
  # pis_latla: the conditional proportions

  m <- length(x)
  # check default argument of nb (number of neighbors to be used)
  nb <- floor(m^(1-eps))
  if (h=='auto'){
    zds=density(d[1,],bw="SJ-ste")
    h=zds$bw
  }

  p.est <- rep(0, m)
  for (i in 1:m){
    kht <- dnorm(d[i,], 0, h)
    kht[i] = 0 # considers only j!=i
    if (nb < m){
      st.nb<- sort(kht)[m-nb]
      kht[which(kht<st.nb)]<-0
    }
    p.est[i] <- sum(kht[which(pval>=tau)])/((1-tau)*sum(kht))
  }
  p.est[which(p.est>1)] <- 1
  return(1-p.est)
}

######## other oracle related functions used in simulation ##########
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
  th <- max(which(mva<=q))
  th <- st.tor[th]

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


awor_asymmetric.func <- function(x, pis, gamma, mu0, mu1, sd0, sd1, q=0.05){
  ## oracle LATLA weight function for asymmetric setting with gamma 
  # controlling the level of asymmetry.
  # x: vector of normal variables
  # pis: conditional probabilities
  # gamma: level of asymmetry
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


#----------------------test func to be deleted
pis_1D_law.func<- function(x, tau=0.1, h=50)
{
  ## pis_est_law.func calculates the conditional proportions pis_law
  ## Arguments
  # x: vector of normal variables 
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis_law: the conditional proportions
  
  m <- length(x)
  s <- 1:m # auxiliary variable
  ds=density(s,from = min(s)-5,to=max(s)+5,n=1000)
  #bandwidth
  h=ds$bw
  pval <- 2*pnorm(-abs(x))
  p.est <-rep(0, m)
  for (i in 1:m) { 
    kht<-dnorm(s-i, 0, h)
    kht[i] = 0 # considers only j!=i
    p.est[i]<-sum(kht[which(pval>=tau)])/((1-tau)*sum(kht))
  }
  p.est[which(p.est>1)] <-1
  return(1-p.est)
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


awor_sideinfo.func <- function(x, pis, mu0, mu1, sd0, sd1, q=0.05){
  ## oracle LATLA weight function for asymmetric setting with gamma 
  # controlling the level of asymmetry.
  # x: vector of normal variables
  # pis: conditional probabilities
  # gamma: level of asymmetry
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
  num <- (1-pis)*dnorm(x,mu0,sd0)
  den <- pis*dnorm(x,mu1,sd1)
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
      den <- pis[i]*dnorm(pt,mu1[i],sd1[i])
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
      den <- pis[i]*dnorm(nt,mu1[i],sd1[i])
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

awor_sideinfo_test.func <- function(x, pis, mu0, mu1, sd0, sd1, q=0.05){
  ## oracle LATLA weight function for asymmetric setting with gamma 
  # controlling the level of asymmetry.
  # x: vector of normal variables
  # pis: conditional probabilities
  # gamma: level of asymmetry
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
  num <- (1-pis)*dnorm(x,mu0,sd0)
  den <- pis*(0.5*dnorm(x,mu1,sd1) + 0.5*dnorm(x,-mu1,sd1))
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
      den <- pis[i]*(0.5*dnorm(pt,mu1[i],sd1[i])+0.5*dnorm(pt,-mu1[i],sd1[i]))
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
      den <- pis[i]*(0.5*dnorm(nt,mu1[i],sd1[i])+0.5*dnorm(nt,-mu1[i],sd1[i]))
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


# version 1.2.2
library(progress)


#' Automatically choose a threshold with desired confidence level
#'
#' This function first ranks all hypotheses according to the weighted p-values
#' and then determines a threshold along the ranking to control the FDR. See
#' algorithm 1 in reference paper for more details.
#'
#' @param pvs Vector of p-values under the null hypotheses.
#' @param pis The sparsity levels.
#' @param ws Weights for the p-values.
#' @param q FDR level.
#' @return A list of three: 1.number of rejections; 2.rejection threshold;
#'   3.indication of whether each index is rejected.
#' @importFrom stats density dnorm pnorm
#' @export
lasla_thres<-function(pvs, pis, ws, q)
{
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



#' Data-driven weights of LASLA
#'
#' The function calculates the oracle-assisted data-driven weights described
#' in the algorithm 2 in the LASLA paper.
#' @param x Vector of the primary statistics.
#' @param d Distance matrix or multidimensional side information. 
#' @param pis Estimated sparsity level.
#' @param mu0 Vector of means under the null hypotheses.
#'            Unless specified, the null means are assumed to be 0s.
#' @param sd0 Vector of standard deviations under the null hypotheses.
#'            Unless specified, the null standard deviations are assumed to be 1s.
#' @param q FDR level, default value is set to be 0.05.
#' @param h Bandwidth for kernel estimation, default is set to 'auto', where
#'          the bandwidth is automatically selected by build-in method 'ste'.
#' @param eps Determines the size of the neighborhood, only m^(1-eps) data
#'            points will be used for kernel estimation. eps should range
#'            from 0 to 1.
#' @param grid_size Step size of grid-searching for the asymmetric threshold
#'                  on primary statistics.
#' @param offset Extend the searching window by amount specified by offset.
#' @param progress Toggle to turn on/off progress bar
#' @return The oracle-assisted weights.
#' @export
############## Data-driven Weights #####################################
lasla_weights <- function(x, d, pis, mu0, sd0,
                          q=0.05, 
                          h='auto', 
                          eps=0.1, 
                          grid_size=0.1,
                          offset=0,
                          progress=TRUE){
  m <- length(x)
  
  # eps should be in [0,1], check validity
  if ((eps < 0) | (eps > 1)){
    stop('Invalid eps, eps should be in [0,1].')
  }
  nb <- floor(m^(1-eps))
  
  # check default argument of bandwidth
  if (h=='auto'){
    zds=density(x,bw="SJ-ste")
    h=zds$bw
  }
  
  tor.est <- rep(0, m)
  for (i in 1:m) {
    vhi <- dnorm(d[i,], 0, h)
    vhi[i] = 0
    kht <- dnorm(x-x[i], 0, h)
    # estimated conditional density evaluated at ti
    fti <- sum(vhi*kht)/sum(vhi)
    tor.est[i] <- (1-pis[i])*dnorm(x[i],mu0[i],sd0[i])/fti
  }
  st.tor <- sort(tor.est)
  
  # calculating the moving average
  mva <- rep(0,m)
  cumulative_sum <- cumsum(st.tor)
  mva <- cumulative_sum / seq_along(st.tor)
  
  # calculating the oracle threshold
  if (sum(mva<=q)==0){
    th <- -Inf
  }
  else{
    th <- max(which(mva<=q))
    th <- st.tor[th]
  }
  
  # calculating the weights
  weights <- rep(0,m)
  pt <- seq(0, max(x)+offset, grid_size)
  nt <- seq(min(x)-offset, 0, grid_size)
  t_star <- rep(0,m)
  n_p <- length(pt)
  n_n <- length(nt)
  
  if (progress){
    pb <- progress_bar$new(total = m)
  }
  
  for(i in 1:m){
    if (progress){pb$tick()}
    
    vhi <- dnorm(d[i,], 0, h)
    vhi[i]=0
    # only consider the test points in the neighborhood
    if (nb < m){
      st.nb<- sort(vhi)[m-nb]
      vhi[which(vhi<st.nb)]<-0
    }
    # adapt to asymmetry
    if(x[i]>=0){
      pt_mat <- matrix(rep(pt,m),m,n_p,byrow = TRUE)
      x_mat <- matrix(rep(x,n_p),m,n_p)
      kht <- dnorm(x_mat - pt_mat, 0, h)
      ft <- c(vhi %*% kht)/sum(vhi)
      tor<-(1-pis[i])*dnorm(pt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- max(x)
      }
      else{
        t_star[i] = pt[min(which(tor<=th))]
      }
      weights[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      nt_mat <- matrix(rep(nt,m),m,n_n,byrow = TRUE)
      x_mat <- matrix(rep(x,n_n),m,n_n)
      kht <- dnorm(x_mat - nt_mat, 0, h)
      ft <- c(vhi %*% kht)/sum(vhi)
      tor<-(1-pis[i])*dnorm(nt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- min(x)
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



#' Estimated sparsity level
#'
#' The is a kernel-based function that estimates the local sparsity levels
#' which are subsequently used to construct oracle-assisted weights.
#' @param x Vector of the primary statistics.
#' @param d Distance matrix.
#' @param pval p-values of the primary statistics
#' @param tau A threshold to approximate the 'null set', roughly speaking,
#'            we deem p-values greater than tau as null points. Default is
#'            set to 0.9, see reference paper for more discussions on
#'            parameter selection.
#' @param h Bandwidth for kernel estimation, default is set to 'auto',
#'          where the bandwidth is automatically selected by build-in
#'          method 'ste'(solve the equation).
#' @param eps Determines the size of the neighborhood, only m^(1-eps) data
#'            points will be used for kernel estimation. eps should range
#'            from 0 to 1.
#' @return The estimated sparsity levels
#' @export
lasla_pis <- function(x, d, pval, tau=0.9, h='auto', eps=0.1){
  m <- length(x)
  # tau should be in (0,1), check validity
  if ((tau <= 0) | (tau >= 1)){
    stop('Invalid tau, tau should be in (0,1).')
  }
  # eps should be in [0,1], check validity
  if ((eps < 0) | (eps > 1)){
    stop('Invalid eps, eps should be in [0,1].')
  }
  nb <- floor(m^(1-eps))
  if (h=='auto'){
    zds=density(x,bw="SJ-ste")
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


#' Normalize distance matrix
#'
#' The is the recommended rule of thumb to normalize the distance matrix by a
#' data-driven scaling such that the spread of entries in the scaled distance matrix 
#' is similar to that of the primary statistics. We use the interquantile range (IQR)
#' to measure to data spread.
#' @param x Vector of the primary statistics.
#' @param d Distance matrix.
#' @return The estimated sparsity levels
#' @export
#'
normalize_distance <- function(x,d){
  a <- IQR(d)/IQR(x)
  return(d/a)
}




#' Data-driven weights of LASLA when side information is spatial locations
#' This specialized implementation is potentially more efficient than the one for
#' generic side information, when testing size is very large. 
#'
#' The function calculates the oracle-assisted data-driven weights described
#' in the algorithm 2 in the LASLA paper for spatial side information.
#' @param x Vector of the primary statistics.
#' @param s spatial locations of primary statistics.
#' @param pis Estimated sparsity level.
#' @param mu0 Vector of means under the null hypotheses.
#'            Unless specified, the null means are assumed to be 0s.
#' @param sd0 Vector of standard deviations under the null hypotheses.
#'            Unless specified, the null standard deviations are assumed to be 1s.
#' @param q FDR level, default value is set to be 0.05.
#' @param h Bandwidth for kernel estimation, default is set to 'auto', where
#'          the bandwidth is automatically selected by build-in method 'ste'.
#' @param eps Determines the size of the neighborhood, only m^(1-eps) data
#'            points will be used for kernel estimation. eps should range
#'            from 0 to 1.
#' @param grid_size Step size of grid-searching for the asymmetric threshold
#'                  on primary statistics.
#' @param offset Extend the searching window by amount specified by offset.
#' @param progress Toggle to turn on/off progress bar
#' @return The oracle-assisted weights.
#' @export
############## Data-driven Weights #####################################
lasla_spatial_weights <- function(x, s, pis, 
                          scale=NULL,
                          mu0 = NULL, 
                          sd0 = NULL,
                          q=0.05, 
                          h='auto', 
                          eps=0.1, 
                          grid_size=0.1,
                          offset=0,
                          progress=TRUE){
  m <- length(x)
  
  # eps should be in [0,1], check validity
  if ((eps < 0) | (eps > 1)){
    stop('Invalid eps, eps should be in [0,1].')
  }
  nb <- floor(m^(1-eps))
  
  # Check null means and standard deviations
  if (is.null(mu0)){mu0 = rep(0,m)}
  if (is.null(sd0)){sd0 = rep(1,m)}
  
  if (is.null(scale)){
    diff <- s - matrix(rep(s[1,], m), nrow = m, byrow = TRUE)
    d <- sqrt(rowSums(diff ^ 2))
    scale <- IQR(d)/IQR(x)
  }
  
  if (h=='auto'){
    zds=density(x,bw="SJ-ste")
    h=zds$bw
  } else if (h=="nrd0"){
    h=bw.nrd0(x)
  } else if (h=="nrd"){
    h=bw.nrd(x)
  } 
    
  tor.est <- rep(0, m)
  for (i in 1:m) {
    diff <- s - matrix(rep(s[i,], m), nrow = m, byrow = TRUE)
    d <- sqrt(rowSums(diff ^ 2))/scale
    vhi <- dnorm(d, 0, h)
    vhi[i] = 0
    kht <- dnorm(x-x[i], 0, h)
    # estimated conditional density evaluated at ti
    fti <- sum(vhi*kht)/sum(vhi)
    tor.est[i] <- (1-pis[i])*dnorm(x[i],mu0[i],sd0[i])/fti
  }
  st.tor <- sort(tor.est)
  
  # Calculating the moving average
  mva <- rep(0,m)
  cumulative_sum <- cumsum(st.tor)
  mva <- cumulative_sum / seq_along(st.tor)
  
  # calculating the oracle threshold
  if (sum(mva<=q)==0){
    th <- -Inf
  }
  else{
    th <- max(which(mva<=q))
    th <- st.tor[th]
  }
  
  # calculating the weights
  weights <- rep(0,m)
  pt <- seq(0, max(x)+offset, grid_size)
  nt <- seq(min(x)-offset, 0, grid_size)
  t_star <- rep(0,m)
  n_p <- length(pt)
  n_n <- length(nt)
  
  if (progress){
    pb <- progress_bar$new(total = m)
  }
  
  for(i in 1:m){
    if (progress){pb$tick()}
    
    diff <- s - matrix(rep(s[i,], m), nrow = m, byrow = TRUE)
    d <- sqrt(rowSums(diff ^ 2))/scale
    vhi <- dnorm(d, 0, h)
    vhi[i]=0
    # only consider the test points in the neighborhood
    if (nb < m){
      st.nb<- sort(vhi)[m-nb]
      vhi[which(vhi<st.nb)]<-0
    }
    # adapt to asymmetry
    if(x[i]>=0){
      pt_mat <- matrix(rep(pt,m),m,n_p,byrow = TRUE)
      x_mat <- matrix(rep(x,n_p),m,n_p)
      kht <- dnorm(x_mat - pt_mat, 0, h)
      ft <- c(vhi %*% kht)/sum(vhi)
      tor<-(1-pis[i])*dnorm(pt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- max(x)
      }
      else{
        t_star[i] = pt[min(which(tor<=th))]
      }
      weights[i] <- 1-pnorm(t_star[i], mu0[i], sd0[i])
    }
    else{
      nt_mat <- matrix(rep(nt,m),m,n_n,byrow = TRUE)
      x_mat <- matrix(rep(x,n_n),m,n_n)
      kht <- dnorm(x_mat - nt_mat, 0, h)
      ft <- c(vhi %*% kht)/sum(vhi)
      tor<-(1-pis[i])*dnorm(nt,mu0[i],sd0[i])/ft
      tor[which(tor>1)]<-1
      if(length(which(tor<=th))==0){
        t_star[i] <- min(x)
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



#' Estimated sparsity level when side information is spatial locations
#' This specialized implementation is potentially more efficient than the one for
#' generic side information, when testing size is very large. 
#'
#' The is a kernel-based function that estimates the local sparsity levels for
#' spatial side information. 
#' @param x Vector of the primary statistics.
#' @param s spatial locations of primary statistics.
#' @param pval p-values of the primary statistics
#' @param scale Scale for normalizing the distance measure.
#' @param tau A threshold to approximate the 'null set', roughly speaking,
#'            we deem p-values greater than tau as null points. Default is
#'            set to 0.9, see reference paper for more discussions on
#'            parameter selection.
#' @param h Bandwidth for kernel estimation, default is set to 'auto',
#'          where the bandwidth is automatically selected by build-in
#'          method 'ste'(solve the equation).
#' @param eps Determines the size of the neighborhood, only m^(1-eps) data
#'            points will be used for kernel estimation. eps should range
#'            from 0 to 1.
#' @return The estimated sparsity levels
#' @export
lasla_spatial_pis <- function(x, s, pval,
                              scale=NULL,
                              tau=0.9, 
                              h='auto', 
                              eps=0.1,
                              progress=FALSE){
  m <- length(x)
  # tau should be in (0,1), check validity
  if ((tau <= 0) | (tau >= 1)){
    stop('Invalid tau, tau should be in (0,1).')
  }
  # eps should be in [0,1], check validity
  if ((eps < 0) | (eps > 1)){
    stop('Invalid eps, eps should be in [0,1].')
  }
  nb <- floor(m^(1-eps))
  
  
  # Compute the scaling factor used to normalize the distance vector
  if (is.null(scale)){
    diff <- s - matrix(rep(s[1,], m), nrow = m, byrow = TRUE)
    d <- sqrt(rowSums(diff ^ 2))
    scale <- IQR(d)/IQR(x)
  }
  
  if (h=='auto'){
    zds=density(x,bw="SJ-ste")
    h=zds$bw
  } else if (h=="nrd0"){
    h=bw.nrd0(x)
  } else if (h=="nrd"){
    h=bw.nrd(x)
  } 
  
  p.est <- rep(0, m)
  
  if (progress){
    pb <- progress_bar$new(total = m)
  }
  
  for (i in 1:m){
    if (progress){pb$tick()}
    diff <- s - matrix(rep(s[i,], m), nrow = m, byrow = TRUE)
    d <- sqrt(rowSums(diff ^ 2))/scale
    kht <- dnorm(d, 0, h)
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
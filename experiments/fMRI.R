setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(R.matlab)
library(abind)
library(adaptMT)
library(mgcv)

ADHD_train <- readMat("./data/ADHD/data_wm_training_(30,36,30)_interpolate.mat")
ADHD_test <- readMat("./data/ADHD/data_wm_testing_(30,36,30)_interpolate.mat")

x_train <- ADHD_train$X
y_train <- ADHD_train$y
x_test <- ADHD_test$X
y_test <- ADHD_test$y

# Combine train and test
x <- abind(x_train,x_test,rev.along=1)
y <- abind(y_train,y_test,along=1)


# Separate the ADHD samples and the control samples
x1 <- x[,,,which(y==1)]
x0 <- x[,,,which(y==0)]
n1 <- dim(x1)[4]
n0 <- dim(x0)[4]

# Compute the two-sample t-statistics
var1 <- apply(x1, 1:3, var)
mean1 <- apply(x1, 1:3, mean)
var0 <- apply(x0, 1:3, var)
mean0 <- apply(x0, 1:3, mean)

t <- mean1 - mean0
non_zero_indices <- which(t!=0, arr.ind = TRUE)
t[non_zero_indices] <- t[non_zero_indices] / sqrt(var1[non_zero_indices]/n1 + var0[non_zero_indices]/n0)
t.vec <- as.vector(t)
pv.vec<-2*pnorm(-abs(t.vec), 0, 1)

# Randomly adjust those p-values equal to 1
#valid_indices <- which(pv.vec != 1)
p1_indices <- which(pv.vec == 1)
pv.vec[p1_indices] <- runif(length(p1_indices), 0.95, 1- 1e-15)
# Get spatial locations 
spatial_locs <- arrayInd(seq_along(t.vec), dim(t))

# Parameters for ADAPT
spatial_df <- as.data.frame(spatial_locs)
colnames(spatial_df) <- c("d1", "d2", "d3")
pi_formula <- mu_formula <- "s(d1, d2, d3)"

# Identify different regions between the two samples
m <- length(t.vec)
mu0<-rep(0,m)
sd0<-rep(1,m)
fdr_level<-c(0.05)
nrep<-length(fdr_level)


bh.nr<-rep(0,nrep)
lasla.dd.nr<-rep(0,nrep)
adapt.nr<-rep(0,nrep)

bh.de<-matrix(rep(0, m*nrep), nrep, m)
lasla.dd.de<-matrix(rep(0, m*nrep), nrep, m)
adapt.de<-matrix(rep(0, m*nrep), nrep, m)


for (i in 1:nrep) {
  cat("Running with FDR level:", fdr_level[i], "\n")
  q<-fdr_level[i]
  
  bh.res<-bh.func(pv.vec, q)
  bh.nr[i]<-bh.res$nr
  bh.de[i,]<-bh.res$de
  
  pis_lasla<-lasla_spatial_pis(t.vec, spatial_locs, scale=6, pval=pv.vec,
                               tau=bh.func(pv.vec,0.8)$th, h="nrd0",progress = TRUE)

  weight<-lasla_spatial_weights(t.vec, spatial_locs, pis_lasla,
                                scale=6, h="nrd0", progress = TRUE)
  lasla.dd.res <- lasla_thres(pv.vec, pis_lasla, weight, q)
  lasla.dd.nr[i]<-lasla.dd.res$nr
  lasla.dd.de[i,]<-lasla.dd.res$de
  
  # library("mgcv")
  
  # adapt.res <- adapt_gam(x=spatial_df, pvals=pv.vec,
  #                        pi_formula = pi_formula,
  #                        mu_formula = mu_formula,
  #                        alphas = c(q),
  #                        nfits=5,
  #                        niter_fit=5,
  #                        verbose = list(print = TRUE,
  #                                       fit = TRUE,
  #                                       ms = TRUE))
  # library("splines")
  # x <- spatial_df[1:subsize,]
  # formulas <- paste0("ns(x, df = ", 6:10, ")")
  # dist <- beta_family()
  # adapt.res <- adapt_glm(x=x, pvals=pv.vec[1:subsize],
  #                        pi_formulas = formulas,
  #                        mu_formulas = formulas,
  #                        dist = dist, nfits = 10)
  # 
  # adapt.res <-  adapt_xgboost(x=spatial_locs ,pvals=pv.vec,
  #                             verbose = list(print = TRUE,
  #                                            fit = TRUE,
  #                                            ms = FALSE),
  #                             piargs = list("nrounds" = 50,
  #                                           "max_depth" = 1,
  #                                           "nthread" = 1,
  #                                           "verbose" = 0),
  #                             muargs = list("nrounds" = 50,
  #                                           "max_depth" = 1,
  #                                           "nthread" = 1,
  #                                           "verbose" = 0),
  #                             alphas = c(q),
  #                             nfits = 10)
  # 
  # adapt.nr[i] <- length(which(adapt.res$qvals <= q))
  # adapt.de[i, which(adapt.res$qvals <= q)] <- 1
}


#######################################
#             Save Results            #
#######################################
data_dir <- "./results"

rejections <- data.frame()
nr = cbind(bh.nr,lasla.dd.nr)
de <- list(bh.de, lasla.dd.de)

method_names <- c("BH", "LASLA")

rejections <- data.frame(Method=method_names,
                         FDR=rep(fdr_level, length(method_names)),
                         nr=nr[1,],
                         de=I(de))

save(rejections, file=sprintf("%s/fMRI.RData", data_dir))



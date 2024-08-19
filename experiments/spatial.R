setwd("C:/Users/liang/Documents/GitHub/r-lasla")
source('./methods/lasla_funcs.R')
source('./methods/utils.R')

library(R.matlab)
library(abind)
library(adaptMT)
library(mgcv)

# Generate a 2-dim x
n <- 900
# x1 <- x2 <- seq(-100, 100, length.out = 30)
x1 <- x2 <- seq(1, 30, length.out = 30)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")

# Generate p-values (one-sided z test)
# Set all hypotheses in the central circle with radius 30 to be
# non-nulls. For non-nulls, z~N(2,1) and for nulls, z~N(0,1).
#theta <- apply(x, 1, function(coord){sum(coord^2) < 1200})
theta <- apply(x, 1, function(coord){coord[1] < 20 & coord[1] >10 & coord[2] <20 & coord[2] > 10})
mu <- ifelse(theta, 2, 0)
set.seed(0)
z <- rnorm(n) + mu
pv <- 1 - pnorm(z)

q <- 0.05

image(matrix(z,30,30), axes = FALSE)


bh.de<-rep(0,n)
lasla.dd.de<-rep(0,n)
adapt.de<-rep(0,n)

formula <- "s(x1, x2)"
# adapt.res <- adapt_gam(x = x, pvals = pv, pi_formulas = formula,
#                  mu_formulas = formula, nfits = 5)
adapt.res <-  adapt_xgboost(as.matrix(x) ,pv,
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
adapt.nr <- length(which(adapt.res$qvals <= q))
adapt.de[which(adapt.res$qvals <= q)] <- 1
print(adapt.nr)
image(matrix(adapt.de, 30, 30), axes=FALSE)
adapt.fdp<-sum((1-theta)*adapt.de)/max(sum(adapt.de), 1)
adapt.ntp<-sum(theta*adapt.de)/sum(theta)

bh.res<-bh.func(pv, q)
bh.de<-bh.res$de
bh.nr<-bh.res$nr
print(bh.nr)
image(matrix(bh.de, 30, 30), axes=FALSE)
bh.fdp<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
bh.ntp<-sum(theta*bh.de)/sum(theta)

x_mat <- as.matrix(x)
pis_lasla<-lasla_spatial_pis(z, x_mat, pval=pv, tau=bh.func(pv,0.8)$th)

weight<-lasla_spatial_weights(z,x_mat,pis_lasla,progress = TRUE)
lasla.dd.res<-lasla_thres(pvs=pv, pis_lasla, ws=weight, q)
lasla.dd.de<-lasla.dd.res$de
image(matrix(lasla.dd.de, 30, 30), axes=FALSE)
lasla.dd.fdp<-sum((1-theta)*lasla.dd.de)/max(sum(lasla.dd.de), 1)
lasla.dd.ntp<-sum(theta*lasla.dd.de)/sum(theta)


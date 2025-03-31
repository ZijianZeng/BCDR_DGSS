library(Rcpp)
library(RcppArmadillo)
library(RcppDist)

####################################
# Data load & basic setting
####################################

# load data 
# data$X is the n-by-q, covariates
# data$Y is the n-by-p, nodes
data = readRDS( '../../data/data.rds'  )

# random seed
set.seed(12345)

# preprocess and load sampler
source('init.R')
sourceCpp('sampler.cpp')

# run sampler
Cpp_sampler( Y_data,  X_data, maxiter,
             Rcpp_beta, Rcpp_b, Rcpp_tau, 
             Rcpp_tilde_tau, sigma2.curr, 
             s2.curr, Rcpp_gammaIJK, Rcpp_gammaIJ,
             xpi.curr,  nodepi.curr,
             t.curr, a.xpi, b.xpi, a.nodepi,  b.nodepi,
             a.sigma,  b.sigma ,
             d)

source('save_bs.R')

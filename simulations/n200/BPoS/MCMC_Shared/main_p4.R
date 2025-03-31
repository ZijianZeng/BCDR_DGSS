library(Rcpp)
library(RcppArmadillo)
library(RcppDist)

####################################
# Data load & basic setting
####################################

cpath = getwd()

for( sim in 31:40){
  if( !file.exists( paste0('./res') )) {
    dir.create( paste0('./res' ))
  }
  
  if( !file.exists( paste0('./res/',sim)) ){
    dir.create( paste0('./res/',sim) )
  }
  
  setwd( paste0('./res/',sim))
  # load data 
  # data$X is the n-by-q, covariates
  # data$Y is the n-by-p, nodes
  data = readRDS( paste0(cpath,'/../../data/',sim,'.rds'  )) 
  # preprocess and load sampler
  source( paste0(cpath,'/init.R'))
  sourceCpp( paste0(cpath,'/sampler.cpp'))
  # run sampler
  Cpp_sampler( Y_data,  X_data, maxiter,
               Rcpp_beta, Rcpp_b, Rcpp_tau, 
               Rcpp_tilde_tau, sigma2.curr, 
               s2.curr, Rcpp_gammaIJK, Rcpp_gammaIJ,
               xpi.curr,  nodepi.curr,
               t.curr, a.xpi, b.xpi, a.nodepi,  b.nodepi,
               a.sigma,  b.sigma ,
               d)
  source( paste0(cpath,'/save_bs.R'))
  
  setwd( cpath)
}








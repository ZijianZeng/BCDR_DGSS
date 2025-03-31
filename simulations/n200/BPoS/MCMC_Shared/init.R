####################################
# Data read 
####################################
Y_data = data$Y
X_data = data$X

# load dimension
N = dim(Y_data)[1]
p = dim(Y_data)[2]
q = dim(X_data)[2]

####################################
# set num of iterations
####################################

# num of iterations
maxiter = 20000
num_update = 100
niter_update = 100

####################################
# prior setting
####################################

a.sigma = 0.1
b.sigma = 0.1
a.nodepi = 1 
b.nodepi = 1
a.xpi = 1
b.xpi = 1
d = 0.05

####################################
# parameter initialization
####################################

# parameters
# b.curr = array( rnorm(p*q*q), c(p,p,q))
# beta.curr = array( rnorm(p*q*q), c( p,p,q))
b.curr = array(0, c(p,p,q))
beta.curr = array(0, c(p,p,q))

tau.curr = array( 1, c(p,p,q))
tilde_tau.curr = array( 1, c(p,p,q))
gammaIJK.curr = array(1, c(p,p,q))
gammaIJ.curr = array( 1, c(p,p,1))

# no diagonal element
for( s in 1:q){
  diag(b.curr[,,s]) = NA
  diag(tau.curr[,,s]) = NA
  diag(tilde_tau.curr[,,s]) = NA
  diag(gammaIJK.curr[,,s]) = NA
  diag(gammaIJ.curr[,,1]) = NA
  diag(beta.curr[,,s]) = NA
}

sigma2.curr =  rep(1, p)
s2.curr = rep(1,q) 
xpi.curr = rep(0.5,q) 
nodepi.curr = rep(0.5,p)
t.curr = 1

####################################
# convert param.init to Cpp input
####################################

Rcpp_beta = matrix(NA, p, (p-1)*q )
Rcpp_b = matrix(NA, p, (p-1)*q )
Rcpp_tau = matrix(NA, p,(p-1)*q )
Rcpp_tilde_tau = matrix(NA, p,(p-1)*q )
Rcpp_gammaIJK = matrix(NA, p,(p-1)*q )
Rcpp_gammaIJ = matrix(NA, p,(p-1) )
for( i in 1:p){
  j_index = 1
  for( j in 1:p){
    if( j != i){
      for( s in 1:q){
        Rcpp_beta[i, (q*(j_index-1)+s)] = beta.curr[i,j,s]
        Rcpp_b[i, (q*(j_index-1)+s)] = b.curr[i,j,s]
        Rcpp_tau[i, (q*(j_index-1)+s)] = tau.curr[i, j, s]
        Rcpp_tilde_tau[i, (q*(j_index-1)+s)] = tilde_tau.curr[ i, j, s]
        Rcpp_gammaIJK[i, (q*(j_index-1)+s)] = gammaIJK.curr[i, j, s]
        Rcpp_gammaIJ[i,j_index] = gammaIJ.curr[i,j,1]
      }
      j_index = j_index + 1
    }
  }
}

# dir check
if( !file.exists( './samples/') ){
  dir.create( './samples/')
  # dir.create('./samples/b/')
  # dir.create('./samples/gammaIJ/')
  # dir.create('./samples/tau/')
  # dir.create('./samples/tilde_tau/')
  # dir.create('./samples/gammaIJK/')
  dir.create('./samples/beta/')
  # dir.create('./samples/sigma2/')
  # dir.create('./samples/s2/')
  # dir.create('./samples/xpi/')
  # dir.create('./samples/nodepi/')
  dir.create('./samples/t/')
}


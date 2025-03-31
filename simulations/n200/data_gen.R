library(huge)
library(pracma)
library(MASS)
library(igraph)

set.seed(123)

# settings
sims = 50
n = 200
d = 25
prob = 0.4
graph = 'random'
nq = 4
q = 10
scale = 2 # scale parameter for diagonal dominance

# function for repeatedly try until valid precision
# thanks Dr. Jingfei Zhang for providing her example codes
# used in High-Dimensional Gaussian Graphical Regression Models with Covariates
# this function is modified based on her PD eigen test part
Y_simulate<-function(n,p,q,tB, Odiag ){
  X = matrix( runif(n*q), n, q)
  X[,1] = 1
  Y<-matrix(0,n,p)
  for(i in 1:n){
    omega<--apply(tB,c(1,2),function(b){b%*%X[i,]})
    diag(omega)<- Odiag
    sigma <- solve(omega)
    if(!all(eigen(omega)$values>0)){ #break if omega(u) is not PD
      break
    }
    if(all(eigen(omega)$values >0)){
      Y[i,] <- mvrnorm(1,rep(0,p),sigma)
    }
  }
  data = list( Y = Y, # observations, 
               X = X, # covariates,
               tB = tB, # coefficient,
               pd_test = apply(abs(Y),1,sum),
               Odiag = Odiag # diagonal
  )
  
  return( data )
}


for( sim in 1:sims){
  L = huge.generator(n = n, d = d, graph = graph, prob = prob,  vis = FALSE)
  edges = as.matrix(L$theta)  
  upper.idx = upper.tri(edges)
  # randomly divided into nq graph
  vec_edges = which( edges[ upper.idx] == 1 )
  coef_idx = split( vec_edges, sample( 1:nq, length(vec_edges), replace = TRUE))
  coef_nidx = lapply( coef_idx, length)
  
  # covariate orders:
  coef_edges = array( 0, c(d,d,q) )
  for( j in 1:nq){
    coef_edges[,,j][upper.idx][ coef_idx[[j]] ] = 1
    coef_edges[,,j] = t( coef_edges[,,j] ) + coef_edges[,,j]
  }
  
  # generate covariates
  X = matrix( runif(n*q), n, q) 
  X[,1] = 1
  
  Omega = edges
  for( i in 1:(d-1)){
    for( j in i:d){
      if(Omega[i,j] == 1){
        Omega[i,j] = Omega[j,i] = (-1)^(rbinom(1,1,0.5)) * runif(1,0.35,0.5)
      }
    }
  }
  coef = array(0, c(d,d,q))
  for( j in 1:nq){ 
    coef[,,j][ upper.idx ][ which( coef_edges[,,j][ upper.idx ] == 1 )] = Omega[ upper.idx ][ which( coef_edges[,,j][ upper.idx] == 1 )]
    
    coef[,,j] = coef[,,j] + t(coef[,,j])
  }
  
  # rescale to ensure diagonal dominance
  tB_temp = array(0, c(d,d,q))
  tB = array(0,c(d,d,q))
  
  for(j in 1:d){
    if( (sum(abs(coef[,j,]))/scale) != 0 ){
      tB_temp[,j,] = coef[,j,] / (sum(abs(coef[,j,]))/scale) #ensure diagonal dominance
    }
  }
  
  for(j in 1:(q)){
    tB[,,j] = (tB_temp[,,j]+t(tB_temp[,,j]))/2
  }
  
  # check eigens
  eigens = matrix(0, n, 1)
  for( i in 1:n){
    tem_Omega = matrix(0, d,d)
    for( j in 1:q){
      tem_Omega = tB[,,j] * X[i,j] + tem_Omega
    }
    eigens[i] = min( eig(tem_Omega) )
  }
  ## generate data ##
  Odiag = 1.2 # this value should be greater than min(eigens) above 
              # it will provide higher chance to generate a valid precision
              # because of diagonal dominance
  pd_test=0
  h = 1
  while( (pd_test==0) ){
    data<-Y_simulate(n,d,q,tB, Odiag)
    pd_test<-min(data[[4]])
    print( paste0( 'the ',h,' try in sim ', sim,'. hdp_test: ',pd_test)) #if this value is greater than 0, it means all covariances are PD
    h = h+1
    
    if( h > 100){ 
      h = 1
      # if too many tries, re-sample beta 
      # it may be the case randomly sampled beta is not good for precision
      Omega = edges
      for( i in 1:(d-1)){
        for( j in i:d){
          if(Omega[i,j] == 1){
            Omega[i,j] = Omega[j,i] = (-1)^(rbinom(1,1,0.5)) * runif(1,0.35,0.5)
          }
        }
      }
      coef = array(0, c(d,d,q))
      for( j in 1:nq){ 
        coef[,,j][ upper.idx ][ which( coef_edges[,,j][ upper.idx ] == 1 )] = Omega[ upper.idx ][ which( coef_edges[,,j][ upper.idx] == 1 )]
        
        coef[,,j] = coef[,,j] + t(coef[,,j])
      }
      
      # rescale to ensure diagonal dominance
      tB_temp = array(0, c(d,d,q))
      tB = array(0,c(d,d,q))
      
      for(j in 1:d){
        if( (sum(abs(coef[,j,]))/scale) != 0 ){
          tB_temp[,j,] = coef[,j,] / (sum(abs(coef[,j,]))/scale) #ensure diagonal dominance
        }
      }
      
      for(j in 1:(q)){
        tB[,,j] = (tB_temp[,,j]+t(tB_temp[,,j]))/2
      }
    }
  }
  
  if( !file.exists( './data') ){
    dir.create( paste0('./data/' ))
  }
  
  saveRDS( data, paste0( './data/',sim,'.rds'))  
  write.table(data$Y, paste0( './data/Z_',sim,'.csv'), row.names = FALSE, col.names = FALSE, sep = ',')
  write.table(data$X[,2:q], paste0( './data/U_',sim,'.csv'), row.names = FALSE, col.names = FALSE, sep = ',')
  
  for( tem_j in 1:nq){
    write.table( data$tB[,,tem_j], paste0( './data/b',tem_j,'_',sim,'.csv'),
                 row.names = FALSE, col.names = FALSE, sep = ',')
  }
  
  rm( data )
}



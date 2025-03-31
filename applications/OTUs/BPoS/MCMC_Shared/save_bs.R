# save beta samples
burnins = 10000

tem = as.matrix( read.table( paste0( './samples/beta/',1,'.data')))
beta = array(NA, c( (maxiter-burnins), dim(tem)) )

for( i in 1:(maxiter-burnins)){
  beta[i,,] = as.matrix( read.table( paste0( './samples/beta/', (burnins+i) ,'.data')))     
}
saveRDS( beta, paste0( './beta_samples.rds'))

tem_res = apply(beta, c(2,3), median)
coef_edges = array(NA, c(p,p,q))
for( tem_i in 1:p){
  j_index = 1
  for( tem_j in 1:p){
    if( tem_j != tem_i){
      for( tem_s in 1:q){
        coef_edges[tem_i,tem_j,tem_s] = ((tem_res[tem_i, (q*(j_index-1)+tem_s)])!=0)*1
      }
      j_index = j_index + 1
    }
  }
}

saveRDS( coef_edges, paste0( './raw_hat_edges.rds'))



library(readr)
library(progress)

method = c('glmnet',
           'GMMReg','BSGSSS', 
           'MCMC_Shared')

n = 200
nsim = 50
table = array(NA, c(length(method), 4, nsim) )

# read regression results
h = 1
for( m in  c('glmnet','GMMReg', 'BSGSSS', 
             'MCMC_Shared')){
  for(tem_d in 1:nsim){
    data = readRDS( paste0( './data/',tem_d,'.rds') )
    p = dim(data$Y)[2]
    q = dim(data$X)[2]
    
    if( sum( m == c("glmnet"))){
      tem_res = readRDS( paste0('./glmnet/res/',tem_d,'.rds'))
      coef_edges = tem_res$edges
    }
    
    if( sum( m == c("GMMReg"))){
      coef_edges = array(NA, c(p,p,q))
      for( tem_j in 1:q){
        coef_edges[,,tem_j] = (as.matrix( read_csv( paste0('./GMMReg/','b',tem_j,'_',tem_d,'.csv'), 
                                                    col_types = cols(),col_names = FALSE))!=0)*1
        diag( coef_edges[,,tem_j]) = NA
      }
    }
    
    if( sum( m == "BSGSSS")){
      tem_res = readRDS( paste0( './BSGSSS/res/res_',tem_d,'.rds'))
      coef_edges = tem_res$edges
    }
    
    if( sum( m == c('MCMC_Shared'))){
      tem_res = readRDS(  paste0( './BPoS/',m,'/res','/',tem_d,'/raw_hat_edges.rds') )
      coef_edges = tem_res
    }
    
    tg = array(0, c(p,p,q))
    pg = array(0, c(p,p,q))
    for( g in 1:q){
      tg[,,g] = (data$tB[,, g]!=0)*1
      pg[,,g] = (coef_edges[,,g]!=0)*1
      
      tg[,,g] = ((tg[,,g] + t(tg[,,g]))!=0) * 1
      pg[,,g] = ((pg[,,g] + t(pg[,,g]))!=0) * 1
      diag(tg[,,g]) = NA
      diag(pg[,,g]) = NA
    }
    
    tg = (apply(tg,3,sum, na.rm = T) != 0) * 1
    pg = (apply(pg,3,sum, na.rm = T) != 0) * 1
    
    TP = sum( (pg!=0)*(tg!=0), na.rm = T) 
    TN = sum( (pg==0)*(tg==0), na.rm = T) 
    FP = sum( (pg!=0)*(tg==0), na.rm = T) 
    FN = sum( (pg==0)*(tg!=0), na.rm = T) 
    
    table[h, 1, tem_d ] = TP / (TP + FN ) 
    table[h, 2, tem_d] = FP / (FP + TN )
    table[h, 3, tem_d] =  TP / ( TP + 1/2*(FP + FN) )
    table[h, 4, tem_d] = (TP * TN - FP * FN ) / ( sqrt(TP+FP) * sqrt(TP + FN) * sqrt( TN + FP ) * sqrt( TN + FN ))
    
  }
  
  
  h = h + 1
  print( paste0('Results of ',m,' loaded;'))
}

# saveRDS(table, 'quick_table_covariates_n200.rds')
# table = readRDS( 'quick_table_covariates_n200.rds')
dim(table)

metric = c( 'tpr', 'fpr', 'F1', 'MCC')

tem_table = matrix( NA, 2*length(metric), length(method))

for( i in 1:length(method)){
  for( j in 1:length(metric)){
    tem_table[ 2*(j-1)+1, i ] = sum( table[i,j,], na.rm = T)  / ( sum( !is.na(table[i,j,])) )
    tem_table[ 2*j, i ] = sd( table[i,j,], na.rm = T) / sqrt(( sum( !is.na(table[i,j,])) ))
  }
}

rownames( tem_table ) = c( 'tpr','sd', 'fpr','sd', 'F1','sd', 'MCC','sd')
colnames( tem_table ) = method

round( tem_table, 3)
library(xtable)
xtable( tem_table, digits = 3)

# glmnet and BSGSSS are expected to have FPR = 1 and fail the MCC = NA
# due to almost always include all covariates;
# meanwhile due to the small sample, there are about 2 out of 100 cases they
# didn't include all covariates, 2 case of the n200 (50 cases), 0 case of the 
# n500 (the other 50 cases)
table
which( !is.na(table[1,4,]))
which( !is.na(table[3,4,]))

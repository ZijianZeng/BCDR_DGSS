library(readr)
library(progress)

method = c('glmnet',
           'GMMReg','BSGSSS', 
           'MCMC_Shared')

n_list = c(200, 500)


for( n in n_list){
  if( n == 200){
    tem_d = 13
  }else if(n == 500){
    tem_d = 10
  }
  
  if( !file.exists( paste0('./n',n,'/') ) ){
    dir.create( paste0('./n',n,'/' ))
  }
  
  data = readRDS( paste0( '../n',n,'/data/',tem_d,'.rds') )
  p = dim(data$Y)[2]
  q = dim(data$X)[2]
  
  tg = array(0, c(p,p,q))
  for( g in 1:q){
    tg[,,g] = (data$tB[,, g]!=0)*1
    tg[,,g] = ((tg[,,g] + t(tg[,,g]))!=0) * 1
    diag(tg[,,g]) = 0
  }
  
  for( g in 1:q){
    png(file=paste0("./n",n,"/sim_S",tem_d,"_B",g,".png",sep=""))
    heatmap( ( tg[,,g] != 0 )*1, col = c('grey30', 'grey85'), symm = TRUE, Rowv=NA, Colv=NA,
             revC=TRUE, labRow = FALSE, labCol = FALSE )
    dev.off()
  }
  
  # read regression results
  for( m in  c('glmnet','GMMReg', 'BSGSSS', 
               'MCMC_Shared')){
    
    if( sum( m == c("glmnet"))){
      tem_res = readRDS( paste0('../n',n,'/glmnet/res/',tem_d,'.rds'))
      coef_edges = tem_res$edges
    }
    
    if( sum( m == c("GMMReg"))){
      coef_edges = array(NA, c(p,p,q))
      for( tem_j in 1:q){
        coef_edges[,,tem_j] = (as.matrix( read_csv( paste0('../n',n,'/GMMReg/','b',tem_j,'_',tem_d,'.csv'), 
                                                    col_types = cols(),col_names = FALSE))!=0)*1
        diag( coef_edges[,,tem_j]) = NA
      }
    }
    
    if( sum( m == "BSGSSS")){
      tem_res = readRDS( paste0( '../n',n,'/BSGSSS/res/res_',tem_d,'.rds'))
      coef_edges = tem_res$edges
    }
    
    if( sum( m == c('MCMC_Shared'))){
      tem_res = readRDS(  paste0( '../n',n,'/BPoS/',m,'/res','/',tem_d,'/raw_hat_edges.rds') )
      coef_edges = tem_res
    }
    
    tg = array(0, c(p,p,q))
    pg = array(0, c(p,p,q))
    for( g in 1:q){
      pg[,,g] = (coef_edges[,,g]!=0)*1
      pg[,,g] = ((pg[,,g] + t(pg[,,g]))!=0) * 1
      diag(pg[,,g]) = 0
      
      png(file=paste0("./n",n,"/model_",m,"_B",g,".png",sep=""))
      heatmap( ( pg[,,g] != 0 )*1, col = c('grey30', 'grey85'), symm = TRUE, Rowv=NA, Colv=NA,
               revC=TRUE, labRow = FALSE, labCol = FALSE )
      dev.off()
    }
    
    print( paste0('Results of ',m,' loaded;'))
  }
  print( paste0('Sample size ',n,' Done;'))
    
}


  


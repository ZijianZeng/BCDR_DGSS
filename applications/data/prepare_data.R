# =================================
#           load data
# =================================
# load selected data
OTU = as.matrix(read.csv("OTU.csv", row.names = 1, header= TRUE, check.names=FALSE))
cyto = as.matrix(read.csv("cyto.csv", row.names = 1, header = TRUE, check.names=FALSE))

# check subject ID
OTU_id = sapply( rownames(OTU), function(x){
  strsplit(x, '_')[[1]][1]
})
cyto_id = sapply( rownames(cyto), function(x){
  strsplit(x, '_')[[1]][1]
})
all.equal( as.vector(OTU_id), as.vector(cyto_id) )

rownames(OTU) = OTU_id
rownames(cyto) = cyto_id

# =================================
#     process OTU: normalization
# =================================
library(compositions)
n = dim(OTU)[1]
p = dim(OTU)[2]

Y = clr(OTU)
tem_Y = matrix( Y, n, p)

# zero sample mean
tem_Y = tem_Y - matrix( rep(   apply(tem_Y, 2, mean) ,n), n, byrow = T)

colnames(tem_Y) = colnames(Y)
rownames(tem_Y) = rownames(Y)


# ======================================
#     process cyto: log scale and [0,1]
# ======================================
tem_X = cyto

for( i in 1:ncol(cyto)){
  # process
  tem = cyto[,i]
  tem[which(tem<1)] = 0
  
  # log scale
  tem = log( tem + 1)
  # normalization [0,1]
  a = min(tem)
  b = max(tem)
  
  tem_X[,i] = (tem - a)/(b-a)
}
tem_X = cbind( rep(1, n), tem_X )
q = ncol(tem_X)

data = list( Y = tem_Y,
             X = tem_X)
# saveRDS( data,  paste0( './data/data.rds' ) )
# write.table(data$Y, paste0( './data/Z.csv'), row.names = FALSE, col.names = FALSE, sep = ',')
# write.table(data$X[,2:q], paste0( './data/U.csv'), row.names = FALSE, col.names = FALSE, sep = ',')


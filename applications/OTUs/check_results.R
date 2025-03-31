library(readr)
library(progress)
library(igraph)
library(xtable)

method = c(
  'GMMReg', 'MCMC_Shared')

data = readRDS('./data/data.rds')
p = dim(data$Y)[2]
q = dim(data$X)[2]
sedge = array( NA, c( length(method), p*(p-1)/2, q))

MCSH_network = array( NA, c(p,p,q))
GMMR_network = array( NA, c(p,p,q))

# read regression results'
h = 1
for( m in  c(
  'GMMReg', 'MCMC_Shared')){
  
  if( sum( m == c('GMMReg')) ){
    coef_edges = array(NA, c(p,p,q))
    for( tem_j in 1:q){
      coef_edges[,,tem_j] = (as.matrix( read_csv( paste0('./GMMReg/results/b',tem_j,'.csv'),
                                                  col_types = cols(), col_names = FALSE))!=0)*1
      diag( coef_edges[,,tem_j]) = NA
    }
  }
  
  if( sum( m == c('MCMC_Shared'))){
    tem_res = readRDS(  paste0( './BPoS/',m,'/raw_hat_edges.rds') )
    coef_edges = tem_res
  }
  
  
  pg = array(0, c(p,p,q))
  for( j in 1:q){
    pg[,,j] = (coef_edges[,,j]!=0)*1
    pg[,,j] = ((pg[,,j] + t(pg[,,j]))!=0) * 1
    diag(pg[,,j]) = NA
    
    sedge[h,,j] = pg[,,j][upper.tri(diag(p))]
    
    if ( m == 'MCMC_Shared'){
      MCSH_network[,,j] = pg[,,j]
    }else if ( m == 'GMMReg'){
      GMMR_network[,,j] = pg[,,j]
    }
  }
  h = h + 1
  print( paste0('Results of ',m,' loaded;'))
}

##################################################

# Table 2 

##################################################

# number of edges selected for covariates cross model given Cyto
detect_cov = matrix(NA, length(method), q)
rownames(detect_cov) = method
colnames(detect_cov) = colnames(data$X)
for( j in 1:q){
  detect_cov[,j] = apply(sedge[,,j],1,sum)  
}
# Table 2
detect_cov


# read Nathan's results
nathan_res = readRDS('./data/nathan_res_ref.rds')
otu_order = c(as.matrix(read.table("./data/otu_order.txt",header = F)))

# graph selected by the methods
GMCSH_network = (apply( MCSH_network, c(1,2), sum )!=0)*1
GMCSH_network = GMCSH_network[otu_order,otu_order]

##################################################
# check graphs
##################################################
# preliminary
library(reshape2)
library(ggplot2)
load("./data/momspi16S_tax.rda")
OTU_names = as.character(c(as.matrix(read.table("./data/otu_names.csv", header = F))))
tax = momspi16S_tax[OTU_names,]
phylum = tax[,2]
unique_phylums = unique(phylum)
order = c()
col = c()
group_size = c()
for(i in 1:length(unique_phylums)) {
  order = c(order,which(phylum == unique_phylums[i]))
  col = c(col,rep(i,length(which(phylum == unique_phylums[i]))))
  group_size = c(group_size,length(which(phylum == unique_phylums[i])))
}
cols = col
cols_network = matrix(1,nrow = 90,ncol = 90)
for(i in 1:90){
  cols_network[which(cols == i),which(cols == i)] = i + 1 
}
colors = c('white','grey','blue','red','orange','yellow','green','grey90', '#ccccff','#ffb2b2','#ffcc99','#ffff99','#ccff99')
col_sel = c(0,1,2,3,4,5,6,100,200,300,400,500,600)


################
# MCMC_Shared
################
# overall network
tem_network = (apply( MCSH_network, c(1,2), sum) != 0)*1
diag(tem_network) = 0
tem_network = tem_network[otu_order, otu_order]
tem_network[upper.tri(tem_network)] = tem_network[upper.tri(tem_network)]*100

ppi.hpHC.hm <- melt((cols_network * as.matrix(tem_network)))
ppi.hpHC.hm[,2]  = sort(rep(1:90,90))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")


if( !file.exists( paste0('./Figures') )) {
  dir.create( paste0('./Figures' ))
}

png( file = paste0('./Figures/MCMC_Shared_network.png'),
     width = 1440, height = 1440, res = 120 )
p1 = ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_text(size = 16,face = 'bold'),axis.title.y=element_text(size = 16,face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  
  ## base of boxes
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[2], y = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[3], y = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[4], y = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[5], y = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
  geom_segment(aes(x = 89.5, y = 89.5, xend = 90.5, yend = 89.5)) +
  ## right side of boxes
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[1], x = c(54.5,68.5,77.5,84.5,89.5)[1], yend = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[2], x = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[3], x = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[4], x = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[5], x = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[5], xend = c(54.5,68.5,77.5,84.5,89.5)[5]),color = 'black') +
  geom_segment(aes(y = 0.5, x = 0.5, yend = 90.5, xend = 90.5)) +
  ##
  geom_text(y = 27.5, x = 26, label = "F i r m i c u t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 61.5, x = 60, label = "A c t i n o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 73, x = 71.5, label = "B a c t e r o i d e t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 81.5, x = 80, label = "P r o t e o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 87.5, x = 86, label = "F u s o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 90.0, x = 89, label = "T M 7",angle = 0,hjust = 1,size = 6,vjust = 0.0) +
  xlab("Microbe") + ylab("Microbe")
print(p1)
dev.off()

################
# plot the shared edges from Nathan
################
# store the shared network
GMCSH_shared = matrix(0, p, p)
GMCSH_shared[upper.tri(diag(p))][which( GMCSH_network[upper.tri(diag(p))] * nathan_res$edge[upper.tri(diag(p))] == 1)] = 1
GMCSH_shared = GMCSH_shared + t(GMCSH_shared)

tem_network = GMCSH_shared
tem_network[upper.tri(tem_network)] = tem_network[upper.tri(tem_network)]*100

ppi.hpHC.hm <- melt((cols_network * as.matrix(tem_network)))
ppi.hpHC.hm[,2]  = sort(rep(1:90,90))
colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")


png( file = paste0('./Figures/MCMC_shared_shared.png'),
     width = 1440, height = 1440, res = 120 )
p1 = ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) + 
  geom_tile(color = 'grey') + 
  scale_fill_manual(values = colors) + 
  # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
  # geom_hline(yintercept = v.lines+.38, colour = 'red') + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_text(size = 16,face = 'bold'),axis.title.y=element_text(size = 16,face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  
  ## base of boxes
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[2], y = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[3], y = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[4], y = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
  geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[5], y = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
  geom_segment(aes(x = 89.5, y = 89.5, xend = 90.5, yend = 89.5)) +
  ## right side of boxes
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[1], x = c(54.5,68.5,77.5,84.5,89.5)[1], yend = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[2], x = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[3], x = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[4], x = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
  geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[5], x = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[5], xend = c(54.5,68.5,77.5,84.5,89.5)[5]),color = 'black') +
  geom_segment(aes(y = 0.5, x = 0.5, yend = 90.5, xend = 90.5)) +
  ##
  geom_text(y = 27.5, x = 26, label = "F i r m i c u t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 61.5, x = 60, label = "A c t i n o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 73, x = 71.5, label = "B a c t e r o i d e t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 81.5, x = 80, label = "P r o t e o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 87.5, x = 86, label = "F u s o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
  geom_text(y = 90.0, x = 89, label = "T M 7",angle = 0,hjust = 1,size = 6,vjust = 0.0) +
  xlab("Microbe") + ylab("Microbe")
print(p1)
dev.off()

# covariates network
for( j in 1:30){
  tem_network = MCSH_network[,,j]
  diag(tem_network) = 0
  tem_network = tem_network[otu_order, otu_order]
  tem_network[upper.tri(tem_network)] = tem_network[upper.tri(tem_network)]*100
  
  ppi.hpHC.hm <- melt((cols_network * as.matrix(tem_network)))
  ppi.hpHC.hm[,2]  = sort(rep(1:90,90))
  colnames(ppi.hpHC.hm) <- c("Column","Row","PPI")
  
  if( j == 1){
    subname = 'base'
  }else{
    subname = colnames(data$X)[j]
  }
  
  png( file = paste0('./Figures/MCMC_Shared_', subname ,'.png'),
       width = 1440, height = 1440, res = 120 )
  p1= ggplot(data = ppi.hpHC.hm, aes(x=Column, y=Row, fill=factor(PPI))) +
    geom_tile(color = 'grey') +
    scale_fill_manual(values = colors[col_sel %in% unique(ppi.hpHC.hm$PPI)]) +
    # geom_vline(xintercept = cumsum(group_size), colour = 'red') +
    # geom_hline(yintercept = v.lines+.38, colour = 'red') +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
          axis.title.x=element_text(size = 16,face = 'bold'),axis.title.y=element_text(size = 16,face = 'bold'),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    
    ## base of boxes
    geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[2], y = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
    geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[3], y = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
    geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[4], y = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
    geom_segment(aes(x = c(0.5,54.5,68.5,77.5,84.5)[5], y = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
    geom_segment(aes(x = 89.5, y = 89.5, xend = 90.5, yend = 89.5)) +
    ## right side of boxes
    geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[1], x = c(54.5,68.5,77.5,84.5,89.5)[1], yend = c(54.5,68.5,77.5,84.5,89.5)[1], xend = c(54.5,68.5,77.5,84.5,89.5)[1]),color = 'black') +
    geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[2], x = c(54.5,68.5,77.5,84.5,89.5)[2], yend = c(54.5,68.5,77.5,84.5,89.5)[2], xend = c(54.5,68.5,77.5,84.5,89.5)[2]),color = 'black') +
    geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[3], x = c(54.5,68.5,77.5,84.5,89.5)[3], yend = c(54.5,68.5,77.5,84.5,89.5)[3], xend = c(54.5,68.5,77.5,84.5,89.5)[3]),color = 'black') +
    geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[4], x = c(54.5,68.5,77.5,84.5,89.5)[4], yend = c(54.5,68.5,77.5,84.5,89.5)[4], xend = c(54.5,68.5,77.5,84.5,89.5)[4]),color = 'black') +
    geom_segment(aes(y = c(0.5,54.5,68.5,77.5,84.5)[5], x = c(54.5,68.5,77.5,84.5,89.5)[5], yend = c(54.5,68.5,77.5,84.5,89.5)[5], xend = c(54.5,68.5,77.5,84.5,89.5)[5]),color = 'black') +
    geom_segment(aes(y = 0.5, x = 0.5, yend = 90.5, xend = 90.5)) +
    ##
    geom_text(y = 27.5, x = 26, label = "F i r m i c u t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
    geom_text(y = 61.5, x = 60, label = "A c t i n o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
    geom_text(y = 73, x = 71.5, label = "B a c t e r o i d e t e s",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
    geom_text(y = 81.5, x = 80, label = "P r o t e o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
    geom_text(y = 87.5, x = 86, label = "F u s o b a c t e r i a",angle = 0,hjust = 1,size = 6,vjust = 0.5) +
    geom_text(y = 90.0, x = 89, label = "T M 7",angle = 0,hjust = 1,size = 6,vjust = 0.0) +
    xlab("Microbe") + ylab("Microbe")
  print(p1)
  dev.off()
}


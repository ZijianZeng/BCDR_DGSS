library(readr)
library(progress)
library(igraph)
library(xtable)
library(ggplot2)
library(reshape2)

data = readRDS('./data/data.rds')
p = dim(data$Y)[2]
q = dim(data$X)[2]
MCSH_edges = matrix( NA, p*(p-1)/2, q)
MCSH_network = array( NA, c(p,p,q))

# read regression results'
tem_res = readRDS(  paste0( './BPoS/MCMC_Shared/raw_hat_edges.rds') )
coef_edges = tem_res

pg = array(0, c(p,p,q))
for( j in 1:q){
  pg[,,j] = (coef_edges[,,j]!=0)*1
  pg[,,j] = ((pg[,,j] + t(pg[,,j]))!=0) * 1
  diag(pg[,,j]) = NA
  
  MCSH_edges[,j] = pg[,,j][upper.tri(diag(p))]
    
  MCSH_network[,,j] = pg[,,j]
}


##################################################

# Covariate_Edge_Table

##################################################

Covariate_Edge_Table = matrix(NA, q, q)

for( i in 1:q){
  for( j in 1:q){
    if( i != j ){
      edge1 = MCSH_edges[,i]
      edge2 = MCSH_edges[,j]
      Covariate_Edge_Table[i,j] = sum((edge1+edge2) == 2)
    }
  }
}
rownames(Covariate_Edge_Table) = colnames(data$X)
colnames(Covariate_Edge_Table) = colnames(data$X)

df <- as.data.frame(as.table(Covariate_Edge_Table[2:q,2:q]))
colnames(df) <- c("Row", "Column", "Value")
df$Row <- factor(df$Row, levels = rev(levels(df$Row)))

# Plot the heatmap
ggplot(df, aes(x = Column, y = Row, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = Value), color = "white", size = 2) +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  ) +
  xlab("")+
  ylab("")+
  ggtitle("Number of Edges Impacted by at least the 2 common Cytokines") 


##################################################

# Edge_Covariate_Table

##################################################

# we have in total 90 nodes, corresponding to 4,005 edges
Edge_mapping_table = matrix(NA, p*(p-1)/2, 3+(q-1))
h = 1
for( i in 1:p){
  for(j in i:p){
    if( j != i){
      Edge_mapping_table[h,1] = h
      Edge_mapping_table[h,2] = as.numeric(colnames(data$Y)[i])
      Edge_mapping_table[h,3] = as.numeric(colnames(data$Y)[j])
      Edge_mapping_table[h,4:(3+(q-1))] = as.numeric(MCSH_edges[h,2:q])
      h = h + 1
    }
  }
}
colnames(Edge_mapping_table) = c('edge_no','node_name', 'node_name', colnames(data$X)[2:q])
write.csv(Edge_mapping_table,'edges_cytokines_table.csv')

edge_covariate_table = MCSH_edges
rownames(edge_covariate_table) = 1:4005
colnames(edge_covariate_table) = colnames(data$X)

df = as.data.frame(as.table(edge_covariate_table[,2:q]))
colnames(df) <- c("Row", "Column", "Value")
df$Row = factor(df$Row, levels = rev(unique(df$Row)))

# Plot the heatmap
ggplot(df, aes(x = Column, y = Row, fill = factor(Value))) +
  geom_tile() +
  scale_fill_manual(
    values = c("lightgray", "black"),  # Adjusted colors
    labels = c("0", "1")
  ) +
  scale_y_discrete(
    breaks = levels(df$Row)[seq(1, length(levels(df$Row)), by = 100)]  
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),  # Rotate column labels
    axis.text.y = element_text(size = 6),                        # Reduce y-axis label size
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "",
    y = "",
    fill = "",
    title = "Edges and Impacting Cytokines"
  )




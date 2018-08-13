library(mclust)
library(dplyr)
library(magrittr)


comparaison_k_means <- function(data, metadata, nb_cluster=2, nb_graph=1, nb_start=50, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_ilr <- data %>% MAP() %>% ilr(base_binaire = base_binaire)
  
  k_data_ilr <- kmeans(data_ilr, nb_cluster, nstart=nb_start)
  k_data <- kmeans(data, nb_cluster, nstart=nb_start)
  
  table1 <- table(metadata, k_data$cluster)
  table2 <- table(metadata, k_data_ilr$cluster)
  
  
  grob <- NULL
  if(nb_graph>0){
    b_coord <- data %>% MAP() %>% biplot()
    grob <- c(graph_biplot_normale(b_coord, k_data$cluster, title="comptage", nb_graph = nb_graph, coord_biplot = TRUE), graph_biplot_normale(b_coord, metadata, title = "correct", nb_graph = nb_graph, coord_biplot = TRUE),  graph_biplot_normale(b_coord, k_data_ilr$cluster, title = "ilr", nb_graph = nb_graph, coord_biplot = TRUE))
  }
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
}


comparaison_hclust <- function(data, metadata, nb_cluster, nb_graph, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_ilr <- data %>% MAP() %>% ilr(base_binaire = base_binaire)
  
  hclust_data <- data %>% dist() %>% as.dist() %>% hclust(method="ward.D") %>% cutree(k=nb_cluster)
  hclust_data_ilr <- data_ilr %>% dist() %>% as.dist() %>% hclust(method="ward.D") %>% cutree(k=nb_cluster)
  
  table1 <- table(metadata, hclust_data)
  table2 <- table(metadata, hclust_data_ilr)
  
  grob <- NULL
  if(nb_graph>0){
    b_coord <- data %>% MAP() %>% biplot()
    grob <- c(graph_biplot_normale(b_coord, hclust_data, title="comptage", nb_graph = nb_graph, coord_biplot = TRUE), graph_biplot_normale(b_coord, metadata, title = "correct", nb_graph = nb_graph, coord_biplot = TRUE),  graph_biplot_normale(b_coord, hclust_data_ilr, title = "ilr", nb_graph = nb_graph, coord_biplot = TRUE))
  }
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
  
}


comparaison_Mclust <- function(data, metadata, nb_cluster, nb_graph, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_ilr <- data %>% MAP() %>% ilr(base_binaire = base_binaire)
  
  Mclust_data <- data %>% Mclust(G=nb_cluster)
  Mclust_data_ilr <- data_ilr %>% Mclust(G=nb_cluster)
  
  table1 <- table(metadata, Mclust_data$classification)
  table2 <- table(metadata, Mclust_data_ilr$classification)
  
  grob <- NULL
  if(nb_graph>0){
    b_coord <- data %>% MAP() %>% biplot()
    grob <- c(graph_biplot_normale(b_coord, Mclust_data$classification, title="comptage", nb_graph = nb_graph, coord_biplot = TRUE), graph_biplot_normale(b_coord, metadata, title = "correct", nb_graph = nb_graph, coord_biplot = TRUE),  graph_biplot_normale(b_coord, Mclust_data_ilr$classification, title = "ilr", nb_graph = nb_graph, coord_biplot =TRUE))
  }
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
  
}




adjusted_rand_index <- function(table1){
  
  part1 <- table1 %>% rowSums() %>% choose(2) %>% sum()
  part2 <- table1 %>% colSums() %>% choose(2) %>% sum()
  part3 <- table1 %>% sum() %>% choose(2)
  
  ARI <- table1 %>% choose(2) %>% sum() - part1*part2/part3
  ARI_denom <- 0.5*( part1+part2 ) - part1*part2/part3
  
  ARI/ARI_denom
}


mutual_information_distance <- function(table1){
  
  P <- rowSums(table1)/sum(table1)
  P <- P+(P==0)*1
  HY <- -P%*%log(P)
  
  P <- colSums(table1)/sum(table1)
  P <- P+(P==0)*1
  HC <- -P%*%log(P)
  
  PC <- t(t(table1)/colSums(table1))
  PC <- PC+(PC==0)*1
  HYC <- sum(colSums((-PC)*log(PC)))
  I <- HY-HYC
  
  HC+HY-2*I
}
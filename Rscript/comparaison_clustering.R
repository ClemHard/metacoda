library(mclust)
library(dplyr)
library(magrittr)



#' perform two k-means clustering, one on a data set of count data, 
#' the second on the same data set but after transforming it to compositionnal 
#' and apply ilr transformation. 
#' 
#'  
#' @param data a dataset (count data)
#' @param metadata a numeric vector containing the name of the group (the group of the samples of the dataset data) 
#' @param nb_cluster the number of clusters use in the kmeans
#' @param nb_graph the number of graph return (show the graphically the results of clusterings and the real group)
#' @param nb_start the number of random sets chosen in kmeans algorithm
#' @param base_binaire binary sequential matrix use for the ilr transformation (default the basis sequence)
#'
#' @return two table contingence (one for the clustering on the count dataset, the other on ilr dataset) 
#' @return list of graphics (showing the results of the clustering on the first axes of an acp)
#' @author Clement Hardy
#' @export

comparaison_k_means <- function(data, metadata, nb_cluster=2, nb_graph=1, nb_start=50, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_ilr <- data %>% MAP() %>% ilr(base_binaire = base_binaire)
  
  k_data_ilr <- kmeans(data_ilr, nb_cluster, nstart=nb_start)
  k_data <- kmeans(data, nb_cluster, nstart=nb_start)
  
  table1 <- table(k_data$cluster, metadata) %>% change_name_table()
  table2 <- table(k_data_ilr$cluster, metadata) %>% change_name_table()
  
  
  grob <- NULL
  if(nb_graph>0){
    b_coord <- data %>% MAP() %>% biplot()
    grob <- c(graph_biplot_normale(b_coord, k_data$cluster, title="comptage", nb_graph = nb_graph, coord_biplot = TRUE), graph_biplot_normale(b_coord, metadata, title = "correct", nb_graph = nb_graph, coord_biplot = TRUE),  graph_biplot_normale(b_coord, k_data_ilr$cluster, title = "ilr", nb_graph = nb_graph, coord_biplot = TRUE))
  }
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
}


#' perform two hierarchical clustering, one on a data set of count data, 
#' the second on the same data set but after transforming it to compositionnal 
#' and apply ilr transformation. 
#' 
#'  
#' @param data a dataset (count data)
#' @param metadata a numeric vector containing the name of the group (the group of the samples of the dataset data) 
#' @param nb_cluster the number of clusters use in the clustering
#' @param nb_graph the number of graph return (show the graphically the results of clusterings and the real group)
#' @param base_binaire binary sequential matrix use for the ilr transformation (default the basis sequence)
#'
#' @return two table contingence (one for the clustering on the count dataset, the other on ilr dataset) 
#' @return list of graphics (showing the results of the clustering on the first axes of an acp)
#' @author Clement Hardy
#' @export

comparaison_hclust <- function(data, metadata, nb_cluster, nb_graph, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_ilr <- data %>% MAP() %>% ilr(base_binaire = base_binaire)
  
  hclust_data <- data %>% dist() %>% as.dist() %>% hclust(method="ward.D") %>% cutree(k=nb_cluster)
  hclust_data_ilr <- data_ilr %>% dist() %>% as.dist() %>% hclust(method="ward.D") %>% cutree(k=nb_cluster)
  
  table1 <- table(hclust_data, metadata)  %>% change_name_table()
  table2 <- table(hclust_data_ilr, metadata)  %>% change_name_table()
  
  grob <- NULL
  if(nb_graph>0){
    b_coord <- data %>% MAP() %>% biplot()
    grob <- c(graph_biplot_normale(b_coord, hclust_data, title="comptage", nb_graph = nb_graph, coord_biplot = TRUE), graph_biplot_normale(b_coord, metadata, title = "correct", nb_graph = nb_graph, coord_biplot = TRUE),  graph_biplot_normale(b_coord, hclust_data_ilr, title = "ilr", nb_graph = nb_graph, coord_biplot = TRUE))
  }
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
  
}


#' perform two gaussian mixture clustering, one on a data set of count data, 
#' the second on the same data set but after transforming it to compositionnal 
#' and apply ilr transformation. 
#'  
#' @param data a dataset (count data)
#' @param metadata a numeric vector containing the name of the group (the group of the samples of the dataset data) 
#' @param nb_cluster the number of gaussian use in the gaussian mixture
#' @param nb_graph the number of graph return (show the graphically the results of clusterings and the real group)
#' @param base_binaire binary sequential matrix use for the ilr transformation (default the basis sequence)
#'
#' @return two table contingence (one for the clustering on the count dataset, the other on ilr dataset) 
#' @return list of graphics (showing the results of the clustering on the first axes of an acp)
#' @author Clement Hardy
#' @export
#' @import mclust

comparaison_Mclust <- function(data, metadata, nb_cluster, nb_graph, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_ilr <- data %>% MAP() %>% ilr(base_binaire = base_binaire)
  
  Mclust_data <- data %>% mclust::Mclust(G=nb_cluster)
  Mclust_data_ilr <- data_ilr %>% mclust::Mclust(G=nb_cluster)
  
  table1 <- table(Mclust_data$classification, metadata) %>% change_name_table()
  table2 <- table(Mclust_data_ilr$classification, metadata) %>% change_name_table()
  
  grob <- NULL
  if(nb_graph>0){
    b_coord <- data %>% MAP() %>% biplot()
    grob <- c(graph_biplot_normale(b_coord, Mclust_data$classification, title="comptage", nb_graph = nb_graph, coord_biplot = TRUE), graph_biplot_normale(b_coord, metadata, title = "correct", nb_graph = nb_graph, coord_biplot = TRUE),  graph_biplot_normale(b_coord, Mclust_data_ilr$classification, title = "ilr", nb_graph = nb_graph, coord_biplot =TRUE))
  }
  
  list( comptage_table=table1, ilr_table=table2, graphics=grob)
  
}

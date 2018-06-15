library(mclust)
library(MASS)

simu_melange_gaussien <- function(n, probability, mean, Sigma){
  
  NB <- rmultinom(1, n, probability)
  
  sample1 <- NULL
  
  for(i in 1:length(probability)){
    
    if(NB[i]>0){
      sample1 <- rbind(sample1, mvrnorm(NB[i], mean[[i]], Sigma[[i]]))
    }
    
  }
  
  list(data= sample1, metadata= factor(rep(1:length(NB), NB)))
}



bootstrap_ilr <- function(data, nb_cluster, nb_sample=nrow(data)){
  
  Mclust_data <- Mclust(data, G=nb_cluster)
  
  probability <- Mclust_data$parameters$pro
  mean_data <- list()
  Sigma_data <- list()
  
  for(i in 1:length(probability)){
    mean_data[[i]] <- Mclust_data$parameters$mean[,i]
    Sigma_data[[i]] <- Mclust_data$parameters$variance$sigma[,,i]
  }

  new_sample <- simu_melange_gaussien(nb_sample, probability, mean_data, Sigma_data)
  
  new_sample
}



bootstrap_pca <- function(data, nb_cluster, nb_sample=nrow(data), nb_axe){
  
  data_biplot <- biplot(data)
  
  Mclust_data <- Mclust(data_biplot$coord[,1:nb_axe], G=nb_cluster)
  probability <- Mclust_data$parameters$pro
  
  mean_data <- list()
  Sigma_data <- list()
  
  for(i in 1:length(probability)){
    mean_data[[i]] <- Mclust_data$parameters$mean[,i]
    Sigma_data[[i]] <- Mclust_data$parameters$variance$sigma[,,i]
  }
  
  new_sample <- simu_melange_gaussien(nb_sample, probability, mean_data, Sigma_data)
  
  D <- ncol(new_sample$data)
  n <- nrow(new_sample$data)
  
  Sigma <- sum(data_biplot$values[(nb_axe+1):length(data_biplot$values)])/length(data_biplot$vector[1,]) 
  Z <- (new_sample$data + mvrnorm(n, rep(0, D), diag(D)*Sigma) ) %*% t(data_biplot$vector[,1:nb_axe])
  
  new_sample$data <- Z %>% ilr_inverse() %>% perturbation(center_data(data))
  
  new_sample
  
}


bootstrap_comptage <- function(data, nb_cluster, nb_sample=nrow(data), PCA=FALSE){
  
  if(!PCA){
    new_sample <- bootstrap_ilr(data %>% MAP() %>% ilr(), nb_cluster, nb_sample)
    data_ilr <- ilr_inverse(new_sample$data)
  }else{
    new_sample <- bootstrap_pca(data %>% MAP(), nb_cluster, nb_sample, nb_axe=15) 
    data_ilr <- new_sample$data
  }
  
  sample_comptage <- matrix(0, nrow=nrow(data_ilr), ncol=ncol(data_ilr))
  
  result_sample <- apply(data_ilr, 1, function(x){
      rmultinom(1, (data %>% sum())/nrow(data) %>% round(), x)
  }) %>%t()
  
  list(data=result_sample, metadata=new_sample$metadata)
}

library(mclust)
library(MASS)
library(capushe)

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



bootstrap <- function(data, nb_axe=NULL, nb_cluster=NULL, nb_sample=nrow(data), type="comptage"){
  
  apprent <- apprentissage(data, nb_axe = nb_axe, nb_cluster = nb_cluster)
  if(type=="comptage"){
    new.sample <- simulation(apprent, nb_sample = nb_sample)
  } else if(type=="MAP"){
    new.sample <- simulation_MAP(apprent, nb_sample = nb_sample)
  }
  new.sample <- list(data=new.sample, d=apprent$d, nb_cluster=apprent$nb_cluster, apprent=apprent)
  new.sample
}



nb_cluster_capushe <- function(data, max_nb_cluster=round(nrow(data)/3)){
  
  nb_cluster <- 1:max_nb_cluster
  
  erreur <- sapply(nb_cluster, function(x){
    Mclust(data, x)$bic 
  })
  
  which.max(erreur)
  
}



nb_axe_capushe <- function(data, data_biplot=FALSE){
  
  new_data <- data
  if(!data_biplot){
    new_data <- new_data %>% biplot()
  }
  D <- ncol(new_data$coord)
  RMSE <- sapply(1:D, function(x){mean((new_data$coord[,x]^2))}) %>% sort() %>% cumsum() %>% sort(decreasing = TRUE) 
  RMSE <- RMSE/sum(RMSE)
  z <- data.frame(1:D,1:D,1:D,RMSE)
  
  i <- nrow(z)
  while(class(try(capushe(z[1:i,]), silent=TRUE))=="try-error"){
    i <- round(i*0.9)
    if(i<=10){
      stop("impossible de fitter une régression")
    }
  }
  as.numeric(capushe(z[1:i,])@DDSE@model)
}



apprentissage_pca <- function(data, nb_cluster=NULL, nb_axe=NULL){
  
  data_biplot <- biplot(data)
  
  if(is.null(nb_axe)){
    nb_axe <- nb_axe_capushe(data_biplot, data_biplot = TRUE)  
  }
  
  if(is.null(nb_cluster)){
    nb_cluster <- nb_cluster_capushe(data_biplot$coord[,1:nb_axe])
  }
  
  Mclust_data <- Mclust(data_biplot$coord[,1:nb_axe], G=nb_cluster)
  probability <- Mclust_data$parameters$pro
  
  mean_data <- list()
  Sigma_data <- list()
  
  for(i in 1:length(probability)){
    mean_data[[i]] <- Mclust_data$parameters$mean[,i]
    Sigma_data[[i]] <- Mclust_data$parameters$variance$sigma[,,i]
  }
  
  Sigma <- sum(data_biplot$values[(nb_axe+1):length(data_biplot$values)])/(length(data_biplot$vector[,1])-nb_axe)
  W <- data_biplot$vector[,1:nb_axe]
  
  result <- list(data=data, W=W, pi_k=probability, mean=mean_data, Sigma=Sigma_data, noise=Sigma, d=nb_axe, nb_cluster=nb_cluster)
  class(result) <- "apprentissage"
  
  result
}


apprentissage <- function(data, nb_axe=NULL, nb_cluster=NULL){
  
  new.data <- data %>% MAP()
  result <- apprentissage_pca(new.data, nb_axe = nb_axe, nb_cluster = nb_cluster)
  result$data <- data
  
  result
}



simulation_ilr <- function(result, nb_sample=nrow(result$data)){
  new_sample <- simu_melange_gaussien(nb_sample, result$pi_k, result$mean, result$Sigma)
  
  Z <- t((result$W) %*% t(new_sample$data))
  
  D <- ncol(Z)
  n <- nrow(Z)
  Z <- Z + mvrnorm(n, rep(0, D), result$noise*diag(D))
  
}



simulation_MAP <- function(result, nb_sample=nrow(result$data)){
  
  Z <- simulation_ilr(result, nb_sample)
  
  new_sample <- Z %>% ilr_inverse() %>% perturbation(center_data(result$data %>% MAP())) 
  new_sample
}



simulation <- function(result, nb_sample=nrow(result$data), type="comptage"){
  
  if(type=="ilr"){
    simu_ilr <- simulation_ilr(result, nb_sample = nb_sample)
    rownames(simu_ilr) <- paste("sample", 1:nrow(simu_ilr))
    colnames(simu_ilr) <- paste("coord", 1:ncol(simu_ilr))
    return(simu_ilr)
  }
  
  simu_MAP <- simulation_MAP(result, nb_sample)
  
  if(type=="MAP"){
    rownames(simu_MAP) <- paste("sample", 1:nrow(simu_MAP))
    colnames(simu_MAP) <- colnames(result$data)
    return(simu_MAP)
  }
  
  deep <- apply(result$data, 1, sum)
  
  result_sample <- apply(simu_MAP, 1, function(x){
    deep_simu <- sample(deep, 1)
    (round(x*deep_simu-1)>0)*rmultinom(1, round(deep_simu), x)
    #(round(x*deep_simu-1)>0)*(round(x*deep_simu-1))
  }) %>%t()
  
  rownames(result_sample) <- paste("sample", 1:nrow(simu_MAP))
  colnames(result_sample) <- colnames(result$data)
  
  result_sample
}

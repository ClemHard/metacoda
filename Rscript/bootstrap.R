library(mclust)
library(MASS)
library(capushe)
library(dplyr)

source("Rscript/coda.R")


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


#' find the most represented value per OTU 
#'
#' @param data a dataset (OTU in column, sample in line)
#'
#' @return a matrix with two column, in first column the most represented value,
#'  in the second the proportion of this value
#'
#' @author Clement Hardy
#' @export


max_abundance_value_OTU <- function(data){
  
  mat <- matrix(0, ncol=2, nrow=ncol(data))
  uni <- apply(data, 2, unique)
  
  if(class(uni)=="matrix") return(NULL)
  
  for(i in 1: length(uni)){
    
    temp <-  sapply(uni[[i]], function(x){
      mean(x==data[,i])})
    
    mat[i, 1] <- max(temp)
    mat[i, 2] <- uni[[i]][which.max(temp)]
  }
  
  colnames(mat) <- c("zero_inflated", "value")
  as.data.frame(mat)
}

#' Zero inflated is the method to force zero in a simulated dataset, this function
#' find the most common value per OTU, it's proporion per OTU, it's proportion per groups
#' and finally compute the coefficient for a bernoulli.
#' 
#' @param data a dataset (OTU in column, sample in line)
#'
#' @return a tibble with 6 column, first column the OTU,second: group of the sample, 
#'                                 third: the most common value, four: the proportion of the value per group, five: the number of sample in the group 
#'                                 six: the proportion of the value per OTU
#' 
#' @author Clement Hardy
#' @export dplyr


zero_inflated <- function(data, classification){
  
  data <- signif(data, digits = 10)
  zero_inflated <- max_abundance_value_OTU(data)
  iid <- which(zero_inflated$zero_inflated>=0)
  
  iid2 <- data.frame(iid=(iid %>% rep(nrow(data))), value=(zero_inflated$value[iid] %>% rep(nrow(data))))
  
  iid2 <- iid2[order(iid2$iid), ]
  iid_cluster <- NULL
  
  nb_sample_cluster <- sapply(1:length(unique(classification)), function(x){sum(x==classification)})
  
  
  if(!(iid %>% is.integer0())){
    pseudo_comptage <- 0.1
    iid_cluster <- data.frame(cluster=classification, data=as.vector(data[, iid]), OTU=iid2$iid, value=iid2$value ,row.names = NULL) %>%
      group_by(OTU, cluster, value) %>% summarise(zero=sum(data==value))
    
    iid_cluster$nb_sample_cluster <- nb_sample_cluster
    iid_cluster$zero <- iid_cluster$zero/ (iid_cluster$nb_sample_cluster + pseudo_comptage)
    iid_cluster$zero_inflated_coeff <- rep(zero_inflated$zero_inflated[iid], rep(length(unique(classification)),length(iid) ))
    
  }
  
  iid_cluster
}


#' produces one or more sample from a gaussian mixture 
#'
#' @param n the number of sample to produces
#' @param probability a numeric vector containing the proportion of each gaussian
#' @param mean a list containing the mean of each gaussian
#' @param Sigma a list containing the covariance matrix of each gaussian
#'
#' @return a list containing a matrix of the produce sample and 
#'         a character vector specify the gaussian from which the sample was produce
#'
#' @author Clement Hardy
#' @export
#' @import MASS
#' @import stats


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


#' produces new sample from a dataset of OTU table,
#' the method is a base on compositional data with isometric log ratio transformation 
#' a dimension reduction and a gaussian mixture is used to produce the news samples
#' the dimension of the latent space could be given in argument or by default find with the heuristic slope.
#' the number of gaussian in the gaussian mixture could be given in argument or default find with the bic criteria
#'
#' @param data the dataset to simulate
#' @param nb_axe the dimension of the latent space
#' @param nb_cluster the number of gaussian in the gaussian mixture
#' @param nb_sample the number of sample to produce
#' @param  base_binaire binary sequential matrix use for the ilr transformation
#'
#' @return a list containing a matrix with the produce sample, an object of class "apprentissage", the dimension of the latent space, the number of gaussian
#' 
#' @author Clement Hardy
#' @export
#' 


bootstrap <- function(data, nb_axe=NULL, nb_cluster=NULL, nb_sample=nrow(data), type="comptage", base_binaire=Base_binary_matrix(ncol(data))){
  
  apprent <- apprentissage(data, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)
  if(type=="comptage"){
    new.sample <- simulation(apprent, nb_sample = nb_sample)
  } else if(type=="MAP"){
    new.sample <- simulation_MAP(apprent, nb_sample = nb_sample)$data
  } else if(type=="ilr"){
    new.sample <- simulation_ilr(apprent, nb_sample = nb_sample)$data
  }
  new.sample <- list(data=new.sample, d=apprent$d, nb_cluster=apprent$nb_cluster, apprent=apprent)
  new.sample
}


#' find the optimal number of cluster to produce a dataset with the criteria BIC
#' 
#' @param data the dataset the gaussian need to learn
#' @param max_nb_cluster the maximun number of gaussian tolerate
#'
#' @return the optimal number of gaussian
#'
#' @author Clement Hardy
#' @export
#' @import mclust

nb_cluster <- function(data, max_nb_cluster=round(nrow(data)/3)){
  
  nb_cluster <- 1:max_nb_cluster
  
  erreur <- sapply(nb_cluster, function(x){
    Mclust(data, x)$bic 
  })
  
  which.max(erreur)
  
}


#' find the optimal dimension of the latent space with the the heuristic slope
#'
#' @param data the original dataset
#' @param data_biplot boolen , TRUE the dataset is already principal component coordinate, FALSE compositionnal data
#' @param  base_binaire binary sequential matrix use for the ilr transformation
#'
#' @return the optimal dimension of the latent space
#'
#' @author Clement Hardy
#' @export
#' @import capushe

nb_axe_capushe <- function(data, data_biplot=FALSE, base_binaire=Base_binary_matrix(ncol(data))){
  
  new_data <- data
  if(!data_biplot){
    new_data <- new_data %>% biplot(base_binaire = base_binaire)
  }
  D <- ncol(new_data$coord)
  
  RMSE <- sapply(1:D, function(x){mean((new_data$coord[,x]^2))}) %>% sort() %>% cumsum() %>% sort(decreasing = TRUE) 
  RMSE <- RMSE/sum(RMSE)
  
  z <- data.frame(1:D,1:D,1:D,RMSE)
  i <- nrow(z)
  while(class(try(capushe(z[1:i,]), silent=TRUE))=="try-error"){
    i <- round(i*0.9)
    if(i<=5){
      warning("impossible de fitter une régression")
      return(D)
    }
  }
  as.numeric(capushe(z[1:i,])@DDSE@model)
}




#'  fit a gaussian mixture on compositional data after doing a reduction dimension
#'  
#' @param a dataset (compositionnal data)
#' @param nb_axe the dimension of the latent space
#' @param nb_cluster the number of gaussian in the gaussian mixture
#' @param  base_binaire binary sequential matrix use for the ilr transformation (default the basis sequence)
#'
#'
#' @return an object of class "apprentissage"
#' @return data the dataset to simulate
#' @return W the matrix of probalistic principal component
#' @return  probability a numeric vector containing the proportion of each gaussian
#' @return mean a list containing the mean of each gaussian
#' @return Sigma a list containing the covariance matrix of each gaussian
#' 
#' @return noise the variance lost by the probalistic principal component
#' @return d the dimension of the latent space
#' @return nb_cluster the number of gaussian in the gaussian mixture
#' @return  base_binaire binary sequential matrix use for the ilr transformation
#'
#' @author Clement Hardy
#' @export
#' @import mclust
#' 


apprentissage_pca <- function(data, nb_cluster=NULL, nb_axe=NULL, base_binaire=Base_binary_matrix(ncol(data))){
  
  data_biplot <- biplot(data, base_binaire=base_binaire)
  if(is.null(nb_axe)){
    nb_axe <- nb_axe_capushe(data_biplot, data_biplot = TRUE)  
  }
  
  if(is.null(nb_cluster)){
    nb_cluster <- nb_cluster(data_biplot$coord[,1:nb_axe])
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
  if(is.na(Sigma) || Sigma<0) Sigma <- 0
  
  W <- data_biplot$vector[,1:nb_axe]
  
  data_ilr <- data %>% center_scale(scale=FALSE) %>% ilr()
  iid_cluster <- zero_inflated(data_ilr, Mclust_data$classification)
  
  result <- list(data=data, classification_data=Mclust_data$classification, W=W, pi_k=probability, mean=mean_data, Sigma=Sigma_data, noise=Sigma, d=nb_axe, nb_cluster=nb_cluster, base_binaire=base_binaire, zero_inflated_ilr=iid_cluster)
  
  class(result) <- "apprentissage"
  
  result
}



#'  fit a gaussian mixture on a dataset (OTU in column, sample in line)
#'  
#' @param data a dataset (count data)
#' @param nb_axe the dimension of the latent space
#' @param nb_cluster the number of gaussian in the gaussian mixture
#' @param  base_binaire binary sequential matrix use for the ilr transformation (default the basis sequence)
#'
#' @return an object of class "apprentissage"
#' @return data the dataset to simulate
#' @return W the matrix of probalistic principal component
#' @return  probability a numeric vector containing the proportion of each gaussian
#' @return mean a list containing the mean of each gaussian
#' @return Sigma a list containing the covariance matrix of each gaussian
#' @return noise the variance lost by the probalistic principal component
#' @return d the dimension of the latent space
#' @return nb_cluster the number of gaussian in the gaussian mixture
#' @return  base_binaire binary sequential matrix use for the ilr transformation
#' 
#' 
#' @author Clement Hardy
#' @export


apprentissage <- function(data, nb_axe=NULL, nb_cluster=NULL, base_binaire=Base_binary_matrix(ncol(data))){
  
  new.data <- data %>% MAP()
  result <- apprentissage_pca(new.data, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)
  result$data <- data
  
  iid_cluster <- zero_inflated(data, result$classification_data)
  
  result$zero_inflated_comptage <- iid_cluster
  
  result
}


#' produce news samples in isometric coordinate with an object of class "apprentissage"
#'
#' @param result an object of class "apprentissage"
#' @param nb_sample the number of sample to simulate
#'
#' @return the news samples
#' @return the gaussian from which each sample belong
#'
#' @author Clement Hardy
#' @export
#' @import MASS


simulation_ilr <- function(result, nb_sample=nrow(result$data)){
  new_sample <- simu_melange_gaussien(nb_sample, result$pi_k, result$mean, result$Sigma)
  
  Z <- (new_sample$data %*% t(result$W))
  
  D <- ncol(Z)
  n <- nrow(Z)
  Z <- Z + mvrnorm(n, rep(0, D), result$noise*diag(D))
  
  list(data=Z, metadata=new_sample$metadata)
}


#' produce news sample (compositionnal data) using an object of class "apprentissage"
#' 
#' @param result an object of class "apprentissage"
#' @param nb_sample the number of sample to simulate
#'
#' @return the news samples
#' @return the gaussian from which each sample belong
#' 
#' @author Clement Hardy
#' @export



simulation_MAP <- function(result, nb_sample=nrow(result$data)){
  
  Z <- simulation_ilr(result, nb_sample)
  new_sample <- Z$data %>% ilr_inverse(base_binaire = result$base_binaire) %>% perturbation(center_data(result$data %>% MAP())) 
  list(data=new_sample, metadata=Z$metadata)
}


#' produce news sample (compositionnal data) using an object of class "apprentissage"
#' after simulated compositionnal data, the method use a zero inflated to force zero.
#' the proportion of zero in news samples is determine by the zero inflated function
#' 
#'
#' @param result an object of class "apprentissage"
#' @param nb_sample the number of sample to simulate
#' @param type the type of samples to simulate (ilr coordinate, compositionnal data, count data)
#'
#' @return the news samples
#'
#' @author Clement Hardy
#' @export
#' @import stats


simulation <- function(result, nb_sample=nrow(result$data), type="comptage"){
  
  if(type=="ilr"){
    simu_ilr <- simulation_ilr(result, nb_sample = nb_sample)
    rownames(simu_ilr$data) <- paste("sample", 1:nrow(simu_ilr$data))
    colnames(simu_ilr$data) <- paste("coord", 1:ncol(simu_ilr$data))
    return(simu_ilr$data)
  }
  
  simu_MAP <- simulation_MAP(result, nb_sample)
  
  if(type=="MAP"){
    rownames(simu_MAP$data) <- paste("sample", 1:nrow(simu_MAP$data))
    colnames(simu_MAP$data) <- colnames(result$data)
    return(simu_MAP$data)
  }
  
  
  deep <- apply(result$data, 1, sum)
  zero_inflated_comptage <- result$zero_inflated_comptage
  
  result_sample <- sapply(1:nrow(simu_MAP$data), function(x){
    
    deep_simu <- sample(deep, 1)
    temp <- zero_inflated_comptage[which(zero_inflated_comptage$cluster==simu_MAP$metadata[x]),]
    coeff <- (1-temp$zero) %*% simu_MAP$data[x, ]
    
    deep_simu <- deep_simu/coeff
    
    if(any(x>0)){
      r <- rmultinom(1, round(deep_simu), simu_MAP$data[x, ])
    }
    else rep(0, length(simu_MAP[x, ]))
  }) %>% t()
  
  
  
  if(!is.null(result$zero_inflated_comptage)){
    for(i in 1:nrow(zero_inflated_comptage)){
      if(zero_inflated_comptage$value[i]==0) {
        e <- rbinom(length(which(simu_MAP$metadata==zero_inflated_comptage$cluster[i])), 1, 1-zero_inflated_comptage$zero[i])
        result_sample[which(simu_MAP$metadata==zero_inflated_comptage$cluster[i]), zero_inflated_comptage$OTU[i]] <- result_sample[which(simu_MAP$metadata==zero_inflated_comptage$cluster[i]), zero_inflated_comptage$OTU[i]] * e
      }
    }
  }
  
  
  rownames(result_sample) <- paste("sample", 1:nrow(simu_MAP$data))
  colnames(result_sample) <- colnames(result$data)
  
  result_sample
}

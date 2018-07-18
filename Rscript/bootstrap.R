library(mclust)
library(MASS)
library(capushe)

is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


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


zero_inflated <- function(data, classification){
  
  data <- signif(data, digits = 10)
  zero_inflated <- max_abundance_value_OTU(data)
  iid <- which(zero_inflated$zero_inflated>0)

  iid2 <- data.frame(iid=(iid %>% rep(nrow(data))), value=(zero_inflated$value[iid] %>% rep(nrow(data))))
 
  iid2 <- iid2[order(iid2$iid), ]
  iid_cluster <- NULL
  
  nb_sample_cluster <- sapply(1:length(unique(classification)), function(x){sum(x==classification)})

  
  if(!(iid %>% is.integer0())){
    pseudo_comptage <- 1
    iid_cluster <- data.frame(cluster=classification, data=as.vector(data[, iid]), OTU=iid2$iid, value=iid2$value) %>%
      group_by(OTU, cluster, value) %>% summarise(zero=sum(data==value))
    
    iid_cluster$nb_sample_cluster <- nb_sample_cluster
    iid_cluster$zero <- iid_cluster$zero / (iid_cluster$nb_sample_cluster + pseudo_comptage)
    iid_cluster$zero_inflated_coeff <- rep(zero_inflated$zero_inflated[iid], rep(length(unique(classification)),length(iid) ))

    iid_cluster <- iid_cluster %>% filter(zero>0)
  }
  iid_cluster
}



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



nb_cluster_capushe <- function(data, max_nb_cluster=round(nrow(data)/3)){
  
  nb_cluster <- 1:max_nb_cluster
  
  erreur <- sapply(nb_cluster, function(x){
    Mclust(data, x)$bic 
  })
  
  which.max(erreur)
  
}



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
    if(i<=10){
      stop("impossible de fitter une régression")
    }
  }
  as.numeric(capushe(z[1:i,])@DDSE@model)
}



apprentissage_pca <- function(data, nb_cluster=NULL, nb_axe=NULL, base_binaire=Base_binary_matrix(ncol(data))){

  data_biplot <- biplot(data, base_binaire=base_binaire)
  
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
  
  data_ilr <- data %>% center_scale(scale=FALSE) %>% ilr()
  
  iid_cluster <- zero_inflated(data_ilr, Mclust_data$classification)

  result <- list(data=data, classification_data=Mclust_data$classification, W=W, pi_k=probability, mean=mean_data, Sigma=Sigma_data, noise=Sigma, d=nb_axe, nb_cluster=nb_cluster, base_binaire=base_binaire, zero_inflated_ilr=iid_cluster)
  
  class(result) <- "apprentissage"
  
  result
}



apprentissage <- function(data, nb_axe=NULL, nb_cluster=NULL, base_binaire=Base_binary_matrix(ncol(data))){

  new.data <- data %>% MAP()
  result <- apprentissage_pca(new.data, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)
  result$data <- data
  
  iid_cluster <- zero_inflated(data, result$classification_data)
  
  result$zero_inflated_comptage <- iid_cluster
  
  result
}


simulation_ilr <- function(result, nb_sample=nrow(result$data)){
  new_sample <- simu_melange_gaussien(nb_sample, result$pi_k, result$mean, result$Sigma)
  
  Z <- (new_sample$data %*% t(result$W))
  
  D <- ncol(Z)
  n <- nrow(Z)
  Z <- Z + mvrnorm(n, rep(0, D), result$noise*diag(D))

  if(!is.null(result$zero_inflated_ilr)){
    temp <- result$zero_inflated_ilr
    for(i in 1:nrow(temp)){
      e <- rbinom(length(which(result$classification_data==temp$cluster[i])), 1, temp$zero[i])
      Z[which(result$classification_data==temp$cluster[i]), temp$OTU[i]] <- Z[which(result$classification_data==temp$cluster[i]), temp$OTU[i]]*(1-e) + temp$value[i]
    }
   }
  # 
  # for(i in 1:ncol(result$zero_inflated)){
  #   Z[, result$zero_inflated[i, 1]] <- Z[, result$zero_inflated[i, 1]] * rbinom(nrow(Z), 1, 1-result$zero_inflated[i, 2])
  # }

  list(data=Z, metadata=new_sample$metadata)
}




simulation_MAP <- function(result, nb_sample=nrow(result$data)){
  
  Z <- simulation_ilr(result, nb_sample)
  new_sample <- Z$data %>% ilr_inverse(base_binaire = result$base_binaire) %>% perturbation(center_data(result$data %>% MAP())) 
  list(data=new_sample, metadata=Z$metadata)
}



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
  
  result_sample <- apply(simu_MAP$data, 1, function(x){
    deep_simu <- sample(deep, 1)
    zix <- x-1/deep_simu
    #zix[zix <= 1 / deep_simu]  <- 0
    zix[zix <=0]  <- 0
    rmultinom(1, deep_simu, x)
    #(round(x*deep_simu-1)>0)*(round(x*deep_simu-1))
  }) %>%t()
  #result_sample <- (result_sample>1)*(result_sample)
  
  
  if(!is.null(result$zero_inflated_comptage)){
    temp <- result$zero_inflated_comptage
    for(i in 1:nrow(temp)){
      e <- rbinom(length(which(simu_MAP$metadata==temp$cluster[i])), 1, temp$zero[i])
      result_sample[which(simu_MAP$metadata==temp$cluster[i]), temp$OTU[i]] <- result_sample[which(simu_MAP$metadata==temp$cluster[i]), temp$OTU[i]] * (e*temp$value[i] + (1-e))
    }
  }
  
  rownames(result_sample) <- paste("sample", 1:nrow(simu_MAP$data))
  colnames(result_sample) <- colnames(result$data)
  
  result_sample
}

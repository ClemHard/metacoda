library(randomForest)
library(class)
library(e1071)
library(dplyr)
library(doParallel)

source("Rscript/bootstrap.R")
source("Rscript/test_bootstrap.R")


table_to_percentage_table <- function(table){
  nb_row <- nrow(table)
  table <- round(table/matrix(rep(apply(table, 2, sum), nb_row), byrow = TRUE, nrow=nb_row)*100, digits=2)
  table
}

list_table_sum <- function(l){
  
  if(length(l)==0){
    warning("list de longeur 0")
    return(NULL)
  }
  sum_table <- l[[1]]
  for(i in 1:length(l)){
    for(j in 1:length(l[[i]])){
      sum_table[[j]] <- sum_table[[j]]+l[[i]][[j]]
    }
  }
  sum_table
}

apprentissage_supervise <- function(data, metadata){
  
  name_group <- unique(metadata)
  apprent <- lapply(name_group, function(x){
    iid <- which(metadata==x)
    a <- apprentissage(data[iid,])
    a$metadata <- x
    class(a) <- "apprentissage supervise"
    a
  })
  apprent
}


simulation_supervise <- function(apprent, nb_sample=NULL){
  
  data <- NULL
  metadata <- NULL
  if(length(nb_sample)!=length(unique(apprent)) || any(nb_sample<0)){
    nb_sample <- NULL
  }
  for(i in 1:length(apprent)){
    if(is.null(nb_sample)){
      data <- rbind(data, simulation(apprent[[i]]))
      metadata <- c(metadata, rep(i, nrow(apprent[[i]]$data)))
    }else{
      data <- rbind(data, simulation(apprent[[i]], nb_sample = nb_sample[i]))
      metadata <- c(metadata, rep(i, nb_sample[i]))
    }
    
  }
  
  list(data=data, metadata=metadata)
}


bootstrap_supervise <- function(data, metadata, nb_sample=NULL){
  apprent <- apprentissage_supervise(data=data, metadata = metadata)
  simulation_supervise(apprent, nb_sample = nb_sample)
}


test_bootstrap_supervise <- function(data, metadata, nb_train=1, type="comptage", nb_sample=NULL){
  
  apprent <- apprentissage_supervise(data=data, metadata = metadata)
  
  if(type=="MAP") data <- data %>% MAP()
  if(type=="ilr") data <- data %>% MAP() %>% ilr()
  
  train <- create_data_frame_train(data=data, metadata = metadata)
  
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parLapply(cl, 1:nb_train, function(x){
    
    simu <- simulation_supervise(apprent, nb_sample = nb_sample)
    
    if(type=="MAP") simu$data <- simu$data %>% MAP()
    if(type=="ilr") simu$data <- simu$data %>% MAP() %>% ilr()
    
    test <- create_data_frame_test(data=simu$data)
    
    find_group(train, test, simu$metadata)
      
   })
  
  stopCluster(cl)

  sum_table <- lapply(list_table_sum(l), table_to_percentage_table)
  list(all=sum_table, all_train=l)
}


find_group <- function(train, test, metadata_test){
  ### random forest
  forest <- randomForest(metadata~., train)
  pred_forest <- predict(forest, test)
  
  ### kNN
  kNN <- knn(train[, -ncol(train)], test, cl=train[, ncol(train)], k=1)
  
  ### Svm
  Svm <- svm(metadata~., train)
  pred_svm <- predict(Svm, test)
  
  list(random_forest=table(pred_forest, metadata_test), kNN=table(kNN, metadata_test), Svm=table(pred_svm, metadata_test))
}




validation_croise <- function(data, metadata, k=nrow(data)){
  
  name_group <- unique(metadata)
  train <- data.frame(data, as.factor(metadata), row.names = NULL)
  colnames(train) <- c(paste("X",1:ncol(data), sep=""), "metadata")
  
  iid <- sample(rep(1:k,length=nrow(data)))
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parLapply(cl, 1:k, function(x){
    
    pred_kNN <- knn(train=data[x!=iid,], test=data[x==iid,], cl=metadata[x!=iid])
    Svm <- svm(metadata~., train[x!=iid,])
    pred_Svm <-  predict(Svm, train[x==iid,])
    
    table1 <- table(metadata[x==iid], pred_kNN)
    table2 <- table(metadata[x==iid], pred_Svm)
    colnames(table1) <- 1:length(name_group)
    colnames(table2) <- 1:length(name_group)
    
    list(table1, table2)
  })
  
  stopCluster(cl)
  
  all <- list_table_sum(l)
  names(all) <- c("kNN", "Svm")
  all
}


find_real <- function(data_real, data_simu_train, data_simu){
  
  data_train <- rbind(data_real, data_simu_train)
  metadata <- c(rep("real", nrow(data_real)), rep("simu", nrow(data_simu_train))) %>% as.factor()

  train <- create_data_frame_train(data=data_train, metadata = metadata)
  test <- create_data_frame_test(data=data_simu)

  
  forest <- randomForest(metadata~., train)
  pred <- predict(forest, test)

  iid <- which(pred=="real")
  
  iid
}



find_nb_sample_cluster <- function(metadata){
  name_group <- unique(metadata)
  sapply(name_group, function(x){sum(x==metadata)})
}



classificateur_group_real_simu <- function(data, metadata, nb_train=1){
  
  
  boot <- bootstrap(data = data)
  apprent <- apprentissage_supervise(data=data, metadata=metadata)
  
  train <- create_data_frame_train(data = data, metadata=metadata)

  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parLapply(cl, 1:nb_train, function(x){
    
    simu <- simulation_supervise(apprent, nb_sample=10*find_nb_sample_cluster(metadata = metadata))
  
    iid_real <- find_real(data_real = data, data_simu_train = boot$data, data_simu = simu$data)
  
    data_real <- simu$data[iid_real, ]
    metadata_real <- simu$metadata[iid_real]

    test <- create_data_frame_test(data=data_real)
  
  
    find_group(train=train, test = test, metadata_test =  metadata_real)
  })
  
  stopCluster(cl)
  
  sum_table <- list_table_sum(l)

  
  sum_table <- lapply(sum_table, table_to_percentage_table)
  list(all=sum_table, all_train=l)
}


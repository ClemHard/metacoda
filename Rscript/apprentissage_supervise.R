library(randomForest)
library(class)
library(e1071)
library(dplyr)
source("Rscript/bootstrap.R")


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
  
  train <- data.frame(data, as.factor(metadata), row.names = NULL)
  colnames(train) <- c(paste("X",1:ncol(data), sep=""), "metadata")
  
  
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parLapply(cl, 1:nb_train, function(x){
    
    simu <- simulation_supervise(apprent, nb_sample = nb_sample)
    
    if(type=="MAP") simu$data <- simu$data %>% MAP()
    if(type=="ilr") simu$data <- simu$data %>% MAP() %>% ilr()
    
    test <- data.frame(simu$data, row.names = NULL)
    colnames(test) <- paste("X",1:ncol(test), sep="")
    
    find_group(train, test, simu$metadata)
      
   })
  
  stopCluster(cl)
  
  
  sum_table <- list(random_forest=0, kNN=0, Svm=0)
  
  for(i in 1:length(l)){
    for(j in 1:length(l[[i]])){
      sum_table[[j]] <- sum_table[[j]]+l[[i]][[j]]
    }
  }
  
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
  table_group_kNN <- matrix(0, nrow=length(name_group), ncol=length(name_group))
  colnames(table_group_kNN) <- 1:length(name_group)
  rownames(table_group_kNN) <- name_group
  
  table_group_Svm <- matrix(0, nrow=length(name_group), ncol=length(name_group))
  colnames(table_group_Svm) <- 1:length(name_group)
  rownames(table_group_Svm) <- name_group
  
  train <- data.frame(data, as.factor(metadata), row.names = NULL)
  colnames(train) <- c(paste("X",1:ncol(data), sep=""), "metadata")
  
  
  iid <- sample(rep(1:k,length=nrow(data)))
  
  
  
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parLapply(cl, 1:k, function(x){
    
    kNN <- knn(train=data[x!=iid,], test=data[x==iid,], cl=metadata[x!=iid])
    Svm <- svm(metadata~., train[x!=iid,])
    pred_Svm <-  predict(Svm, train[x==iid,])
    
    mat <- matrix(0, nrow=length(name_group), ncol=length(name_group))
    mat2 <- matrix(0, nrow=length(name_group), ncol=length(name_group))
    
    for(i in 1:length(kNN)){
      mat[which(name_group==kNN[i]), which(name_group==metadata[x==iid][i])] <- mat[which(name_group==kNN[i]), which(name_group==metadata[x==iid][i])] + 1
      mat2[which(name_group==pred_Svm[i]), which(name_group==metadata[x==iid][i])] <- mat2[which(name_group==pred_Svm[i]), which(name_group==metadata[x==iid][i])] + 1
    }
    
    list(mat, mat2)
  })
  
  stopCluster(cl)
  
  
  for(i in 1:length(l)){
    table_group_kNN <- table_group_kNN+l[[i]][[1]]
    table_group_Svm <- table_group_Svm+l[[i]][[2]]
  }
  
  list(kNN=table_group_kNN, Svm=table_group_Svm)
}


find_real <- function(data_real, data_simu_train, data_simu){
  
  data_train <- rbind(data_real, data_simu_train)
  metadata <- c(rep("real", nrow(data_real)), rep("simu", nrow(data_simu_train))) %>% as.factor()
  train <- data.frame(data_train, metadata ,row.names = NULL)
  colnames(train) <- c(paste("X",1:ncol(data_train), sep=""), "metadata")
  
  test <- data.frame(data_simu, row.names = NULL)
  colnames(test) <- paste("X",1:ncol(test), sep="")

  
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
  
  apprent <- apprentissage_supervise(data, metadata)
  
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parLapply(cl, 1:nb_train, function(x){
    simu <- simulation_supervise(apprent, nb_sample=10*find_nb_sample_cluster(metadata = metadata))
  
    s <- sample(1:nrow(simu$data), nrow(data))
    data_simu_train <- simu$data[s, ]
    data_simu <- simu$data[-s, ]
    metadata_simu <- simu$metadata[-s]
  
    iid_real <- find_real(data_real = data, data_simu_train = data_simu_train, data_simu = data_simu)
  
    data_real <- data_simu[iid_real, ]
    metadata_real <- metadata_simu[iid_real]


    train <- data.frame(data, as.factor(metadata), row.names = NULL)
    colnames(train) <- c(paste("X",1:ncol(data), sep=""), "metadata")
  
    test <- data.frame(data_real, row.names = NULL)
    colnames(test) <- paste("X",1:ncol(test), sep="")
  
  
    find_group(train=train, test = test, metadata_test =  metadata_real)
  })
  
  stopCluster(cl)
  
  
  sum_table <- list(random_forest=0, kNN=0, Svm=0)
  
  for(i in 1:length(l)){
    for(j in 1:length(l[[i]])){
      sum_table[[j]] <- sum_table[[j]]+l[[i]][[j]]
    }
  }
  
  list(all=sum_table, all_train=l)
}


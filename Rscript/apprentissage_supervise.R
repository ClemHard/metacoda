library(randomForest)
library(class)
library(e1071)
source("Rscript/bootstrap.R")


apprentissage_supervise <- function(data, metadata){
  
  name_group <- unique(metadata)
  apprent <- lapply(name_group, function(x){
    iid <- which(metadata==x)
    apprentissage(data[iid,])
  })
  apprent
}


simulation_supervise <- function(apprent){
  data <- NULL
  metadata <- NULL
  for(i in 1:length(apprent)){
    data <- rbind(data, simulation(apprent[[i]]))
    metadata <- c(metadata, rep(i, nrow(apprent[[i]]$data)))
  }
  
  list(data=data, metadata=metadata)
}


bootstrap_supervise <- function(data, metadata){
  apprent <- apprentissage_supervise(data=data, metadata = metadata)
  simulation_supervise(apprent)
}



test_bootstrap_supervise <- function(data, metadata, nb_train=1, type="comptage"){
  
  apprent <- apprentissage_supervise(data=data, metadata = metadata)
  
  if(type=="MAP") data <- data %>% MAP()
  if(type=="ilr") data <- data %>% MAP() %>% ilr()
  
  train <- data.frame(data, as.factor(metadata), row.names = NULL)
  colnames(train) <- c(paste("X",1:ncol(data), sep=""), "metadata")
  
  l <- lapply(1:nb_train, function(x){
    
    simu <- simulation_supervise(apprent)
    
    if(type=="MAP") simu$data <- simu$data %>% MAP()
    if(type=="ilr") simu$data <- simu$data %>% MAP() %>% ilr()
    
    test <- data.frame(simu$data, row.names = NULL)
    colnames(test) <- paste("X",1:ncol(test), sep="")
    
    
    ### random forest
    forest <- randomForest(metadata~., train)
    pred_forest <- predict(forest, test)
    
    ### kNN
    kNN <- knn(data, simu$data, cl=metadata, k=1)
    
    ### Svm
    Svm <- svm(metadata~., train)
    pred_svm <- predict(Svm, test)
    
    list(random_forest=table(pred_forest, simu$metadata), kNN=table(kNN, simu$metadata), Svm=table(pred_svm, simu$metadata))
  })
  
  
  sum_table <- list(random_forest=0, kNN=0, Svm=0)
  
  for(i in 1:length(l)){
    for(j in 1:length(l[[i]])){
      sum_table[[j]] <- sum_table[[j]]+l[[i]][[j]]
    }
  }
  
  list(all=sum_table, all_train=l)
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
  clusterEvalQ(cl, {library(randomForest); library(e1071); library(dplyr); source("Rscript/apprentissage_supervise.R")})
  
  
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


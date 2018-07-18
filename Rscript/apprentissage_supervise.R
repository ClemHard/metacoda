library(randomForest)

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



test_bootstrap_supervise <- function(data, metadata, nb_train=1){
  
 apprent <- apprentissage_supervise(data=data, metadata = metadata)
 
 lapply(1:nb_train, function(x){
   
    simu <- simulation_supervise(apprent)
    train <- data.frame(data, as.factor(metadata), row.names = NULL)
    colnames(train) <- c(paste("X",1:ncol(data), sep=""), "metadata")
    
    test <- data.frame(simu$data, row.names = NULL)
    colnames(test) <- paste("X",1:ncol(test), sep="")
    
    
    ### random forest
    forest <- randomForest(metadata~., train)
    pred <- predict(forest, test)
    
    ### kNN
    kNN <- knn(data, simu$data, cl=metadata, k=1)
    
    
    list(random_forest=table(pred, simu$metadata), kNN=table(kNN, simu$metadata))
 })
  
}



validation_croise <- function(data, metadata, k=nrow(data)){
  
  name_group <- unique(metadata)
  table_group <- matrix(0, nrow=length(name_group), ncol=length(name_group))
  colnames(table_group) <- 1:length(name_group)
  rownames(table_group) <- name_group
  
  iid <- sample(rep(1:k,length=nrow(data)))
  
  l <- lapply(1:k, function(x){
    
    kNN <- knn(train=data[x!=iid,], test=data[x==iid,], cl=metadata[x!=iid])
    mat <- matrix(0, nrow=length(name_group), ncol=length(name_group))
    
    for(i in 1:length(kNN)){
      mat[which(name_group==kNN[i]), which(name_group==metadata[x==iid][i])] <- mat[which(name_group==kNN[i]), which(name_group==metadata[x==iid][i])] + 1
    }
    
    mat
    })
  
  
  for(i in 1:length(l)){
    table_group <- table_group+l[[i]]
  }
  table_group
}



leave_one_out <- function(data, metadata, nb_train=nrow(data)){
  
  
  lapply(1:nb_train, function(x){
    knn(train=data[-x,], test=data[x,], cl=metadata[-x])==metadata[x]  
  })

  }
library(randomForest)
library(class)
library(e1071)
library(dplyr)
library(doParallel)

source("Rscript/bootstrap.R")
source("Rscript/test_bootstrap.R")

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
  k <- best_k_kNN(train)
  kNN <- knn(train[, -ncol(train)], test, cl=train[, ncol(train)], k=k)
  
  ### Svm
  Svm <- svm(metadata~., train)
  pred_svm <- predict(Svm, test)
  
  ### neuronal
 # nn <- nnet(metadata~., train, size=5, MaxNWts=1000000, maxit=500)
 # pred_nn <- predict(nn, test, type="class")
  
  list(random_forest=table(pred_forest, metadata_test), kNN=table(kNN, metadata_test), Svm=table(pred_svm, metadata_test))
}




validation_croise <- function(data, metadata, n=nrow(data)){
  
  name_group <- unique(metadata)
  train <- create_data_frame_train(data, metadata)
  
  iid <- sample(rep(1:n,length=nrow(data)))
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  l <- parLapply(cl, 1:n, function(x){
    
    k <- best_k_kNN(train=train)
    pred_kNN <- knn(train=data[x!=iid,], test=data[x==iid,], cl=metadata[x!=iid], k=k)
    Svm <- svm(metadata~., train[x!=iid,])
    pred_Svm <-  predict(Svm, train[x==iid,-ncol(train)])
    
    table1 <- table(metadata[x==iid], pred_kNN)
    table2 <- table(metadata[x==iid], pred_Svm)
    colnames(table1) <- 1:length(name_group)
    colnames(table2) <- 1:length(name_group)
    
    list(table1, table2)
  })
  
  stopCluster(cl)
  
  all <- list_table_sum(l)
  names(all) <- c("kNN", "Svm")
  all$confusion_random_Forest <- randomForest(metadata~., create_data_frame_train(data, metadata))$confusion[,-(length(unique(metadata))+1)]
  colnames(all$confusion_random_Forest) <- 1:(length(unique(metadata)))
  all
}



find_real <- function(data_real, data_simu_train, data_simu, algo="randomForest"){
  
  data_train <- rbind(data_real, data_simu_train)
  metadata <- c(rep("real", nrow(data_real)), rep("simu", nrow(data_simu_train))) %>% as.factor()

  train <- create_data_frame_train(data=data_train, metadata = metadata)
  test <- create_data_frame_test(data=data_simu)

  pred <- NULL
  
  
  if(algo=='randomForest'){
    forest <- randomForest(metadata~., train)
    pred <- predict(forest, test)
  }

  if(algo=="kNN"){
    k <- best_k_kNN(train)
    pred <- knn(train[,-ncol(train)], 
                    test, 
                    cl=train[,ncol(train)], 
                    k=train)
  }
  
  if(algo=="logistic"){
    logi <- glm(metadata~., family = binomial, data=train, control = list(maxit = 200))
    pred <- predict(logi, test)>0.5
  }
  
  if(algo=="neural"){
    nn <- nnet(metadata~., train, size=5, MaxNWts=1000000, maxit=500)
    pred <- predict(nn, test, type="class")
  }
  
  iid <- which(pred=="real")
  
  iid
}



find_nb_sample_cluster <- function(metadata){
  name_group <- unique(metadata)
  sapply(name_group, function(x){sum(x==metadata)})
}



classificateur_group_real_simu <- function(data, metadata, nb_train=1, type="comptage"){
  
  
  boot <- bootstrap(data = data)
  
  apprent <- apprentissage_supervise(data=data, metadata=metadata)
  
  if(type=="ilr"){
    data <- data %>% MAP() %>% ilr()
    boot$data <- boot$data %>% MAP() %>% ilr()
  }
  
  if(type=="MAP"){
    data <- data %>% MAP()
    boot$data <- boot$data %>% MAP()
  }
  
  
  train <- create_data_frame_train(data = data, metadata=metadata)

  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <-  parLapply(cl, 1:nb_train, function(x){
    
    simu <- simulation_supervise(apprent, nb_sample=10*find_nb_sample_cluster(metadata = metadata))
  
    if(type=="ilr"){
      simu$data <- simu$data %>% MAP() %>% ilr()
    }
    
    if(type=="MAP"){
      simu$data <- simu$data %>% MAP()
    }
    
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


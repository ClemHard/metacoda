library(randomForest)
library(class)
library(e1071)
library(dplyr)
library(doParallel)

source("Rscript/bootstrap.R")
source("Rscript/test_bootstrap.R")


#'  fit a gaussian mixture on each group of a dataset (OTU in column, sample in line)
#'  
#' @param data a dataset (count data)
#' @param nb_cluster the number of gaussian in the gaussian mixture
#' @param  base_binaire binary sequential matrix use for the ilr transformation (default the basis sequence)
#'
#' @return an object of class "apprentissage supervise"
#' 
#' @author Clement Hardy
#' @export

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


#' produce news sample (compositionnal data) using an object of class "apprentissage supervise"
#' after simulated compositionnal data of each group, the method use a zero inflated to force zero.
#' the proportion of zero in news samples is determine by the zero inflated function
#' 
#'
#' @param apprent an object of class "apprentissage supervise"
#' @param nb_sample a numerice vector containing the number of sample of each group to simulate
#' 
#' @return the news samples
#' @return the group of each sample
#'
#' @author Clement Hardy
#' @export
#' @import stats


simulation_supervise <- function(apprent, nb_sample=NULL){
  
  data <- NULL
  metadata <- NULL
  if(length(nb_sample)!=length(unique(apprent)) || any(nb_sample<0)){
    nb_sample <- NULL
  }
  for(i in 1:length(apprent)){
    if(is.null(nb_sample)){
      data <- rbind(data, simulation(apprent[[i]]))
      metadata <- c(metadata, rep(apprent[[i]]$metadata, nrow(apprent[[i]]$data)) %>% as.character())
    }else{
      data <- rbind(data, simulation(apprent[[i]], nb_sample = nb_sample[i]))
      metadata <- c(metadata, rep(apprent[[i]]$metadata, nb_sample[i]) %>% as.character())
    }
    
  }
  
  list(data=data, metadata=metadata)
}

#' produces news samples (with different group) from a dataset of OTU table,
#' for each group of the dataset the following method is apply
#' the method is a base on compositional data with isometric log ratio transformation 
#' a dimension reduction and a gaussian mixture is used to produce the news samples
#' the dimension of the latent space could be given in argument or by default find with the heuristic slope.
#' the number of gaussian in the gaussian mixture could be given in argument or default find with the bic criteria
#'
#' @param data the dataset to simulate
#' @param metadata a numeric vector containing the name of the group 
#' @param nb_sample a numeric vector containing the number of sample to produce for each group
#'
#' @return a list containing a matrix with the produce sample, an object of class "apprentissage", the dimension of the latent space, the number of gaussian
#' 
#' @author Clement Hardy
#' @export


bootstrap_supervise <- function(data, metadata, nb_sample=NULL){
  apprent <- apprentissage_supervise(data=data, metadata = metadata)
  simulation_supervise(apprent, nb_sample = nb_sample)
}



#' test the performance of the simulator for a count dataset 
#' Firstly, the function learn the density of the different groups (of the dataset) with a gaussian mixture and 
#' different classification algorithm like randomForest are train in the dataset
#' 
#' Then, nb_train new dataset are simulated ; for each new dataset the different classification algoritm predict the class
#' of each new sample (of the new dataset) , the result of the prediction is then compare to the real class to create
#' a contingence table.
#' 
#' @param data the dataset to simulate for the test 
#' @param metadata a numeric vector containing the name of the group 
#' @param nb_sample a numeric vector containing the number of sample to produce for each group (default nrow(data))
#'
#' @return a list containing the mean of the nb_train contingence table (for each incompatible) and a list of all contigence table produces (for each algorithm)
#' 
#' @author Clement Hardy
#' @export
#' @import doParallel

test_bootstrap_supervise <- function(data, metadata, nb_train=1, type="comptage", nb_sample=NULL){
  
  apprent <- apprentissage_supervise(data=data, metadata = metadata)
  
  if(type=="MAP") data <- data %>% MAP()
  if(type=="ilr") data <- data %>% MAP() %>% ilr()
  
  train <- create_data_frame_train(data=data, metadata = metadata)
  
  
  no_cores <- parallel::detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- parallel::makeCluster(no_cores)
  parallel::clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <- parallel::parLapply(cl, 1:nb_train, function(x){
    
    simu <- simulation_supervise(apprent, nb_sample = nb_sample)
    
    if(type=="MAP") simu$data <- simu$data %>% MAP()
    if(type=="ilr") simu$data <- simu$data %>% MAP() %>% ilr()
    
    test <- create_data_frame_test(data=simu$data)
    
    find_group(train, test, simu$metadata)
    
  })
  
  parallel::stopCluster(cl)
  
  sum_table <- lapply(list_table_sum(l), table_to_percentage_table)
  list(all=sum_table, all_train=l)
}

#' test the performance of the simulator for a count dataset 
#' Firstly, the function learn the density of the different groups (of the dataset) with a gaussian mixture and 
#' different classification algorithm like randomForest are train in the dataset
#' 
#' Then, nb_train new dataset are simulated ; for each new dataset the different classification algoritm predict the class
#' of each new sample (of the new dataset) , the result of the prediction is then compare to the real class to create
#' a contingence table.
#' 
#' @param data the dataset to simulate for the test 
#' @param metadata a numeric vector containing the name of the group 
#' @param nb_sample a numeric vector containing the number of sample to produce for each group (default nrow(data))
#'
#' @return a list containing the mean of the nb_train contingence table (for each incompatible) and a list of all contigence table produces (for each algorithm)
#' 
#' @author Clement Hardy
#' @export
#' @import randomForest
#' @import e1071
#' @import class

find_group <- function(train, test, metadata_test){
  ### random forest
  forest <- randomForest::randomForest(metadata~., train)
  pred_forest <- predict(forest, test)
  
  ### kNN
  k <- best_k_kNN(train)
  kNN <- class::knn(train[, -ncol(train)], test, cl=train[, ncol(train)], k=k)
  
  ### Svm
  Svm <- e1071::svm(metadata~., train)
  pred_svm <- predict(Svm, test)
  
  ### neuronal
  # nn <- nnet(metadata~., train, size=5, MaxNWts=1000000, maxit=500)
  # pred_nn <- predict(nn, test, type="class")
  
  list(random_forest=table(pred_forest, metadata_test) %>% change_name_table(), kNN=table(kNN, metadata_test) %>% change_name_table(), Svm=table(pred_svm, metadata_test) %>% change_name_table())
}



#' perform the cross validation for tow algorithms (kNN, SVM) and return the confusion matrix for the random Forest
#' 
#' @param data the dataset use for the cross validation
#' @param metadata a numeric vector containing the name of the group 
#' @param n the dataset is separated in n groups for the cross validation (by default n=nrow(data), leave one out is perform)
#'
#' @return contigences table of cross validation (for each algorithm)
#' 
#' @author Clement Hardy
#' @export
#' @import parallel
#' @import e1071
#' @import class
#' @import randomForest



validation_croise <- function(data, metadata, n=nrow(data)){
  
  name_group <- unique(metadata)
  train <- create_data_frame_train(data, metadata)
  
  iid <- sample(rep(1:n,length=nrow(data)))
  
  no_cores <- parallel::detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- parallel::makeCluster(no_cores)
  parallel::clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  l <- parallel::parLapply(cl, 1:n, function(x){
    
    k <- best_k_kNN(train=train)
    pred_kNN <- class::knn(train=data[x!=iid,], test=data[x==iid,], cl=metadata[x!=iid], k=k)
    Svm <- e1071::svm(metadata~., train[x!=iid,])
    pred_Svm <-  predict(Svm, train[x==iid,-ncol(train)])
    
    table1 <- table(pred_kNN, metadata[x==iid])
    table2 <- table(pred_Svm, metadata[x==iid])
    
    
    list(table1, table2)
  })
  
  parallel::stopCluster(cl)
  
  all <- list_table_sum(l)
  names(all) <- c("kNN", "Svm")
  all$confusion_random_Forest <- randomForest::randomForest(metadata~., create_data_frame_train(data, metadata))$confusion[,-(length(unique(metadata))+1)]
  all <- lapply(all, change_name_table)
  all
}

#' find the simulated sample pass in argument which are considered as real sample 
#' 
#' @param data_real the dataset of real sample use to train the classification algorithm
#' @param data_simu_train a dataset of simulated sample use to train the classification algorithm
#' @param data_simu a dataset of simulated sample (the function return the index of the sample considered as real of this dataset)
#' @param algo the classification algorithm use to find the real sample (by defautl randomForest), could be randomForest, logistic, kNN or neural (nnet)
#' 
#' @return a vector of the index of the sample considere as real
#' 
#' @author Clement Hardy
#' @export
#' @import randomForest


find_real <- function(data_real, data_simu_train, data_simu, algo="randomForest"){
  
  data_train <- rbind(data_real, data_simu_train)
  metadata <- c(rep("real", nrow(data_real)), rep("simu", nrow(data_simu_train))) %>% as.factor()
  
  train <- create_data_frame_train(data=data_train, metadata = metadata)
  test <- create_data_frame_test(data=data_simu)
  
  pred <- NULL
  
  
  if(algo=='randomForest'){
    forest <- randomForest::randomForest(metadata~., train)
    pred <- predict(forest, test)
  }
  
  if(algo=="kNN"){
    k <- best_k_kNN(train)
    pred <- class::knn(train[,-ncol(train)], 
                test, 
                cl=train[,ncol(train)], 
                k=train)
  }
  
  if(algo=="logistic"){
    logi <- glm(metadata~., family = binomial, data=train, control = list(maxit = 200))
    pred <- predict(logi, test)>0.5
  }
  
  if(algo=="neural"){
    nn <- nnet::nnet(metadata~., train, size=5, MaxNWts=1000000, maxit=500)
    pred <- predict(nn, test, type="class")
  }
  
  iid <- which(pred=="real")
  
  iid
}


#' return the number of sample in each group of a dataset 
#' 
#' @param metadata a numeric vector containing the name of the group 
#' 
#' @return a vector containing the number of sample for each group
#' 
#' @author Clement Hardy
#' @export


find_nb_sample_cluster <- function(metadata){
  name_group <- unique(metadata)
  sapply(name_group, function(x){sum(x==metadata)})
}



#' test the performance of the simulator for a count dataset 
#' The different with the function test_bootstrap supervise is that this time only the simulated sample considered as real
#' are keep for the test of classification
#' 
#' Firstly, the function learn the density of the different groups (of the dataset) with a gaussian mixture and 
#' different classification algorithm like randomForest are train in the dataset
#' 
#' Then, nb_train new dataset are simulated ; for each new dataset the different classification algoritm predict the class
#' of each new sample considered as real by the function find_real (of the new dataset) , the result of the prediction is then compare to the real class to create
#' a contingence table.
#' 
#' @param data the dataset to simulate for the test 
#' @param metadata a numeric vector containing the name of the group 
#' @param nb_sample a numeric vector containing the number of sample to produce for each group (default nrow(data))
#'
#' @return a list containing the mean of the nb_train contingence table (for each incompatible) and a list of all contigence table produces (for each algorithm)
#' 
#' @author Clement Hardy
#' @export
#' @import parallel



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
  
  no_cores <- parallel::detectCores()-1
  if(no_cores==0) no_cores <- 1
  
  cl <- parallel::makeCluster(no_cores)
  parallel::clusterEvalQ(cl, {source("Rscript/apprentissage_supervise.R")})
  
  
  l <-  parallel::parLapply(cl, 1:nb_train, function(x){
    
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
  
  parallel::stopCluster(cl)
  
  sum_table <- list_table_sum(l)
  
  
  sum_table <- lapply(sum_table, table_to_percentage_table)
  list(all=sum_table, all_train=l)
}

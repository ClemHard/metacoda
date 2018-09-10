library(randomForest)
library(doParallel)
library(dplyr)
library(class)
library(e1071)
library(nnet)


source("Rscript/bootstrap.R")


#' make sure the two table are compatible(same rownames, same colnames)
#'In case the rownames, colnames aren't the same, the function add column, row to the tables
#'
#' @param table1,table2 the two tables to check
#' 
#' @return a list containing the tow tables check and modify (if necessary)
#'
#' @author Clement Hardy
#' @export
#' 


check_table <- function(table1, table2){

  c <- colnames(table1) %in% colnames(table2)
  if(any(c==FALSE)){
    name <- colnames(table1)[which(c==FALSE)]
    temp <- matrix(0, ncol=length(name), nrow=nrow(table1))
    colnames(temp) <- name
    table2 <- cbind(table2, temp) %>% as.table()
  }
  
  c <- colnames(table2) %in% colnames(table1)
  if(any(c==FALSE)){
    name <- colnames(table2)[which(c==FALSE)]
    temp <- matrix(0, ncol=length(name), nrow=nrow(table2))
    colnames(temp) <- name
    table1 <- cbind(table1, temp) %>% as.table()
  }
  
  c <- rownames(table1) %in% rownames(table2)
  if(any(c==FALSE)){
    name <- rownames(table1)[which(c==FALSE)]
    temp <- matrix(0, nrow=length(name), ncol=ncol(table1))
    rownames(temp) <- name
    table2 <- rbind(table2, temp) %>% as.table()
  }
  
  c <- rownames(table2) %in% rownames(table1)
  if(any(c==FALSE)){
    name <- rownames(table2)[which(c==FALSE)]
    temp <- matrix(0, nrow=length(name), ncol=ncol(table2))
    rownames(temp) <- name
    table1 <- rbind(table1, temp) %>% as.table()
  }
  
  row_name <- (rownames(table1))
  
  if(ncol(table1)>1){
    table1 <- table1[ ,order(colnames(table1))]
    table2 <- table2[ ,order(colnames(table2))]
  }
  
  if(class(table1)=="table" && class(table2)=="table"){
    table1 <- table1[order(rownames(table1)), ]
    table2 <- table2[order(rownames(table2)), ]
  }
  
  if(class(table1)=="integer"){
    table1 <- table1 %>% as.matrix %>% t()
    table2 <- table2 %>% as.matrix() %>% t()
    rownames(table1) <- row_name
    rownames(table2) <- row_name
  }
  
  list(table1, table2)
}



#' Modify a table to a table with percentage
#'
#' @param table the table to modify to a percentage table
#' 
#' @return the percentage table
#'
#' @author Clement Hardy
#' @export


table_to_percentage_table <- function(table){
  nb_row <- nrow(table)
  table <- round(table/matrix(rep(apply(table, 2, sum), nb_row), byrow = TRUE, nrow=nb_row)*100, digits=2)
  table
}


#' add predict or observe to the colnames and rownames of a table 
#'
#' @param table 
#' 
#' @return the table with the colnames and rownames modify
#'
#' @author Clement Hardy
#' @export


change_name_table <- function(table){
  table <- table[,order(colnames(table))]
  table <- table[order(rownames(table)), ]
  colnames(table) <- colnames(table) %>% paste("(obs)", sep="")
  rownames(table) <- rownames(table) %>% paste("(pred)", sep="")
  table
}



list_table_sum <- function(l){
  
  if(class(l)!="list") stop("l isn't a list")
    
  if(length(l)==0){
    warning("list is empty")
    return(NULL)
  }
  if(length(l)==1){
    return(l[[1]])
  }
  
  sum_table <- l[[1]]
  
  for(i in 2:length(l)){
    for(j in 1:length(l[[i]])){
      temp <- check_table(sum_table[[j]], l[[i]][[j]])
      sum_table[[j]] <- temp[[1]]
      l[[i]][[j]] <- temp[[2]]
      sum_table[[j]] <- sum_table[[j]]+l[[i]][[j]]
    }
  }
  
  sum_table
}


#' create a test data frame (to make sure the colnames are compatible with the test function) with the data (matrix) pass in argument
#'
#' @param data a matrix containing the test data
#' 
#' @return a data.frame 
#'
#' @author Clement Hardy
#' @export
#' 

create_data_frame_test <- function(data){
  test <- data.frame(data, row.names = NULL)
  colnames(test) <- paste("X",1:ncol(test), sep="")
  test
}


#' create a train data frame (to make sure the colnames are compatible with the test function) with the data (matrix) pass in argument
#' the metadata is in the last column of the data frame
#' 
#' @param data a matrix containing the test data
#' 
#' @return a data.frame 
#'
#' @author Clement Hardy
#' @export
#' 


create_data_frame_train <- function(data, metadata){
  train <- create_data_frame_test(data)
  train$metadata <- metadata %>% as.factor()
  train
}


#' create a test data frame (to make sure the colnames are compatible with the test function) with the data (matrix) pass in argument
#' 
#' @param data1 a matrix containing the real sample
#' @param apprent an object of class "apprentissage" (create by the function apprentissage)
#' @param test_real_data a matrix containing the real sample use for the test
#' @param type comptage or MAP
#' @param zero_inflated boolean, if TRUE a zero inflated method is apply, FALSE otherwise
#'  
#'    
#' @return train a training data frame
#' @return test a test data frame
#' @return metadata_test the real metadata of the test sample (necessary to create a contingence table after applying a classification algorithm)
#'
#' @author Clement Hardy
#' @export


create_data_frame_simu <- function(data1, apprent, test_real_data=NULL, type, zero_inflated){
  
  boot <- simulation(apprent, nb_sample = nrow(data1), type = "comptage", zero_inflated = zero_inflated)
  nb_sample <- 1000
  boot_test <- simulation(apprent, nb_sample = nb_sample, type = "comptage", zero_inflated = zero_inflated)
  
  
  if(type=="MAP") {
    boot_test <- boot_test %>% MAP()
    boot <- boot %>% MAP()
  }
  
  if(type=="ilr"){
    boot_test <- boot_test %>% MAP() %>% ilr()
    boot <- boot %>% MAP() %>% ilr()
  }
  
  
  data_train <- rbind(data1, boot)
  metadata <- c(rep("real", nrow(data1)), rep("simu", nrow(boot)))
  
  train <- create_data_frame_train(data=data_train, metadata = metadata)
  
  test <- create_data_frame_test(rbind(boot_test, test_real_data))
  colnames(test_real_data) <- paste("X",1:ncol(test), sep="")

  metadata_test <- c(rep("simu", nrow(boot_test)), rep("real", nrow(test_real_data)))

  list(train=train, test=test, metadata_test=metadata_test) 
}


#' select the optimal number of nearest neighbors in the kNN  with a 10 folds cross validation
#' 
#' @param train a data frame created by the function create_data_frame_train
#' 
#' @return the optimal k
#'
#' @author Clement Hardy
#' @export
#' @import class
#' 


best_k_kNN <- function(train){
  
  metadata <- train[,ncol(train)]
  train <- train[,-ncol(train)]
  
  k_max <- (nrow(train)/5) %>% ceiling()
  n <- (nrow(train)/10) %>% ceiling()
  
  iid <- sample(rep(1:n,length=nrow(train)))
  
  l <- sapply(1:n, function(x){
                          pred_kNN <- class::knn(train=train[x!=iid,], 
                                          test=train[x==iid,], 
                                          cl=metadata[x!=iid])
                         sum(pred_kNN!=metadata[x==iid]) 
                })
  which.min(l)
}


#' three differents classification algorithm which are random Forest, logistic regression and k nearest neighbors
#' predict the type (real or simulated) of the sample pass in argument. Firstly the classificators are train in a dataset
#' containing real and simulated sample
#' 
#' @param train a data frame created by the function create_data_frame_train
#' @param test a data frame created by the function create_data_frame_test
#' @return list containing prediction (a vector) of each classification algorithm
#'
#' @author Clement Hardy
#' @export
#' @import class
#' @import randomForest



pred_real_simu <- function(train, test){
  
  rf <- randomForest::randomForest(metadata~. , train)
  pred_forest <- predict(rf, test)
  
  k <- best_k_kNN(train)
  pred_kNN <- class::knn(train[,-ncol(train)], 
                  test, 
                  cl=train[,ncol(train)], 
                  k=k)
  
  logi <- glm(metadata~., family = binomial, data=train, control = list(maxit = 200))
  pred_logi <- (predict(logi, test)>0.5)
  
  #nn <- nnet(metadata~., train, size=5, MaxNWts=1000000, maxit=500)
  #pred_nnet <- predict(nn, test, type="class")
  
  
  list(random_forest=pred_forest, kNN=pred_kNN, logi=pred_logi)
}

#' test the performance of the simulator for a count dataset 
#' three differents classification algorithm are use in the goal to predict the type (real or simualted) of simulated sample
#' A good simulation lead to an impossibility for the classification algorithm to see
#' the different beetween a simulated or a real sample.
#' 
#' For i in 1,2,..., nb_train, this function simulate two new dataset. One is with the real dataset (some samples aren't use in the training dataset but are use for the test) use to train the classificators (to learn the differents beetween the real and simulated sample)
#' The classificators are then use to predict the type of the sample of the second simulated dataset (some real samples are also use in the test dataset), a table contingence are finally created.
#' 
#' 
#' @param data1 the dataset (containing the real samples) to simulate for the test 
#' @param metadata a numeric vector containing the name of the group (the group of the samples of the dataset data1) 
#' @param nb_cluster the number of gaussian in the gaussian mixture
#' @param nb_axe the dimension of the latent space (dimension reduction)
#' @param nb_sample the number of sample to produce
#' @param  base_binaire binary sequential matrix use for the ilr transformation
#' @param zero_inflated boolean, if TRUE a zero inflated method is apply (on the simulated dataset), FALSE otherwise
#' @param proportion_real_data the proportion of real samples use in the test dataset (should be a number beetween 0 and 1) 
#' 
#' @return a list containing the mean of the nb_train contingence table (for each incompatible) and a list of all contigence table produces (for each algorithm)
#' 
#' @author Clement Hardy
#' @export
#' @import parallel


test_bootstrap_all <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_train=1, proportion_real_data=0.1, type="comptage", zero_inflated=TRUE, base_binaire=Base_binary_matrix(ncol(data1))){
  
  
  data2 <- data1

  if(type=="ilr") data2 <- data1 %>% MAP() %>% ilr(base_binaire = base_binaire)
  if(type=="MAP") data2 <- data1 %>% MAP()
  
  no_cores <- parallel::detectCores()-1
  if(no_cores==0) no_cores <- 1
    
  cl <- parallel::makeCluster(no_cores)
  parallel::clusterEvalQ(cl, {source("Rscript/test_bootstrap.R")})
  parallel::clusterExport(cl, list("data1", "data2", "nb_cluster", "nb_axe", "base_binaire"), envir = environment())
  
  l <- parallel::parLapply(cl, 1:nb_train, function(x){
    
                            rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*proportion_real_data))
                            test_real_data <- data2[rand_real_data,]
                            
                            data3 <- data1[-rand_real_data,]
                            data_train <- data2[-rand_real_data,]
                            apprent <- apprentissage(data3, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)

                            
                            data_test <- create_data_frame_simu(data_train, apprent=apprent, test_real_data = test_real_data, type=type, zero_inflated=zero_inflated)
                    
                            
                            pred <- pred_real_simu(data_test$train ,data_test$test)

                            
                            lapply(pred, function(x){x %>% table(data_test$metadata) %>% change_name_table()})
                            
                            })
  
  parallel::stopCluster(cl)
  print(l)
  all <- l %>% list_table_sum()
  
  all <- lapply(all, table_to_percentage_table)
  
  list(all=all, all_train=l)
}




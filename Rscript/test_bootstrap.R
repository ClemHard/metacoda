library(randomForest)
library(doParallel)
library(dplyr)
library(class)
library(e1071)
library(nnet)


source("Rscript/bootstrap.R")



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



table_to_percentage_table <- function(table){
  nb_row <- nrow(table)
  table <- round(table/matrix(rep(apply(table, 2, sum), nb_row), byrow = TRUE, nrow=nb_row)*100, digits=2)
  table
}


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


create_data_frame_test <- function(data){
  test <- data.frame(data, row.names = NULL)
  colnames(test) <- paste("X",1:ncol(test), sep="")
  test
}


create_data_frame_train <- function(data, metadata){
  train <- create_data_frame_test(data)
  train$metadata <- metadata %>% as.factor()
  train
}



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


best_k_kNN <- function(train){
  
  metadata <- train[,ncol(train)]
  train <- train[,-ncol(train)]
  
  k_max <- (nrow(train)/5) %>% ceiling()
  n <- (nrow(train)/10) %>% ceiling()
  
  iid <- sample(rep(1:n,length=nrow(train)))
  
  l <- sapply(1:n, function(x){
                          pred_kNN <- knn(train=train[x!=iid,], 
                                          test=train[x==iid,], 
                                          cl=metadata[x!=iid])
                         sum(pred_kNN!=metadata[x==iid]) 
                })
  which.min(l)
}


pred_real_simu <- function(train, test){
  
  rf <- randomForest(metadata~. , train)
  pred_forest <- predict(rf, test)
  
  k <- best_k_kNN(train)
  pred_kNN <- knn(train[,-ncol(train)], 
                  test, 
                  cl=train[,ncol(train)], 
                  k=k)
  
  logi <- glm(metadata~., family = binomial, data=train, control = list(maxit = 200))
  pred_logi <- (predict(logi, test)>0.5)
  
  #nn <- nnet(metadata~., train, size=5, MaxNWts=1000000, maxit=500)
  #pred_nnet <- predict(nn, test, type="class")
  
  
  list(random_forest=pred_forest, kNN=pred_kNN, logi=pred_logi)
}


test_bootstrap_all <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_train=1, proportion_real_data=0.1, type="comptage", zero_inflated=TRUE, base_binaire=Base_binary_matrix(ncol(data1))){
  
  
  data2 <- data1

  if(type=="ilr") data2 <- data1 %>% MAP() %>% ilr(base_binaire = base_binaire)
  if(type=="MAP") data2 <- data1 %>% MAP()
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
    
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("Rscript/test_bootstrap.R")})
  clusterExport(cl, list("data1", "data2", "nb_cluster", "nb_axe", "base_binaire"), envir = environment())
  
  l <- parLapply(cl, 1:nb_train, function(x){
    
                            rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*proportion_real_data))
                            test_real_data <- data2[rand_real_data,]
                            
                            data3 <- data1[-rand_real_data,]
                            data_train <- data2[-rand_real_data,]
                            apprent <- apprentissage(data3, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)

                            
                            data_test <- create_data_frame_simu(data_train, apprent=apprent, test_real_data = test_real_data, type=type, zero_inflated=zero_inflated)
                    
                            
                            pred <- pred_real_simu(data_test$train ,data_test$test)

                            
                            lapply(pred, function(x){x %>% table(data_test$metadata) %>% change_name_table()})
                            
                            })
  
  stopCluster(cl)
  
  all <- l %>% list_table_sum()
  
  all <- lapply(all, table_to_percentage_table)
  
  list(all=all, all_train=l)
}




test_bootstrap_variable <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_train=1, proportion_real_data=0.1, type="comptage", base_binaire=Base_binary_matrix(ncol(data1))){
  
  
  if(type=="ilr") data1 <- data1 %>% MAP() %>% center_scale(scale=FALSE) %>% ilr(base_binaire = base_binaire)
  
  if(type=="MAP") data1 <- data1 %>% MAP()
  
  rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*proportion_real_data))
  test_real_data <- data1[rand_real_data,]
  data1 <- data1[-rand_real_data,]
  apprent <- apprentissage(data1, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)
  
  data_test <- create_data_frame_simu(data1, apprent=apprent, test_real_data = test_real_data, type=type)
  
  variables <- colnames(data_test$test)
  
  
  no_cores <- detectCores()-1
  if(no_cores>1){
    
    cl <- makeCluster(no_cores)
    clusterEvalQ(cl, {library(randomForest); library(dplyr)})
    clusterExport(cl, list("data_test"), envir = environment())
    
      each <- lapply(variables, function(x){
                                              rf <- randomForest(as.formula(paste("metadata~",x)), data_test$train)
                                              predict(rf, data_test$test) %>% table(data_test$metadata_test)
                                              })
      
    stopCluster(cl)
    
  }else{
    each <- lapply(variables, function(x){
                                      rf <- randomForest(as.formula(paste("metadata~",x)), data_test$train)
                                      predict(rf, data_test$test) %>% table(data_test$metadata_test)
                                       })
  }
  
  misclassification_each <- sapply(each, function(x){x[1, 2]})
  
  names(misclassification_each) <- colnames(data1)
  names(each) <- colnames(data1)
  
  
  list(table_each=each, misclassification=misclassification_each)
}




test_bootstrap <- function(data1, nb_cluster=NULL, nb_axe=NULL, proportion_real_data=0.1, nb_train_all=1, type="comptage", base_binaire=Base_binary_matrix(ncol(data1))){
  
  all <- test_bootstrap_all(data1, nb_train = nb_train_all, type=type, base_binaire=base_binaire, proportion_real_data = proportion_real_data)
  each <- test_bootstrap_variable(data1, type=type, base_binaire=base_binaire, proportion_real_data = proportion_real_data)
  
  list(all=all, table_each=each$table_each, misclassification=each$misclassification)
}



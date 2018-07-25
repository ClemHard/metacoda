library(randomForest)
library(doParallel)
library(dplyr)

source("Rscript/bootstrap.R")


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



create_data_frame_simu <- function(data1, apprent, test_real_data=NULL, type){
  
  
  boot <- simulation(apprent, nb_sample = nrow(data1), type = type)

  data_train <- rbind(data1, boot)
  metadata <- c(rep("real", nrow(data1)), rep("simu", nrow(boot)))
  
  nb_sample <- 1000
  boot_test <- simulation(apprent, nb_sample = nb_sample, type = type)
  
  train <- create_data_frame_train(data=data_train, metadata = metadata)
  
  test <- create_data_frame_test(rbind(boot_test, test_real_data))
  colnames(test_real_data) <- paste("X",1:ncol(test), sep="")

  metadata_test <- c(rep("simu", nrow(boot_test)), rep("real", nrow(test_real_data)))

  list(train=train, test=test, metadata_test=metadata_test) 
}



test_bootstrap_all <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1), apprent=NULL, nb_train=1, proportion_real_data=0.1, type="comptage", base_binaire=Base_binary_matrix(ncol(data1))){
  
  
  data2 <- data1
  
  if(type=="ilr") data2 <- data1 %>% MAP() %>% center_scale(scale=FALSE) %>% ilr(base_binaire = base_binaire)
  if(type=="MAP") data2 <- data1 %>% MAP()
  
  no_cores <- detectCores()-1
  if(no_cores==0) no_cores <- 1
    
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(randomForest); library(dplyr); source("Rscript/bootstrap.R"); source("Rscript/test_bootstrap.R")})
  
  all <- parLapply(cl, 1:nb_train, function(x){
    
                            rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*proportion_real_data))
                            test_real_data <- data2[rand_real_data,]
                            
                            data3 <- data1[-rand_real_data,]
                            data_train <- data2[-rand_real_data,]
                            apprent <- apprentissage(data3, nb_axe = nb_axe, nb_cluster = nb_cluster, base_binaire = base_binaire)

                            
                            data_test <- create_data_frame_simu(data_train, apprent=apprent, test_real_data = test_real_data, type=type)
                
                            rf <- randomForest(metadata~. , data_test$train)
                            pred <- predict(rf, data_test$test) %>% table(data_test$metadata_test)
                            list(pred=pred, rf=rf, data=data_test$train)
                            })
  
  stopCluster(cl)
  sum_table <- 0
  
  for(i in 1:length(all)){
    sum_table <- sum_table+all[[i]][[1]]
  }
  
  sum_importance <- 0
  
  for(i in 1:length(all)){
    sum_importance <- sum_importance + all[[i]][[2]]$importance
  }
  
  sum_table <- sum_table/matrix(c(sum(sum_table[,1]), sum(sum_table[,1]), sum(sum_table[,2]), sum(sum_table[,2])), nrow = 2) * 100
  
  list(all=sum_table, all_importance=sum_importance, all_train=all)
}




test_bootstrap_variable <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1), apprent=NULL, nb_train=1, proportion_real_data=0.1, type="comptage", base_binaire=Base_binary_matrix(ncol(data1))){
  
  
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
      each <- parLapply(cl, variables, function(x){
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




test_bootstrap <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1), proportion_real_data=0.1, nb_train_all=1, type="comptage", base_binaire=Base_binary_matrix(ncol(data1))){
  
  all <- test_bootstrap_all(data1, apprent = apprent, nb_train = nb_train_all, type=type, base_binaire=base_binaire, proportion_real_data = proportion_real_data)
  each <- test_bootstrap_variable(data1, apprent = apprent, type=type, base_binaire=base_binaire, proportion_real_data = proportion_real_data)
  
  list(all=all, table_each=each$table_each, misclassification=each$misclassification)
}



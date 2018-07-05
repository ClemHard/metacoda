library(randomForest)
library(doParallel)
library(dplyr)

source("Rscript/bootstrap.R")



create_data_frame_test <- function(data1, apprent, test_real_data=NULL, type){
  
  
  boot <- simulation(apprent,type = type)

  data_train <- rbind(data1, boot)
  metadata <- c(rep("real", nrow(data1)), rep("simu", nrow(boot))) %>% as.factor()
  
  nb_sample <- 1000
  boot_test <- simulation(apprent, nb_sample = nb_sample, type = type)
  
  train <- data.frame(data_train, metadata ,row.names = NULL)
   
  
  colnames(train) <- c(paste("X",1:ncol(data_train), sep=""), "metadata")
  
  test <- data.frame(boot_test, row.names = NULL)
  metadata_test <- c(rep("simu", nrow((test))), rep("real", nrow(test_real_data)))
  test <- rbind(test, data.frame(test_real_data, row.names = NULL))
  
  colnames(test) <- paste("X",1:ncol(test), sep="")
  
 
  list(train=train, test=test, metadata_test=metadata_test) 
}



test_bootstrap_all <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1), apprent=NULL, nb_train=1, proportion_real_data=0.1, type="comptage"){
  
  if(is.null(apprent)) apprent <- apprentissage(data1, nb_axe = nb_axe, nb_cluster = nb_cluster)
  
  if(type=="ilr") data1 <- data1 %>% MAP() %>% ilr()
  
  if(type=="MAP") data1 <- data1 %>% MAP()
    
  all <- lapply(1:nb_train, function(x){
    
                            rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*proportion_real_data))
                            test_real_data <- data1[rand_real_data,]
                            data1 <- data1[-rand_real_data,]
                            
                            
                            data_test <- create_data_frame_test(data1, apprent=apprent, test_real_data = test_real_data,type=type)
                            rf <- randomForest(metadata~. , data_test$train)
                            pred <- predict(rf, data_test$test) %>% table(data_test$metadata_test)
                            pred
                            })
  sum_table <- 0
  
  for(i in 1:length(all)){
    sum_table <- sum_table+all[[i]]
  }
  sum_table <- sum_table/matrix(c(sum(sum_table[,1]), sum(sum_table[,1]), sum(sum_table[,2]), sum(sum_table[,2])), nrow = 2) * 100
  
  list(all=sum_table, all_train=all)
}




test_bootstrap_variable <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1), apprent=NULL, nb_train=1, proportion_real_data=0.1, type="comptage"){
  
  
  if(is.null(apprent)) apprent <- apprentissage(data1, nb_axe = nb_axe, nb_cluster = nb_cluster)
  
  if(type=="ilr") data1 <- data1 %>% MAP() %>% center_scale(scale=FALSE) %>% ilr()
  
  if(type=="MAP") data1 <- data1 %>% MAP()
  
  rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*proportion_real_data))
  test_real_data <- data1[rand_real_data,]
  data1 <- data1[-rand_real_data,]
  
  data_test <- create_data_frame_test(data1, apprent=apprent, test_real_data = test_real_data, type=type)
  
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




test_bootstrap <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1), nb_train_all=1, type="comptage"){
  
  rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*0.10))
  
  apprent <- apprentissage(data1, nb_cluster = nb_cluster, nb_axe = nb_axe)
  all <- test_bootstrap_all(data1, apprent = apprent, nb_train = nb_train_all, type=type)
  each <- test_bootstrap_variable(data1, apprent = apprent, type=type)
  
  list(all=all, table_each=each$table_each, misclassification=each$misclassification)
}



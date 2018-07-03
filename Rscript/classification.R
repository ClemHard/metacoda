library(randomForest)
library(parallel)

source("Rscript/bootstrap.R")



test_bootstrap <- function(data1, nb_cluster=NULL, nb_axe=NULL, nb_sample=nrow(data1)){
  
  rand_real_data <- sample(1:nrow(data1), floor(nrow(data1)*0.10))
  
  boot <- bootstrap(data1, PCA=TRUE, nb_cluster = nb_cluster, nb_axe = nb_axe)
  data_train <- rbind(data1, boot$data)
  metadata <- c(rep("real", nrow(data1)), rep("simu", nrow(boot$data))) %>% as.factor()
  
  nb_sample <- 1000
  boot_test <- simulation(boot$apprent, nb_sample = nb_sample)
  
  train <- data.frame(data_train, metadata ,row.names = NULL)
  deleted_data_train <- train[rand_real_data, 1:ncol(data_train)]
  train <- train[-rand_real_data,]
  
  colnames(train) <- c(paste("X",1:ncol(data_train), sep=""), "metadata")
  
  test <- data.frame(boot_test, row.names = NULL)
  test <- rbind(test, deleted_data_train)
  
  colnames(test) <- paste("X",1:ncol(test), sep="")
  
  
  rf <- randomForest(metadata~. , train)
  all <- predict(rf, test) %>% table(c(rep("simu", nb_sample), rep("real", length(rand_real_data))))
  
  
  variables <- colnames(test)
  
  no_cores <- detectCores()-1
  if(no_cores<1) no_cores <- 1
  
  # cl <- makeCluster(no_cores)
  # 
  # each <- parSapply(cl, variables, function(x){
  #   rf <- randomForest(as.formula(paste("metadata~",x)), train)
  #   predict(rf, test) %>% summary()
  # })
  # stopCluster(cl)
  
  each <- lapply(variables, function(x){
    rf <- randomForest(as.formula(paste("metadata~",x)), train)
    predict(rf, test) %>% table(c(rep("simu", nb_sample), rep("real", length(rand_real_data))))
  })

  misclassification_each <- sapply(each, function(x){x[1, 2]})

  names(misclassification_each) <- colnames(data1)
  names(each) <- colnames(data1)
  
  list(all=all, table_each=each, misclassification=misclassification_each)
}



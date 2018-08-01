library(ggplot2)
library(gridExtra)
library(dplyr)
library(magrittr)

source("Rscript/test_groupe.R")

Graph_proportion_evolution<-function(data, abscisse=1:nrow(data)){
  
  data <- norm_data(data)
  time=abscisse
  D <- ncol(data)
  
  somme <- sum(data[1,])
  data <- data/somme
  
  data <- as.data.frame.table(data)
  data$Var1 <- rep(time, D)
  
  
  new.data <- group_by(data, Var2, Freq)
  
  g<- new.data %>% ggplot(aes(x=Var1, y=Freq))+
    geom_line(aes(group=Var2, color=Var2))+
    labs(y="proportions", x="alpha")+
    scale_color_discrete(name=NULL)+
    theme_bw()
  g
}





Graph_cumulative_evolution<-function(data, abscisse=1:nrow(data)){
  
  data <- norm_data(data)
  time=abscisse
  D <- ncol(data)
  
  data <- closure(data)
  data <- t(apply(data,1,cumsum))
  
  for(i in D:2){
    data[,i] <- data[,i]-data[,i-1]
  }
  
  data <- as.data.frame.table(data)
  data$Var1 <- rep(time, D)
  
  new.data <- group_by(data, Var2, Freq)
  
  g<- new.data %>% ggplot(aes(x=Var1, y=Freq, fill=Var2))+
    geom_area()+
    labs(y="proportions", x="alpha")+
    scale_fill_discrete(name=NULL)+
    theme_bw()
  g
  
}



graph_biplot_normale <- function(data, metadata_group, nb_graph=1, title=NULL, legend_title="group", coord_biplot=FALSE, base_binaire=Base_binary_matrix(ncol(data))){
  
  b_data <- 0
  
  if(!coord_biplot){
    data_MAP <- MAP(data)
    b_data <- biplot(data_MAP, base_binaire = base_binaire);
    
  }
  else{
    b_data <-  data
  }
  
  if(ncol(b_data$coord)<nb_graph-1){
    stop("number of graph incorrect")
  }
  
  metadata_group <- metadata_group %>% as.character()
  
  
  m <- data.frame(b_data$coord, group=metadata_group, stringsAsFactors = FALSE)
  
  name_group <- nth(summarise(group_by(m,group)),1)
  
  ellipse_confiance <- list()
  for(i in 1:nb_graph) ellipse_confiance[[i]] <- data.frame(x = numeric(), y = numeric(), group = character(), stringsAsFactors = FALSE)
  
  for(i in name_group){
    for(j in 1:nb_graph){
      if(nrow((filter(m,group==i)[,j:(j+1)]))>5){
        if(((filter(m,group==i)[,j:(j+1)]) %>% as.matrix() %>%ilr_inverse() %>% Raduis_test())$Watson<0.187){
          temp <- as.matrix(filter(m,group==i)[,j:(j+1)])
          el <- intervalle_confiance(temp, alpha=0.05,1,moy=apply(temp,2,mean),var=var(temp))
          
          ellipse_confiance[[j]] <- bind_rows(ellipse_confiance[[j]], 
                                              data.frame(x = el[ , 1], 
                                                         y = el[, 2], 
                                                         group = i, 
                                                         stringsAsFactors = FALSE))
          
        }
      }
    }
  }
  
  # check_type_data <- function(data_list){
  #   data_list[,1] <- data_list[,1] %>% as.character() %>% as.numeric()
  #   data_list[,2] <- data_list[,2] %>% as.character() %>% as.numeric()
  #   data_list[,3] <- data_list[,3] %>% as.character()
  #   
  #   data_list
  # }
  # 
  # for(i in 1:length(ellipse_confiance)){
  #   if(nrow(ellipse_confiance[[i]])!=0){
  #     ellipse_confiance[[i]] %<>% check_type_data()
  #   }
  # }
  
  create_graph <- function(data, ellipse){
    
    graph_bc <- list()
    for(i in 1:nb_graph){
      g <- ggplot(m)+geom_point(aes_string(x=names(m)[i], y=names(m)[i+1],col="group")) +
        labs(title=title, y=paste("comp",i+1), x=paste("comp",i)) +
        scale_color_discrete(name=legend_title) + 
        theme(plot.title = element_text(hjust=0.5))
      if(nrow(ellipse_confiance[[i]])!=0){
        data_ellipse <- ellipse_confiance[[i]]
        g <- g + geom_path(data = data_ellipse, aes(x=x, y=y, group=group, col=group))
      }
      graph_bc[[i]] <- g
    }
    graph_bc
  }
  
  create_graph(m,ellipse_confiance)
  
}
library(phyloseq)
library(dplyr)
library(ape)

tree_binaire_matrice <- function(tree){
  
  nb_otus <- length(tree$tip.label)
  tree1 <- ape::reorder.phylo(tree, "postorder")
  mat <- matrix(0,nrow=tree$Nnode, ncol=nb_otus)
  
  fait <- rep(0, length(tree1$edge[,1])+1)
  
  for(i in 1:length(tree$edge[,1])){
    
    if(fait[tree1$edge[i,1]]==0){
      
      if(tree1$edge[i,2]<=nb_otus){
        mat[tree1$edge[i,1]-nb_otus,tree1$edge[i,2]] <- 1
      }else{
        mat[tree1$edge[i,1]-nb_otus,] <- abs(mat[tree1$edge[i,2]-nb_otus,])
      }
      fait[tree1$edge[i,1]] <- 1
      
    }else{
      
      if(tree1$edge[i,2]<=nb_otus){
        mat[tree1$edge[i,1]-nb_otus,tree1$edge[i,2]] <- -1
      }else{
        mat[tree1$edge[i,1]-nb_otus,] <- -abs(mat[tree1$edge[i,2]-nb_otus,]) + mat[tree1$edge[i,1]-nb_otus,]
      }
    }
  }
  mat
}



random_binary_base <- function(D){
  rtree(D) %>% tree_binaire_matrice()
}


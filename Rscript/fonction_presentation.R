library(MASS)
library(randomForest)

source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")


ternary_diagram3 <- function(){
  u0=0.2;
  v0=0.2;
  A=c(u0+0.5,v0+sqrt(3)/2)
  B=c(u0,v0)
  C=c(u0+1,v0)
  
  data <- ligne(c(0.5, 1.1,3),c(1.1, 1.2, 1))
  
  k <- 1
  coord=(1./k)*(sapply(data[,1],mult,t=A)+sapply(data[,2],mult,t=B)+sapply(data[,3],mult,t=C))
  par(mar=c(rep(0,4)))
  plot(coord[1,],coord[2,], ylim=c(v0-0.03,v0+sqrt(3)/2+0.11), type='l', lwd=2, col='red', axes=FALSE, asp=1, xlab="",ylab="")
  
  data <- list(ligne(c(1.2, 1.1,3),c(1.1, 1.2, 1)), ligne(c(2.8, 1.1,3),c(1.1, 1.2, 1)), ligne(c(9, 1.1,3),c(1.1, 1.2, 1)))
  
  lapply(data, function(x){ 
    par(mar=c(rep(0,4)))
    coord=(1./k)*(sapply(x[,1],mult,t=A)+sapply(x[,2],mult,t=B)+sapply(x[,3],mult,t=C))
    lines(coord[1,], coord[2,], lwd=2, col='red')
  })
  
  
  data <- list(ligne(c(1.1, 0.5,3),c(1.2, 1.1, 1)), ligne(c(1.1, 1.2,3),c(1.2, 1.1, 1)), ligne(c(1.1, 2.8,3),c(1.2, 1.1, 1)), ligne(c(1.1, 9,3),c(1.2, 1.1, 1)))
  
  lapply(data, function(x){ 
    coord=(1./k)*(sapply(x[,1],mult,t=A)+sapply(x[,2],mult,t=B)+sapply(x[,3],mult,t=C))
    lines(coord[1,], coord[2,], lwd=2, col='green')
  })
  
  
  segments(B[1],B[2],A[1],A[2])
  segments(C[1],C[2],A[1],A[2])
  segments(B[1],B[2],C[1],C[2])
  text(A[1], A[2], "A", pos=3)
  text(B[1], B[2], "B", pos=2)
  text(C[1], C[2], "C", pos=4)
  
}


ternary_diagram2 <- function(){
  u0=0.2;
  v0=0.2;
  A=c(u0+0.5,v0+sqrt(3)/2)
  B=c(u0,v0)
  C=c(u0+1,v0)
  
  data <- matrix(c(0.3,0.5,0.2), byrow = TRUE, nrow=1)
  k <- 1
  coord=(1./k)*(sapply(data[,1],mult,t=A)+sapply(data[,2],mult,t=B)+sapply(data[,3],mult,t=C))
  plot(coord[1,],coord[2,],xlim=c(u0-0.2,u0+1.2),ylim=c(v0-0.2,v0+1.1), col='red', axes=FALSE, asp=1, pch=16, xlab="",ylab="")
  
  
  segments(B[1],B[2],A[1],A[2])
  segments(C[1],C[2],A[1],A[2])
  segments(B[1],B[2],C[1],C[2])
  
  segments(coord[1,],coord[2,], coord[1,], v0)
  segments(coord[1,],coord[2,], 0.4, 0.55)
  segments(coord[1,],coord[2,], 0.93, 0.68)
  points(coord[1,], coord[2,], col='red', cex=1.3, pch=16)
  text(A[1], A[2], "A", pos=3)
  text(B[1], B[2], "B", pos=2)
  text(C[1], C[2], "C", pos=4)
  
  text(0.6, 0.33,"x1", cex=0.8)
  text(0.46, 0.47,"x2", cex=0.8)
  text(0.75, 0.53,"x3", cex=0.8)
  text(coord[1,], coord[2,]+0.05, "x",cex=1.5)
}



ternary_diagram1 <- function(data ,colour, style, add_line, colour_line, 
                             add_circle, colour_circle, xlab="",ylab="",
                             ...){
  
  data <- norm_data(data)
  u0=0.2;
  v0=0.2;
  A=c(u0+0.5,v0+sqrt(3)/2)
  B=c(u0,v0)
  C=c(u0+1,v0)
  
  k=sum(data[1,])
  coord=(1./k)*(sapply(data[,1],mult,t=A)+sapply(data[,2],mult,t=B)+sapply(data[,3],mult,t=C))
   plot(coord[1,],coord[2,], ylim=c(v0-0.03,v0+sqrt(3)/2+0.11), col=colour, xlab="", ylab="", pch=style, ...)
  segments(B[1],B[2],A[1],A[2])
  segments(C[1],C[2],A[1],A[2])
  segments(B[1],B[2],C[1],C[2])
  text(A[1], A[2], "A", pos=3)
  text(B[1], B[2], "B", pos=2)
  text(C[1], C[2], "C", pos=4)
  
  coord_circle <- (1./k)*(sapply(add_circle[,1],mult,t=A) + sapply(add_circle[,2],mult,t=B) + sapply(add_circle[,3],mult,t=C))
  lines(coord_circle[1,], coord_circle[2,], col=colour_circle, lwd=2)
  
  coord_line <- (1./k)*(sapply(add_line[,1],mult,t=A) + sapply(add_line[,2],mult,t=B) + sapply(add_line[,3],mult,t=C))
  
  lines(coord_line[1, ], coord_line[2, ], col=colour_line, lwd=2)
}



coord_ternary_diagram <- function(data){
  
  data <- norm_data(data)
  u0=0.2;
  v0=0.2;
  A=c(u0+0.5,v0+sqrt(3)/2)
  B=c(u0,v0)
  C=c(u0+1,v0)
  
  k=sum(data[1,])
  coord=(1./k)*(sapply(data[,1],mult,t=A)+sapply(data[,2],mult,t=B)+sapply(data[,3],mult,t=C))
  coord
}




### melange gaussian 3D
density_gaussienne <- function(x, y, moy, sigma){
  xy <- c(x,y)-moy
  1/(2*pi*det(sigma)^(1/2))*exp(-0.5*xy%*%solve(sigma)%*%xy)
}


density_melange <- function(x, y, pro, moy, sigma){
  if(length(moy)!=length(sigma)){
    stop("length different")
  }
  
  d <- expand.grid(x=x,y=y)
  z <- 0
  
  for(i in 1:length(moy)){
    z <- z + pro[i]*apply(d,1,function(x){density_gaussienne(x[1],x[2],moy[[i]],sigma[[i]])}) %>% matrix(nrow = length(x))
  }
  z
}

library(MASS)
library(randomForest)
library(plotly)

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
  plot(coord[1,],coord[2,],xlim=c(u0-0.2,u0+1.2),ylim=c(v0-0.2,v0+1.1), type='l', lwd=2, col='red', axes=FALSE, asp=1, xlab="",ylab="")
  
  data <- list(ligne(c(1.2, 1.1,3),c(1.1, 1.2, 1)), ligne(c(2.8, 1.1,3),c(1.1, 1.2, 1)), ligne(c(9, 1.1,3),c(1.1, 1.2, 1)))
  
  lapply(data, function(x){ 
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
  plot(coord[1,],coord[2,],xlim=c(u0-0.2,u0+1.2),ylim=c(v0-0.2,v0+1.1), col=colour, xlab="", ylab="", pch=style, ...)
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


c <- matrix(c(c(1,2,7)/10, ilr_inverse(c(-1.2,0)), c(1,8,1, 10/3, 10/3, 10/3)/10), byrow=TRUE, nrow=4)
l <- ligne(1.2,c(0.2,1,3))
y <- seq(0, 2*pi, length=1000)
circl <- cbind(1*cos(y)-1.2, 1*sin(y))

#png(filename = "ternary_diagram.png", res = 150, width = 10, height = 10, units = "in")
ternary_diagram1(c, colour=c(11, 2, 4, 1), style=c(16, 16, 16, 3),add_line = l, colour_line = 15, 
                add_circle = (circl%>% ilr_inverse()), colour_circle = 14, 
                asp = 1, axes = FALSE)
#dev.off()

c_ilr <- ilr(c)
l_ilr <- ilr(l)
plot(c_ilr, col=c(11, 2, 4, 1), pch=c(16, 16, 16, 3), cex=1.5, xlim = c(-2.2,2), ylim=c(-1,1.5), asp=1, xlab = "", ylab="")
lines(l_ilr, col=15, lwd=3)
lines(circl, col=14, lwd=2)


comp <- c(2, 5, 12, 0, 1, 0)

print(comp %>% count_to_proportion() %>% fractions())
print(comp %>% MAP %>% fractions() )


boot <- bootstrap_comptage(ravel, nb_axe=15, nb_cluster = 15, PCA=TRUE)
boot1 <- bootstrap_comptage(ravel, nb_axe=15, nb_cluster = 15, PCA=TRUE)

data <- rbind(ravel, boot$data)
metadata <- c(rep("reel", nrow(ravel)), rep("simu", nrow(boot$data))) %>% as.factor()

data_rf <- data.frame(data,metadata=metadata)
rf <- randomForest(metadata~.,data=data_rf)

boot3 <- lapply(colnames(data), function(x){rf <- randomForest(as.formula(paste("metadata~",x)),data=data_rf)
                                    sum(predict(rf, newdata = boot1$data)=="real")}) %>%unlist()

#hist(ravel[,1], breaks=50, main="données réels", xlab="comptage")
#hist(boot$data[,1], breaks=50, main="données simulées", xlab="comptage")

ggplot(data_rf, aes(x=OTU_163,fill=metadata, color=metadata)) + geom_histogram(position = "dodge")





set.seed(20160919)
x <- rnorm(100)
noise <- matrix(rnorm(200, 0.1), ncol = 2)
mu <- c(-1, 1)
data.perfect <- (x %*% t(c(x = 1, y = 2)) + matrix(rep(mu, each =
                                                         length(x)), ncol = 2)) %>% as.data.frame()
segment.data <- cbind(data.perfect, data.perfect + noise)
names(segment.data) <- c("x", "y", "xend", "yend")
p.template <- ggplot(data = NULL, aes(x, y)) + geom_abline(slope = 2,
                                                           intercept = 3, color = "darkred") + scale_x_continuous(limits = c(-5,
                                                                                                                             3)) + scale_y_continuous(limits = c(-4, 5))
grid.arrange(ggplot(data.frame(x = x, y = 0), aes(x, y)) + geom_point()
             + geom_hline(yintercept = 0, color = "darkred") +
               ggtitle(expression(Latent~Space~(W))),
             p.template + geom_point(data = data.perfect) +
               ggtitle(expression(Parameter~Space~(mu+WB))),
             p.template + geom_point(data = data.perfect, color =
                                       "red") + ggtitle(expression(Observation~Space~(mu+WB))),
             p.template + geom_point(data = data.perfect, color = "red") +
               geom_point(data = data.perfect + noise) +
               geom_segment(data = segment.data, aes(xend = xend, yend
                                                     = yend), color = "black", alpha = 0.2) +
               ggtitle(expression(Observation~Space~(mu+WB+E))),
             ncol = 2)



x <- seq(-5, 6, length=10000)
y1 <- dnorm(x, mean=-1.8)*0.3
y2 <- dnorm(x, mean=1.6, sd=1.4)*0.5
y3 <- dnorm(x, mean=-3, sd=0.4)*0.2

plot(x, y1, type='l', col='red', ylim=c(0,0.3), xlim=c(-4,5), ylab="y")
lines(x, y2, col='blue')
lines(x,y3, col="green")
lines(x, y1+y2+y3)
legend("topright", lty=c(1,1,1), col=c("red", "blue", "green", "black"), legend=c("G1","G2", "G3","G1+G2+G3"))



n <- 100000
r <- rnorm(n*0.4, -1.8)
r1 <- rnorm(n*0.6, 1.6, 1.4)

hist(c(r1, r), proba=TRUE, breaks = 50, main="mélange gaussien", xlab = "x")
lines(x, y1, col="green")
lines(x, y2, col="blue")
lines(x, y1+y2, col='red')

legend("topright", lty=c(1,1,1), col=c("green", "blue", "red"), legend=c("G1","G2","G1+G2"))



g <- graph_biplot_normale(ravel, "data")
g[[1]] + geom_density2d(aes(X1, X2))




#####graph droites -> points

droites <- function(x_ini, y_ini, x_end, y_end){
  coeff <- (y_end-y_ini)/(x_end-x_ini)
  x <- seq(x_ini, x_end, length=1000)
  y <- y_ini + (x-x_ini)*coeff
  data.frame(x=x,y=y)
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


x <- seq(0,10,length=100)
y <- seq(0,10, length=100)
d <- expand.grid(x=x,y=y)

mu <- list(c(3,2.5), c(7,2), c(5, 6))
sig <- list(matrix(c(1.4,0,0,0.5),nrow=2), matrix(c(2,0,0,0.7),nrow=2), matrix(c(1,0,0,3),nrow=2))
prob=c(0.3, 0.2, 0.5)

z <- density_melange(x, y, prob, mu, sig)
plot_ly(x=x, y=y, z=z, type="surface")                                  

# %>% layout(scene=list(camera = list(eye = list(x = 1,y =1.4, z=0),
                                    # center = list(x = 0, y = 0, z = 1)

# t <- biplot(ravel %>% MAP())
# y <- mvrnorm(n=394, rep(0, 3), matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)) %*% t(t$vector[1:10,1:3])
# y <- y + mvrnorm(394, rep(0, 10), diag(10))


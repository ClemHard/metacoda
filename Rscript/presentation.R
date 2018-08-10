library(MASS)
library(randomForest)

source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")


c <- matrix(c(1,2,7,5,2,3, 1,8,1, 10/3, 10/3, 10/3)/10, byrow=TRUE, nrow=4)
l <- ligne(1.2,c(0.2,1,3))
y <- seq(0, 2*pi, length=1000)
circl <- cbind(1*cos(y)-1.2, 1*sin(y))

#png(filename = "ternary_diagram.png", res = 150, width = 10, height = 10, units = "in")
ternary_diagram(c, colour=c(11, 2, 4, 1), style=c(16, 16, 16, 3),add_line = l, colour_line = 15, 
                add_circle = (circl%>% ilr_inverse()), colour_circle = 14, 
                asp = 1, axes = FALSE)
#dev.off()

print(c)


c_ilr <- ilr(c)
l_ilr <- ilr(l)
plot(c_ilr, col=c(11, 2, 4, 1), pch=c(16, 16, 16, 3), cex=1.5, xlim = c(-2.2,2), ylim=c(-1,1.5), asp=1, xlab = "x1", ylab="x2")
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
y1 <- dnorm(x, mean=-1.8)*0.4
y2 <- dnorm(x, mean=1.6, sd=1.4)*0.6

plot(x, y1, type='l', col='red', main="mélange gaussien", ylim=c(0,0.2), xlim=c(-4,5), ylab="y")
lines(x, y2, col='blue')
lines(x, y1+y2)
legend("topright", lty=c(1,1,1), col=c("red", "blue", "black"), legend=c("G1","G2","G1+G2"))


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
density_melange <- function(x, y, moy, sigma){
  xy <- c(x,y)-moy
  1/(2*pi*det(sigma)^(1/2))*exp(-0.5*xy%*%solve(sigma)%*%xy)
}



x <- seq(0,10,length=100)
y <- seq(0,10, length=100)
d <- expand.grid(x=x,y=y)

mu <- c(5,5)
sig <- matrix(c(2,0.1,0.1,2),nrow=2)

z <- apply(d,1,function(x){density_melange(x[1],x[2],mu,sig)}) %>% matrix(nrow = length(x))
plot_ly(x=x,y=y,z=z) %>% add_surface()


t <- biplot(ravel %>% MAP())
y <- mvrnorm(n=394, rep(0, 3), matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)) %*% t(t$vector[1:10,1:3])
y <- y + mvrnorm(394, rep(0, 10), diag(10))


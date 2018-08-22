library(MASS)
library(randomForest)
library(plotly)

source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")
source("Rscript/fonction_presentation.R")




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

n <- 20
data <- matrix(c(runif(n,1,7),runif(n,0.1,1),runif(n,1,7)),nrow=n)
par(mar=c(rep(0,4)))
ternary_diagram(data %>% closure(), colour="red", type_point = 16)
ternary_diagram(data %>% closure() %>% center_scale(), colour="red", type_point = 16)


comp <- c(2, 5, 12, 0, 1, 0)

print(comp %>% count_to_proportion() %>% fractions())
print(comp %>% MAP %>% fractions() )


boot <- bootstrap(ravel, nb_axe=15, nb_cluster = 15)
boot1 <- bootstrap(ravel, nb_axe=15, nb_cluster = 15)

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



n <- 1000000
r1 <- rnorm(n*0.3, -1.8)
r2 <- rnorm(n*0.5, 1.6, 1.4)
r3 <- rnorm(n*0.2, -3, 0.4)

hist(c(r1, r2, r3), proba=TRUE, breaks = 50, main="melange gaussien", xlab = "x")
lines(x, y1, col="red")
lines(x, y2, col='blue')
lines(x,y3, col="green")
lines(x, y1+y2+y3)
legend("topright", lty=c(1,1,1), col=c("red", "blue", "green", "black"), legend=c("G1","G2", "G3","G1+G2+G3"))



g <- graph_biplot_normale(ravel, "data")
g[[1]] + geom_density2d(aes(X1, X2))



x <- seq(0,10,length=100)
y <- seq(0,10, length=100)
d <- expand.grid(x=x,y=y)

mu <- list(c(3,2.5), c(7,2), c(5, 6))
sig <- list(matrix(c(1.4,0,0,0.5),nrow=2), matrix(c(2,0,0,0.7),nrow=2), matrix(c(1,0,0,3),nrow=2))
prob=c(0.3, 0.2, 0.5)

z <- density_melange(x, y, prob, mu, sig)

scene = list(camera = list(eye = list(x = 0.5, y = 1, z = 1)), xaxis=axis, yaxis=axis, zaxis=axis)
plot_ly(x=x, y=y, z=z, type="surface") %>% hide_colorbar()  %>%
  layout(scene = scene)

# %>% layout(scene=list(camera = list(eye = list(x = 1,y =1.4, z=0),
                                    # center = list(x = 0, y = 0, z = 1)

y <- mvrnorm(n=1000, rep(0, 3), matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)) %*% t(t$vector[1:10,1:3])
y <- y + mvrnorm(1000, rep(0, 10), diag(10))



x <- matrix(c(0.1,0.4,0.5,0.2,0.3,0.5,0.4,0.4,0.2,0.5,0.3,0.2),byrow=TRUE, nrow=4)
x1 <- coord_ternary_diagram(x)
cat("euclidienne:",dist(x[1:2,]), "simplexe:", dist_simplex(x[1:2,]))
cat("euclidienne:",dist(x[3:4,]), "simplexe:", dist_simplex(x[3:4,]))
ternary_diagram(x, type=16, colour = 1:4, cex=1.3)
lines(x1[1,1:2],x1[2,1:2])
lines(x1[1,3:4],x1[2,3:4])
x2 <- x %>% ilr()
plot(x2, col=1:4, pch=16, xlim=c(-0.5,1), asp=1, xlab="x", ylab="y", cex=1.3)
lines(x2[1:2,1],x2[1:2,2])
lines(x2[3:4,1],x2[3:4,2])

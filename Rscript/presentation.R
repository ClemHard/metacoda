library(MASS)
library(randomForest)
library(plotly)
library(ellipse)
library(compositions)

source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")
source("Rscript/fonction_presentation.R")



c <- matrix(c(c(3.5,4,3)/10.5, ilr_inverse(c(-1.2,0)), c(1,8,1, 10/3, 10/3, 10/3)/10), byrow=TRUE, nrow=4)

t <- c(2,1,3)
ternary_diagram_vide()
points(c[1,] %>% coord_ternary_diagram(), pch=16, col=4, cex=2)
points(c[2,] %>% coord_ternary_diagram(), pch=16, col='red',cex=2)
points(c[3,] %>% coord_ternary_diagram(), pch=16,cex=2)

l <- ligne(c[3,],t)
l %>% coord_ternary_diagram() %>% lines(col=15, lwd=2)

l <- ligne(c[2,],t)
l %>% coord_ternary_diagram() %>% lines(col=15, lwd=2)

l <- ligne(c[1,],t)
l %>% coord_ternary_diagram() %>% lines(col=15, lwd=2)


t <- c(1,4,3)

l <- ligne(c[3,],t)
l %>% coord_ternary_diagram() %>% lines(col=3, lwd=2)

l <- ligne(c[2,],t)
l %>% coord_ternary_diagram() %>% lines(col=3, lwd=2)

l <- ligne(c[1,],t)
l %>% coord_ternary_diagram() %>% lines(col=3, lwd=2)




c <- matrix(c(c(1,2,7)/10, ilr_inverse(c(-1.2,0)), c(1,8,1, 10/3, 10/3, 10/3)/10), byrow=TRUE, nrow=4)
y <- seq(0, 2*pi, length=1000)
circl <- (cbind(1*cos(y)-1.2, 1*sin(y))) %>% ilr_inverse()

#png(filename = "ternary_diagram.png", res = 150, width = 10, height = 10, units = "in")
ternary_diagram1(c, colour=c(11, 2, 4, 1), style=c(16, 16, 16, 3),add_line = l, colour_line = 15, 
                add_circle = (circl), colour_circle = 14, 
                asp = 1, axes = FALSE)
#dev.off()

c_ilr <- ilr(c)
l_ilr <- ilr(l)
circl_ilr <- circl %>% ilr
plot(c_ilr, col=c(11, 2, 4, 1), pch=c(16, 16, 16, 3), cex=1.5, xlim = c(-2.2,2), ylim=c(-1,1.5), asp=1, xlab = "", ylab="")
lines(l_ilr, col=15, lwd=3)
lines(circl_ilr, col=14, lwd=2)

n <- 20
data <- matrix(c(runif(n,1,7),runif(n,0.1,1),runif(n,1,7)),nrow=n)
par(mar=c(rep(0,4)))
ternary_diagram(data %>% closure(), colour="red", type_point = 16)
ternary_diagram(data %>% closure() %>% center_scale(), colour="red", type_point = 16)


comp <- c(2, 5, 12, 0, 1, 0)

print(comp %>% count_to_proportion() %>% fractions())
print(comp %>% MAP %>% fractions() )


boot <- bootstrap(liver_500, nb_axe=12, nb_cluster = 12, zero_inflated = FALSE)
data <- rbind(liver_500, boot$data)
metadata <- c(rep("reel", nrow(liver_500)), rep("simu", nrow(boot$data))) %>% as.factor()

data_rf <- data.frame(data,metadata=metadata)
rf <- randomForest(metadata~.,data=data_rf)

#hist(ravel[,1], breaks=50, main="données réels", xlab="comptage")
#hist(boot$data[,1], breaks=50, main="données simulées", xlab="comptage")

ggplot(data_rf, aes(x=CAG00214_hs_9.9,fill=metadata, color=metadata)) + geom_histogram(position = "dodge") + xlim(c(-1, 50))





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

y <- matrix(0, ncol=10, nrow=1000)
y[,1:3] <- mvrnorm(n=1000, rep(0, 3), matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)) 
y <- y + mvrnorm(1000, rep(0, 10), seq(2,0.1,length=10)^7*diag(10))
z <- cov(y) %>% eigen()
RMSE <- sapply(1:D, function(x){mean((new_data$coord[,x]^2))}) %>% sort() %>% cumsum() %>% sort(decreasing = TRUE) 
RMSE <- RMSE/sum(RMSE)
plot(z$values, ylab="y", xlab="Index", pch=16)
abline(v=4, col='red')



x <- matrix(c(0.1,0.4,0.5,0.2,0.3,0.5,0.4,0.4,0.2,0.5,0.3,0.2),byrow=TRUE, nrow=4)
x1 <- coord_ternary_diagram(x)
cat("euclidienne:",dist(x[1:2,]), "simplexe:", dist_simplex(x[1:2,]))
cat("euclidienne:",dist(x[3:4,]), "simplexe:", dist_simplex(x[3:4,]))
ternary_diagram(x, type=16, colour = 1:4, cex=1.3)
lines(x1[1:2,1],x1[1:2,2])
lines(x1[3:4,1],x1[3:4,2])
x2 <- x %>% ilr()
plot(x2, col=1:4, pch=16, xlim=c(-0.5,1), asp=1, xlab="x", ylab="y", cex=1.3)
lines(x2[1:2,1],x2[1:2,2])
lines(x2[3:4,1],x2[3:4,2])



ravel_boot <- bootstrap_presentation(liver_500, nb_cluster = 12, nb_axe = 12, type="comptage")
data <- data.frame(data=c(ravel_boot$data %>% apply(1, sum), liver %>% apply(1, sum)), metadata=c(rep("simu",nrow(liver_500)),rep("real", nrow(liver_500))))
ggplot(data, aes(x=data, color=metadata, fill=metadata)) + geom_histogram(position="dodge", bins=30)


#####courbe loi normale
mat1 <- matrix(c(0.3,0,0,0.3), nrow=2) *0.1
mat2 <- matrix(c(0.4,-0.5,-0.5,2), nrow=2)*0.1
mat3 <- matrix(c(0.8,0.9,0.9,1.5), nrow=2)*0.1
u1 <- c(0,0)*0.5
u2 <- c(3,-2)*0.5
u3 <- c(-1.8,-0.5)*0.5
col1="red"
col2="green"
col3="gray"

ellipse(mat2, centre = u2) %>% plot(type="l", asp=1, xlim=c(-2,2), ylim=c(-2,1), col=col2)
ellipse(mat2, centre = u2, level=0.8) %>% lines(col=col2)
ellipse(mat2, centre = u2, level=0.4) %>% lines(col=col2)

ellipse(mat1, centre =u1) %>% lines(col=col1)
ellipse(mat1, centre =u1, level=0.8) %>% lines(col=col1)
ellipse(mat1, centre =u1, level=0.4) %>% lines(col=col1)

ellipse(mat3, centre =u3) %>% lines(col=col3)
ellipse(mat3, centre =u3, level=0.8) %>% lines(col=col3)
ellipse(mat3, centre =u3, level=0.4) %>% lines(col=col3)


ternary_diagram_vide()
ellipse(mat2, centre = u2) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col2)
ellipse(mat2, centre = u2, level=0.8) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col2)
ellipse(mat2, centre = u2, level=0.4) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col2)


ellipse(mat1, centre = u1) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col1)
ellipse(mat1, centre = u1, level=0.8) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col1)
ellipse(mat1, centre = u1, level=0.4) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col1)


ellipse(mat3, centre = u3) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col3)
ellipse(mat3, centre = u3, level=0.8) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col3)
ellipse(mat3, centre = u3, level=0.4) %>% ilr_inverse %>% coord_ternary_diagram() %>% lines(col=col3)




#####
x <- runif(3, 0.1,1.2) %>% round(digits = 2) %>% norm_data()
x1 <- runif(3, 0.1,1.2) %>% round(digits = 2) %>% norm_data()

perturbation(x,x1) %>% ilr
(x %>% ilr) + (x1%>% ilr)

power(x,2.3) %>% ilr
(x %>% ilr) * 2.3

data <- 
ternary_diagram_vide()
Hongite %>%

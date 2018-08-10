source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")
source("Rscript/test_bootstrap.R")
source("Rscript/tree_phyloseq.R")
source("Rscript/apprentissage_supervise.R")


chaillou_p <- count_to_proportion(chaillou)
mach_500_p <- count_to_proportion(mach)
ravel_p <- count_to_proportion(ravel)
vacher_p <- count_to_proportion(vacher)


chaillou_MAP <- MAP(chaillou)
mach_500_MAP <- MAP(mach_500)
ravel_MAP <- MAP(ravel)
vacher_MAP <- MAP(vacher)


### reduction de dimension
bc <- biplot(chaillou_MAP)
bm <- biplot(mach_500_MAP)
br <- biplot(ravel_MAP)
bv <- biplot(vacher_MAP)


#chaillou
grid.arrange(grobs=graph_biplot_normale(chaillou, metadata_chaillou$EnvType, 4, "Chaillou", "EnvType"), ncol=2)


#ravel
grid.arrange(grobs=graph_biplot_normale(ravel, metadata_ravel$CST, 4, "Ravel", "CST"), ncol=2)


#mach
grid.arrange(grobs=graph_biplot_normale(mach_500, metadata_mach$Weaned, 4, "Mach", "Weaned"), ncol=2)


#vacher
grid.arrange(grobs=graph_biplot_normale(vacher, metadata_vacher$tree, 4, "vacher", "pmInfection"), ncol=2)


#liver
grid.arrange(grobs=graph_biplot_normale(liver_500, metadata_liver$status, 4, "liver", "status"), ncol=2)


### test de normalite

##test global

# chaillou
marginal_univariate_distributions(chaillou_MAP)
Bivariate_angle_distribution(chaillou_MAP)
Raduis_test(chaillou_MAP)

# mach
marginal_univariate_distributions(mach_500_MAP)
Bivariate_angle_distribution(mach_500_MAP)
Raduis_test(mach_500_MAP)

# ravel
marginal_univariate_distributions(ravel_MAP)
Bivariate_angle_distribution(ravel_MAP)
Raduis_test(ravel_MAP)



### test normalite des donnees en groupe

#ravel
val <- c("I", "II", "III", "IV", "V")

sapply(val,function(x){
  (Raduis_test(ravel_MAP[which(metadata_ravel$CST==x),]))
})

#mach
val <- c("TRUE", "FALSE")

sapply(val,function(x){
  (Raduis_test(mach_500_MAP[which(metadata_mach$Weaned==x),]))
})

#chaillou
val <- c("BoeufHache", "Crevette", "DesLardons", "FiletCabillaud", "FiletSaumon", "SaucisseVolaille", "SaumonFume", "VeauHache")

sapply(val,function(x){
  (Raduis_test(chaillou_MAP[which(metadata_chaillou$EnvType==x),]))
})



### test normalite donnees biplot


#chaillou
val <- c("BoeufHache", "Crevette", "DesLardons", "FiletCabillaud", "FiletSaumon", "SaucisseVolaille", "SaumonFume", "VeauHache")

o <- sapply(val,function(x){
  (Raduis_test(ilr_inverse(bc$coord[which(metadata_chaillou$EnvType==x), 1:2])))
})


#ravel
val <- c("I", "II", "III", "IV", "V")

sapply(val,function(x){
  (Raduis_test(ilr_inverse(br$coord[which(metadata_ravel$CST==x), 1:2])))
})


#mach 500
val <- c("TRUE", "FALSE")

sapply(val,function(x){
  (Raduis_test(ilr_inverse(bm$coord[which(metadata_mach$Weaned==x), 1:2])))
})
############### simulation


#### dans R
D <- 5

#parametres
# Sigma1 <- matrix(0, ncol=D, nrow=D)

# val <- runif(0.5*D*(D-1), -2, 2)
# Sigma1[lower.tri(Sigma1)] <- val
# diag(Sigma1) <- runif(D, 0, 2)
# Sigma1 <- Sigma1%*%t(Sigma1)

Sigma1 <- diag(1, D, D)

mu1 <- runif(D, -1, 1)
mu2 <- mu1+((runif(D,-1,1)>0)*2-1)*0.1 # perturbation de la moyenne de +-0.1


#data
N=round(10^(seq(1.5,3,length=30)))

total <- 0
NB <- 600
for(i in 1:NB){
  data1 <- mvrnorm(N[length(N)], mu1, Sigma1) %>%ilr_inverse()
  data2 <- mvrnorm(N[length(N)], mu2, Sigma1) %>%ilr_inverse()
  
  total <- total+sapply(N, function(x){
    
    testing(data1[1:x,], data2[1:x,], 0.05, 3)$result
    
  })/NB
}

m <- data.frame(N=N, power=1-total)
g <- ggplot(m, aes(N, power))+ geom_line()+labs(title="Evolution de la puissance")+theme(panel.background = element_rect(fill="white"),plot.title = element_text(hjust=0.5), axis.line = element_line(colour="black"))
g



#### dans N
D <- 50
nb_sample <- 200

prob1 <- runif(D, 1, 2)
prob2 <- prob1+((runif(D,-1,1)>0)*2-1)*0.05 # perturbation de +-0.1

prob1 <- prob1/sum(prob1)
prob2 <- prob2/sum(prob2)


N=round(10^(seq(1,4,length=3)))

total <- 0
NB <- 500
for(i in 1:NB){
  
  total <- total+sapply(N, function(x){
    
    data1 <- t(rmultinom(nb_sample, size=x, prob1)) %>% MAP()
    data2 <- t(rmultinom(nb_sample, size=x, prob2)) %>% MAP()
    
    
    testing(data1, data2, 0.05, 3)$result
    
  })/NB
}

m <- data.frame(N=N, power=1-total)
g <- ggplot(m, aes(N, power))+ geom_line()+labs(title="Evolution de la puissance", x="N multinomiale")+theme(panel.background = element_rect(fill="white"),plot.title = element_text(hjust=0.5), axis.line = element_line(colour="black"))
g



###test variation D
nb_sample <- 200
D <- round(2^(seq(2,3,length=3)))
N=500

total <- 0
NB <- 500
for(i in 1:NB){
  
  total <- total+sapply(D, function(x){
    
    
    prob1 <- runif(x, 1, 2)
    
    prob1 <- prob1/sum(prob1)
    
    data1 <- t(rmultinom(nb_sample, size=N, prob1)) %>% MAP()
    
    
    mean(eigen(var(ilr(data1)))$values)
    
  })/NB
}

m <- data.frame(D=D, power=total)
g <- ggplot(m, aes(D, power))+ geom_line()+labs(title="Evolution de la puissance", x="N multinomiale")+theme(panel.background = element_rect(fill="white"),plot.title = element_text(hjust=0.5), axis.line = element_line(colour="black"))
g













## idee

abondance <- function(data, i){
  sum(data<=i)
}

x <- 0:8
y <- sapply(x, abondance, data=chaillou)
y1 <- sapply(x, abondance, data=(chaillou_boot$data))
plot(x, y, type='l', ylim=c(min(y, y1), max(y, y1)))
lines(x, y1 ,col='red')


u <- t_chaillou$all_importance
r <- order(u)
plotdata <- reshape2::melt(chaillou, value.name = "count", 
                           varnames = c("sample", "OTU")) %>% 
  mutate(Status = "real") %>% 
  bind_rows(reshape2::melt(chaillou_boot$data, value.name = "count", 
                           varnames = c("sample", "OTU")) %>% 
              mutate(Status = "simu")) %>% 
  group_by(OTU, count, Status) %>% 
  summarize(n_sample = n())

plotdata <- plotdata %>% 
  inner_join(tibble(OTU = colnames(chaillou), 
                    importance = u[, 1], 
                    VIP = rank(-u[ , 1])), 
             by = "OTU")
    

ggplot(data = plotdata %>% filter(VIP >= 500), 
       aes(x = count, y = n_sample, fill = Status)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~OTU, scales = "free")

for(i in 490:508){hist((chaillou %>% MAP() %>% ilr())[,r[i]], main=i, breaks=500)}
for(i in 1:10){hist((chaillou %>% MAP() %>% ilr())[,r[i]], main=i, breaks=500)}



x <- 0:8
y <- sapply(x, abondance, data=mach_500)
y1 <- sapply(x, abondance, data=(mach_boot$data))
plot(x, y, type='l', ylim=c(min(y, y1), max(y,y1)))
lines(x, y1 ,col='red')




u <- t_ravel$all_importance
r <- order(u)
plotdata <- reshape2::melt(ravel, value.name = "count", 
                           varnames = c("sample", "OTU")) %>% 
  mutate(Status = "real") %>% 
  bind_rows(reshape2::melt(ravel_boot$data, value.name = "count", 
                           varnames = c("sample", "OTU")) %>% 
              mutate(Status = "simu")) %>% 
  group_by(OTU, count, Status) %>% 
  summarize(n_sample = n())

plotdata <- plotdata %>% 
  inner_join(tibble(OTU = colnames(ravel), 
                    importance = u[, 1], 
                    VIP = rank(-u[ , 1])), 
             by = "OTU")


ggplot(data = plotdata %>% filter(VIP >=240), 
       aes(x = count, y = n_sample, fill = Status)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~OTU, scales = "free")  


u <-  t_liver_500$all_importance
r <- order(u)
for(i in 485:500){hist(liver_500[,r[i]], main=i, breaks=5000)}
for(i in 1:15){hist(liver_500[,r[i]], main=i, breaks=500, xlim=c(0,30))}




x <- 0:8
y <- sapply(x, abondance, data=vacher)
y1 <- sapply(x, abondance, data=(vacher_boot$data))
plot(x, y, type='l', ylim=c(min(y, y1), max(y,y1)))
lines(x, y1 ,col='red')



u <- t_vacher$all_importance
r <- order(u)
for(i in 100:114){hist(vacher[,r[i]], main=i, breaks=500, xlim=c(0,30))}
for(i in 1:15){hist(vacher[,r[i]], main=i, breaks=500, xlim=c(0,30))}



x <- 0:8
y <- sapply(x, abondance, data=ravel)
y1 <- sapply(x, abondance, data=(ravel_boot$data))
plot(x, y, type='l', ylim=c(min(y, y1), max(y,y1)))
lines(x, y1 ,col='red')


u <- t_ravel$all_importance
r <- order(u)
for(i in 230:247){hist(ravel[,r[i]], main=i, breaks=500, xlim=c(0,30))}
for(i in 1:15){hist(ravel[,r[i]], main=i, breaks=500, xlim=c(0,30))}




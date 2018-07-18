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
grid.arrange(grobs=graph_biplot_normale(vacher, (metadata_vacher$pmInfection>0)*1, 4, "vacher", "pmInfection"), ncol=2)


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









########clustering



##ravel
k_ravel <- comparaison_k_means(ravel, metadata_ravel$CST, 5, 1)
grid.arrange(grobs=k_ravel$graphics, ncol=1)


hclust_ravel <- comparaison_hclust(ravel, metadata_ravel$CST, 5, 1)
grid.arrange(grobs=hclust_ravel$graphics, ncol=1)

Mclust_ravel <- comparaison_Mclust(ravel, metadata_ravel$CST, 5, 1)
grid.arrange(grobs=Mclust_ravel$graphics, ncol=1)


##mach 500
k_mach_500 <- comparaison_k_means(mach_500, metadata_mach$Weaned, 2, 1)
grid.arrange(grobs=k_mach_500$graphics, ncol=1)

hclust_mach_500 <- comparaison_hclust(mach_500, metadata_mach$Weaned, 2, 1)
grid.arrange(grobs=hclust_mach_500$graphics, ncol=1)

Mclust_mach_500 <- comparaison_Mclust(mach_500, metadata_mach$Weaned, 2, 1)
grid.arrange(grobs=Mclust_mach_500$graphics, ncol=1)



##chaillou
k_chaillou <- comparaison_k_means(chaillou, metadata_chaillou$EnvType, 8, 1)
grid.arrange(grobs=k_chaillou$graphics, ncol=1)

hclust_chaillou <- comparaison_hclust(chaillou, metadata_chaillou$EnvType, 8, 1)
grid.arrange(grobs=hclust_chaillou$graphics, ncol=1) 

Mclust_chaillou <- comparaison_Mclust(chaillou, metadata_chaillou$EnvType, 8, 1)
grid.arrange(grobs=Mclust_chaillou$graphics, ncol=1)



##mach
k_mach <- comparaison_k_means(mach, metadata_mach$Weaned, 2, 1, 1)
grid.arrange(grobs=k_mach$graphics, ncol=1)

hclust_mach <- comparaison_hclust(mach, metadata_mach$Weaned, 2, 1)
grid.arrange(grobs=hclust_mach$graphics, ncol=1)

Mclust_mach <- comparaison_Mclust(mach, metadata_mach$Weaned, 2, 1)
grid.arrange(grobs=Mclust_mach$graphics, ncol=1)



##vacher
k_vacher <- comparaison_k_means(vacher, (metadata_vacher$pmInfection>0)*1, 2, 1)
grid.arrange(grobs=k_vacher$graphics, ncol=1)

hclust_vacher <- comparaison_hclust(vacher, (metadata_vacher$pmInfection>0)*1, 2, 1)
grid.arrange(grobs=hclust_vacher$graphics, ncol=1)

Mclust_vacher <- comparaison_Mclust(vacher, (metadata_vacher$pmInfection>0)*1, 2, 1)
grid.arrange(grobs=Mclust_vacher$graphics, ncol=1)



##liver
k_liver <- comparaison_k_means(liver_500, metadata_liver$status, 2, 1)
grid.arrange(grobs=k_liver$graphics, ncol=1)

hclust_liver <- comparaison_hclust(liver_500, metadata_liver$status, 2, 1)
grid.arrange(grobs=hclust_liver$graphics, ncol=1)

Mclust_liver <- comparaison_Mclust(liver_500, metadata_liver$status, 2, 1)
grid.arrange(grobs=Mclust_liver$graphics, ncol=1)





### bootstrap

#chaillou

chaillou_boot <- bootstrap(chaillou, nb_cluster = 15, nb_axe = 16, type="ilr")
data <- rbind(chaillou %>% MAP() %>% center_scale(scale = FALSE) %>% ilr(), chaillou_boot$data)
#metadata <- c(as.character(metadata_chaillou$EnvType), as.character(chaillou_boot$metadata)) %>% as.factor()
metadata <- c(as.character(rep("real", nrow(chaillou))), rep("simu", nrow(chaillou_boot$data))) %>% as.factor()

grid.arrange(grobs=graph_biplot_normale(data, metadata, 4, "chaillou", "data"), ncol=2)



#ravel

ravel_boot <- bootstrap(ravel, nb_cluster = 4, nb_axe = 4, type="comptage")
data <- rbind(ravel, ravel_boot$data)
# metadata <- c(as.character(metadata_ravel$CST), as.character(ravel_boot$metadata)) %>% as.factor()
metadata <- c(rep("real", nrow(ravel)), rep("simu", nrow(ravel_boot$data))) %>% as.factor()
grid.arrange(grobs=graph_biplot_normale(data, metadata, 4, "ravel", "data"), ncol=2)



#mach500

mach_boot <- bootstrap(mach_500, type="comptage", nb_cluster = 4, nb_axe = 10)
data <- rbind(mach_500, mach_boot$data)
#metadata <- c(as.character(metadata_mach$Weaned), as.character(mach_boot$metadata)) %>% as.factor()
metadata <- c(rep("real", nrow(mach_500)), rep("simu", nrow(mach_boot$data))) %>% as.factor()
grid.arrange(grobs=graph_biplot_normale(data, metadata, 4, "mach 500", "data"), ncol=2)




#vacher

vacher_boot <- bootstrap(vacher, nb_cluster = 4, nb_axe = 11, type="comptage")
data <- rbind(vacher, vacher_boot$data)
# metadata <- c(as.character(metadata_ravel$CST), as.character(ravel_boot$metadata)) %>% as.factor()
metadata <- c(rep("real", nrow(vacher)), rep("simu", nrow(vacher_boot$data))) %>% as.factor()
grid.arrange(grobs=graph_biplot_normale(data, metadata, 4, "vacher", "data"), ncol=2)


#liver
liver_boot <- bootstrap(liver_500, nb_cluster = 8, nb_axe=7, nb_sample = 1000)
data <- rbind(liver_500, liver_boot$data)
# metadata <- c(as.character(metadata_ravel$CST), as.character(ravel_boot$metadata)) %>% as.factor()
metadata <- c(rep("real", nrow(liver)), rep("simu", nrow(liver_boot$data))) %>% as.factor()
grid.arrange(grobs=graph_biplot_normale(data, metadata, 4, "liver", "data"), ncol=2)







### classification

#chaillou
u <- apply(chaillou, 2, function(x){sum(x>10)}) %>% order(decreasing = TRUE)
c <- chaillou[,u]
set.seed(1)
t_chaillou <- test_bootstrap_all(chaillou, nb_cluster = 15, nb_axe = 16, type = "comptage", nb_train=1)
t_chaillou$all


abondance_otus <- apply(chaillou %>% MAP(), 2 , function(x){sum(x>1e-4)})
r <- order(abondance_otus, decreasing = TRUE)
plot(t_chaillou$misclassification[r])




#mach_500
u <- apply(mach_500, 2, function(x){sum(x>10)}) %>% order(decreasing = TRUE)
c <- mach_500[,u]

t_mach_500 <- test_bootstrap_all(mach_500, nb_cluster = 4, nb_axe = 10, nb_train = 1, type="comptage")
t_mach_500$all

abondance_otus <- apply(mach_500 %>% MAP(), 2 , function(x){sum(x>1e-4)})
r <- order(abondance_otus, decreasing = TRUE)
plot(t_mach_500$misclassification[r])



#vacher
t_vacher <- test_bootstrap_all(vacher, nb_cluster = 4, nb_axe = 11, nb_train=1, type="comptage")
t_vacher$all


abondance_otus <- apply(vacher %>% MAP(), 2 , function(x){sum(x>1e-4)})
r <- order(abondance_otus, decreasing = TRUE)
plot(t_vacher$misclassification[r])




#liver
t_liver <- test_bootstrap_all(liver, nb_cluster = 7, nb_axe = 8, nb_train = 1)
t_liver$all


abondance_otus <- apply(liver %>% MAP(), 2 , function(x){sum(x>1e-4)})
r <- order(abondance_otus, decreasing = TRUE)
plot(t_liver$misclassification[r])



#liver 500
t_liver_500 <- test_bootstrap_all(liver_500, nb_cluster = 7, nb_axe = 8, nb_train=1)
t_liver_500$all



abondance_otus <- apply(liver_500 %>% MAP(), 2 , function(x){sum(x>1e-4)})
r <- order(abondance_otus, decreasing = TRUE)
plot(t_liver_500$misclassification[r])




#ravel
t_ravel <- test_bootstrap_all(ravel, nb_cluster = 4, nb_axe = 11, nb_train = 1, type="comptage")
t_ravel$all



##### bootstrap supervise
c_super <- test_bootstrap_supervise(chaillou, metadata_chaillou$EnvType)
c_super
r_super <- test_bootstrap_supervise(ravel, metadata_ravel$CST)
r_super
l_super <- test_bootstrap_supervise(liver_500, metadata_liver$status)
l_super
m_super <- test_bootstrap_supervise(mach_500, metadata_mach$Weaned)
m_super


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



#################################################
#################################################

##apprentissage supervise



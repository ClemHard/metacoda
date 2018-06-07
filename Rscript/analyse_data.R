source("Rscript/read_metagenomic_data.R")
source("Rscript/script1.R")

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
grid.arrange(grobs=graph_biplot_normale(chaillou, metadata_chaillou$EnvType,4, "Chaillou", "EnvType"), ncol=2)

#ravel
grid.arrange(grobs=graph_biplot_normale(ravel, metadata_ravel$CST, 4, "Ravel", "CST"), ncol=2)

#mach
grid.arrange(grobs=graph_biplot_normale(mach_500, metadata_mach$Weaned, 4, "Mach", "Weaned"), ncol=2)

#vacher


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
k_ravel <- comparaison_k_means(ravel, metadata_ravel$CST, 5, 3)
grid.arrange(grobs=k_ravel$graphics, ncol=3)




##mach 500
k_mach_500 <- comparaison_k_means(mach_500, metadata_mach$Weaned, 2, 3)
grid.arrange(grobs=k_mach_500$graphics, ncol=3)




##chaillou
k_chaillou <- comparaison_k_means(chaillou, metadata_chaillou$EnvType, 8, 3)
grid.arrange(grobs=k_chaillou$graphics, ncol=3)


##mach
k_mach <- comparaison_k_means(mach, metadata_mach$Weaned, 2, 3, 1)
grid.arrange(grobs=k_mach$graphics, ncol=3)



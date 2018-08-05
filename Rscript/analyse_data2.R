source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")
source("Rscript/test_bootstrap.R")
source("Rscript/tree_phyloseq.R")
source("Rscript/apprentissage_supervise.R")




### classification
##donne comptage

#chaillou
t_chaillou <- test_bootstrap_all(chaillou, nb_cluster = 15, nb_axe = 16, nb_train=50)
t_chaillou$all


#mach_500
t_mach_500 <- test_bootstrap_all(mach_500, nb_cluster = 4, nb_axe = 10, nb_train=50)
t_mach_500$all


#vacher
t_vacher <- test_bootstrap_all(vacher, nb_cluster = 4, nb_axe = 11, nb_train=50)
t_vacher$all


#liver 500
t_liver_500 <- test_bootstrap_all(liver_500, nb_cluster = 13, nb_axe = 17, nb_train=50)
t_liver_500$all



#ravel
t_ravel <- test_bootstrap_all(ravel, nb_cluster = 12, nb_axe = 12, nb_train=50)
t_ravel$all




##donne ilr
#chaillou
t_chaillou_ilr <- test_bootstrap_all(chaillou, nb_cluster = 15, nb_axe = 16, nb_train=5, type="ilr")
t_chaillou_ilr$all


#mach_500
t_mach_500_ilr <- test_bootstrap_all(mach_500, nb_cluster = 4, nb_axe = 10, nb_train=5, type="ilr")
t_mach_500_ilr$all


#vacher
t_vacher_ilr <- test_bootstrap_all(vacher, nb_cluster = 4, nb_axe = 11,nb_train=5, type="ilr")
t_vacher_ilr$all


#liver 500
t_liver_500_ilr <- test_bootstrap_all(liver_500, nb_cluster = 13, nb_axe = 17, nb_train=5, type="ilr")
t_liver_500_ilr$all



#ravel
t_ravel_ilr <- test_bootstrap_all(ravel, nb_cluster = 12, nb_axe = 12, nb_train=5, type="ilr")
t_ravel_ilr$all


##### bootstrap supervise

## donne comptage
c_super <- test_bootstrap_supervise(chaillou, metadata_chaillou$EnvType, type="comptage", nb_train = 5)
c_super$all
r_super <- test_bootstrap_supervise(ravel, metadata_ravel$CST, type="comptage", nb_train = 5)
r_super$all
l_super <- test_bootstrap_supervise(liver_500, metadata_liver$status, type="comptage", nb_train = 5)
l_super$all
m_super <- test_bootstrap_supervise(mach_500, metadata_mach$Weaned, type="comptage", nb_train = 5)
m_super$all


## donne ilr
c_super_ilr <- test_bootstrap_supervise(chaillou, metadata_chaillou$EnvType, type="ilr", nb_train = 5)
c_super_ilr$all
r_super_ilr <- test_bootstrap_supervise(ravel, metadata_ravel$CST, type="ilr", nb_train = 5)
r_super_ilr$all
l_super_ilr <- test_bootstrap_supervise(liver_500, metadata_liver$status, type="ilr", nb_train = 5)
l_super_ilr$all
m_super_ilr <- test_bootstrap_supervise(mach_500, metadata_mach$Weaned, type="ilr", nb_train = 5)
m_super_ilr$all


### validation croise

#donne de comptage
v_chaillou <- validation_croise(chaillou, metadata_chaillou$EnvType)
v_chaillou
v_ravel <- validation_croise(ravel, metadata_ravel$CST)
v_ravel
v_liver_500 <- validation_croise(liver_500, metadata_liver$status)
v_liver_500
v_mach_500 <- validation_croise(mach_500, metadata_mach$Weaned)
v_mach_500



#donne ilr
v_chaillou_ilr <- validation_croise(chaillou %>% MAP() %>% ilr(), metadata_chaillou$EnvType)
v_chaillou_ilr
v_ravel_ilr <- validation_croise(ravel %>% MAP() %>% ilr(), metadata_ravel$CST)
v_ravel_ilr
v_liver_500_ilr <- validation_croise(liver_500 %>% MAP() %>% ilr(), metadata_liver$status)
v_liver_500_ilr
v_mach_500_ilr <- validation_croise(mach_500 %>% MAP() %>% ilr(), metadata_mach$Weaned)
v_mach_500_ilr


#### classificateur supervise apres classificateur real/simu

# donne comptage
c_chaillou <- classificateur_group_real_simu(chaillou, metadata_chaillou$EnvType, nb_train = 5)
c_chaillou

c_ravel <- classificateur_group_real_simu(ravel, metadata_ravel$CST, nb_train = 5)
c_ravel

c_liver <- classificateur_group_real_simu(liver_500, metadata_liver$status, nb_train = 5)
c_liver

c_mach <- classificateur_group_real_simu(mach_500, metadata_mach$Weaned, nb_train = 5)
c_mach

# donne ilr
c_chaillou_ilr <- classificateur_group_real_simu(chaillou, metadata_chaillou$EnvType, type="ilr", nb_train = 5)
c_chaillou_ilr

c_ravel_ilr <- classificateur_group_real_simu(ravel, metadata_ravel$CST, type="ilr", nb_train = 5)
c_ravel_ilr

c_liver_ilr <- classificateur_group_real_simu(liver_500, metadata_liver$status, type="ilr", nb_train = 5)
c_liver_ilr

c_mach_ilr <- classificateur_group_real_simu(mach_500, metadata_mach$Weaned, type="ilr", nb_train = 5)
c_mach_ilr



# save(t_chaillou, t_mach_500, t_vacher, t_ravel, t_liver_500,
#      t_chaillou_ilr, t_mach_500_ilr, t_vacher_ilr, t_ravel_ilr, t_liver_500_ilr,
#      c_super, r_super, l_super, m_super,
#      c_super_ilr, r_super_ilr, m_super_ilr, l_super_ilr,
#      v_chaillou, v_ravel, v_liver_500, v_mach_500,
#      c_chaillou, c_ravel, c_liver, c_mach,
#      c_chaillou_ilr, c_ravel_ilr, c_liver_ilr, c_mach_ilr,
#      file="results.RData"
#      )


# base binaire aléatoire
plot_aleatoire <- function(l){
  misclass_simu_forest <- sapply(l, function(x){x$random_forest[1,2]})
  misclass_real_forest <- sapply(l, function(x){x$random_forest[2,1]})
  misclass_simu_kNN <- sapply(l, function(x){x$kNN[1,2]})
  misclass_real_kNN <- sapply(l, function(x){x$kNN[2,1]})
  
  data <- data.frame(x=1:length(misclass_real_forest), misclass_real_forest=misclass_real_forest, misclass_simu_forest=misclass_simu_forest, misclass_real_kNN=misclass_real_kNN, misclass_simu_kNN=misclass_simu_kNN)
  
  ggplot(data) + geom_line(aes(x, misclass_real_forest), linetype=2) +
    geom_line(aes(x, misclass_simu_forest)) +
    geom_line(aes(x, misclass_real_kNN), colour=2, linetype=2) +
    geom_line(aes(x, misclass_simu_kNN), colour=2) + labs(y="misclass", x= "index")
}


aleatoire_vacher <- lapply(1:20, function(x){
                r <- random_binary_base(ncol(vacher))
                t_vacher <- test_bootstrap_all(vacher, nb_cluster = 4, nb_axe = 11, nb_train=20, type="comptage", base_binaire = r)
                t_vacher$all
              })

aleatoire_chaillou <- lapply(1:20, function(x){
  r <- random_binary_base(ncol(chaillou))
  t_chaillou <- test_bootstrap_all(chaillou, nb_cluster = 15, nb_axe = 16, nb_train=20, type="comptage", base_binaire = r)
  t_chaillou$all
})

save(aleatoire_chaillou, aleatoire_vacher, file = "aleatoire.RData")

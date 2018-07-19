source("Rscript/read_metagenomic_data.R")
source("Rscript/coda.R")
source("Rscript/graph.R")
source("Rscript/comparaison_clustering.R")
source("Rscript/bootstrap.R")
source("Rscript/test_bootstrap.R")
source("Rscript/tree_phyloseq.R")
source("Rscript/apprentissage_supervise.R")



##### bootstrap supervise

## donne comptage
c_super <- test_bootstrap_supervise(chaillou, metadata_chaillou$EnvType, type="comptage", nb_train = 5)
c_super$all
r_super <- test_bootstrap_supervise(ravel, metadata_ravel$CST, type="comptage")
r_super$all
l_super <- test_bootstrap_supervise(liver_500, metadata_liver$status, type="comptage")
l_super$all
m_super <- test_bootstrap_supervise(mach_500, metadata_mach$Weaned, type="comptage")
m_super$all


## donne ilr
c_super_ilr <- test_bootstrap_supervise(chaillou, metadata_chaillou$EnvType, type="ilr")
c_super_ilr$all
r_super_ilr <- test_bootstrap_supervise(ravel, metadata_ravel$CST, type="ilr")
r_super_ilr$all
l_super_ilr <- test_bootstrap_supervise(liver_500, metadata_liver$status, type="ilr")
l_super_ilr$all
m_super_ilr <- test_bootstrap_supervise(mach_500, metadata_mach$Weaned, type="ilr")
m_super_ilr$all


# validation croise
v_chaillou <- validation_croise(chaillou, metadata_chaillou$EnvType)
v_chaillou
v_ravel <- validation_croise(ravel, metadata_ravel$CST)
v_ravel
v_liver_500 <- validation_croise(liver_500, metadata_liver$status)
v_liver_500
v_mach_500 <- validation_croise(mach_500, metadata_mach$Weaned)
v_mach_500




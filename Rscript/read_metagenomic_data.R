library(biomformat)
library(phyloseq)
library(ape)
library(dplyr)
library(magrittr)

source("Rscript/coda.R")
## Short template on how to import metagenomic data

## Count data tables (otus in line, samples in column)

## Read a biom file
chaillou <- read_biom("data/chaillou/chaillou.biom") %>% biom_data() %>% as.matrix()
mach <- read_biom("data/mach/kinetic.biom") %>% biom_data() %>% as.matrix()
ravel <- read_biom("data/ravel/ravel.biom") %>% biom_data() %>% as.matrix()



### Read RDS file
liver <- readRDS("data/liver_qin_nan/LC_24042018.mgs.matrix.integer.RDS")[[1]]
metadata_liver <- readRDS("data/liver_qin_nan/metadata.RDS")


## Read a table
load("data/vacher/oaks.RData")
vacher <- Data$count %>% t()
metadata_vacher <- Data$covariates



## Read metadata file
metadata_chaillou <- read.table(file = 'data/chaillou/sample_metadata.tsv', sep = '\t', header = TRUE)
metadata_ravel <- read.table(file = 'data/ravel/sample_metadata.tsv', sep = '\t', header = TRUE)
metadata_mach<- read.table(file = 'data/mach/kinetic_sample_metadata.tsv', sep = '\t', header = TRUE)



###### transpose
liver <- t(liver)
chaillou <- t(chaillou)
mach <- t(mach)
ravel <- t(ravel)
vacher <- t(vacher)


##mise des donnees dans le meme ordre
chaillou <- chaillou[order(names(chaillou[,1])),]
metadata_chaillou<- metadata_chaillou[order(metadata_chaillou[,1]),]

mach <- mach[order(names(mach[,1])),]
metadata_mach<- metadata_mach[order(metadata_mach[,1]),]

ravel <- ravel[order(names(ravel[,1])),]
metadata_ravel<- metadata_ravel[order(metadata_ravel[,1]),]

liver <- liver[order(rownames(liver)),]
metadata_liver <- metadata_liver[order(rownames(metadata_liver)),]


#### selection otus les plus abondants
mach_sum <- apply(mach, 2, function(x){sum(x>0)})

mach_500 <- mach[,which(mach_sum>=nrow(mach)*0.05)]

liver_sum <- apply(liver, 2, function(x){sum(x>0)})

liver_500 <- liver[,which(liver_sum>=nrow(liver)*0.05)]



###### phyloseq
delete_tip_tree <- function(tree, data){
  p <- tree$tip.label[-which(colnames(data) %in% tree$tip.label)]
  drop.tip(tree, p)
}

tree_chaillou <- read_tree("data/chaillou/tree.nwk") %>% delete_tip_tree(chaillou)
tree_mach <- read_tree("data/mach/tree.nwk") %>% delete_tip_tree(mach)
tree_mach_500 <- read_tree("data/mach/tree.nwk") %>% delete_tip_tree(mach_500)



library(biomformat)
library(phyloseq)
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
liver <- (readRDS("data/liver_qin_nan/microbiota.abundance.2018.RDS")*1e8) %>% round()
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


#### 500 otus les plus abondants
mach_MAP <- MAP(mach)

mean_otus_mach <- rbind(apply(mach_MAP, 2, mean), 1:ncol(mach))
mean_otus_mach <- mean_otus_mach[, order(mean_otus_mach[1,], decreasing = TRUE) ]
mach_500 <- mach[,mean_otus_mach[2,1:500]]

liver_MAP <- MAP(liver)

mean_otus_liver <- rbind(apply(liver_MAP, 2, mean), 1:ncol(liver))
mean_otus_liver <- mean_otus_liver[, order(mean_otus_liver[1,], decreasing = TRUE) ]
liver_500 <- liver[,mean_otus_liver[2,1:500]]



###### phyloseq

tree_chaillou <- read_tree("data/chaillou/tree.nwk")
tree_mach <- read_tree("data/mach/tree.nwk")



library(biomformat)
library(magrittr)
## Short template on how to import metagenomic data

## Count data tables (otus in line, samples in column)

## Read a biom file
chaillou <- read_biom("data/chaillou/chaillou.biom") %>% biom_data() %>% as.matrix()
mach <- read_biom("data/mach/kinetic.biom") %>% biom_data() %>% as.matrix()
ravel <- read_biom("data/ravel/ravel.biom") %>% biom_data() %>% as.matrix()

## Read a table
load("data/vacher/oaks.RData")
vacher <- Data$count %>% t()



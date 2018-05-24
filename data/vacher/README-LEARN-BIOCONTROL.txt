% README
% 2017-05-16
% PAUVERT Charlie

Objective: Reconstruct grapevine leaf microbial networks

This microbial network reconstruction will give us a deeper insight into the biotic environment of pathogen species --i.e the pathobiome -- and will address the following questions:

1. Which fungi interacts with Erysiphe necator (causing grapevine powdery mildew) ?
2. Which fungi interacts with Plasmopara viticola (causing grapevine downy mildew) ?
3. Are these interactions impacted by the treatment (with or without cover crops in the vineyard) or the date?


In order to answer these questions, 288 leaf samples were collected and sequenced.
Among them, two negative control samples and two positive control samples.

After sequencing, raw sequences were processed as follows:

* Forward and reverse paired end merging
* Stringent sequence quality filter
* ITS1 region extraction
* Sequence clustering at 97% similarity
* Chimeric OTU detection (and removal) with two methods: de novo and reference based
* Taxonomic assignation with UNITE v6 and RDP Classifier
* OTU table generation

Controls samples were used to filter out OTU or to threshold their abundance to 0 according to Galan et al. 2016.
In short, negative control samples allows us to calibrate detection threshold and cross contamination.
Positive control samples allows to determine the ratio of wrongly assigned sequences to samples and hence to level down some OTU abundance to 0.

OTU with less than 10 sequences in total were also removed following Brown et al. 2015.
OTU independantly assigned to plants with BLAST were discarded.
OTU that were not classified by either methods (RDP or BLAST) were considered spurious and removed.



These stringent filters on OTU led some environmental samples to be without any or few sequences. 
They were kept for clarity but a resequencing is underway to investigate these low DNA samples. 




Hence, you will find the OTU table (739 OTU x 284 samples), a sample metadata table, and the OTU taxonomy assignments.


# OTU table

OTU identifiers are row names and sample names are column names.

The two most abundant OTU are the two pathogens and their OTU identifiers are the following:

* Erysiphe necator: 

* Plasmopara viticola: 


# Sample metadata

The sample table contains (1) the sampling design variable (Treatment, Date and Replicate), and (2)  two types of covariates corresponding to two types of ecological processes (selection by the plant and by weather conditions).
Hence, these columns corresponds to the following variables:

* Sampling design:
      - Treatment: CoverCrop or NoCoverCrop
      - Date: corresponding to the three sampling dates (2016-06-01, 2016-06-23 and 2016-07-05)
      - Replicate: 1 to 4
* Selection exerted by the plant (leaf traits):
      - Leaf Age in days
      - SLA: Specific Leaf Area (in meter square / kg) computed as the Leaf surface / Dry weight ratio
      - LeafWaterContent in % computed as one minus the ratio Dry weight / Fresh weight
* Selection exerted by weather conditions:
      - Mean and Variance Temperature (in °C) from leaf flush to sampling date
      - Mean and Variance Vapor Pressure Deficit (VPD in Pascal)  from leaf flush to sampling date
      - Number of rainy days from leaf flush to sampling date
* Sample sequence number before filtering in order to account for diffence in sampling effort


# Taxonomy table

You will also find enclosed a taxonomic classification table of the OTUs with the 7 ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species).


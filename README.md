# MangroveMicrobiomeProject
This repository contains the code run in the publication [[to be announced - in revision process]].

The publication deals with microbial community structures present in a mangrove in Paranaguá Bay, Southern Brazil. We unravel the most common community structures present, then establish the environmental factors that correlate with such structures, and later compare them with community from other biomes such as sea sediment, estuarine water and forest soils.

Files:
- functions.R: contains code for the functions common to both parts of the study.
-- Part 1 - 16S rRNA gene based
--- environmental.tsv: contains the chemical composition of the sediments.
--- metadata.tsv: contains metadata information.
--- functions.R: contains code for the functions exclusive to this part of the study.
--- 16S.R: contains the code to reproduce the results and the plots present at the publication.
--- otu_table.tsv: contains the OTU table, outputted by the Qiime pipeline mentioned in the publication.
--- table_mc33942_sorted_L6.tsv: contains the proportion of each genera per sample, outputted by the Qiime pipeline mentioned in the publication.
--- rep_set.tre: contains the phylogenetic tree, outputted by the Qiime pipeline mentioned in the publication.
-- Part 2 - WGS based
--- functions.R: contains code for the functions exclusive to this part of the study.
--- pre-WGS.R: contains code to retrieve the BIOM (JSON-like) files containing count information for all metagenomes addressed in this part of the study, from the MG-RAST webserver.
--- WGS.R: contains the code to reproduce the results and the plots present at the publication.
--- metadata.tsv: contains metadata information.
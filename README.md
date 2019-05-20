# MangroveMicrobiomeProject
This repository contains the code run in the publication \[\[to be announced - in revision process\]\].

The publication deals with microbial community structures present in a mangrove in Paranaguá Bay, Southern Brazil. We unravel the most common community structures present, then establish the environmental factors that correlate with such structures, and later compare them with community from other biomes such as sea sediment, estuarine water and forest soils.

Files:

\- functions.R: contains code for the functions common to both parts of the study.<br>
\-\- Part 1 - 16S rRNA gene based<br>
\-\-\- environmental.tsv: contains the chemical composition of the sediments.<br>
\-\-\- metadata.tsv: contains metadata information.<br>
\-\-\- functions.R: contains code for the functions exclusive to this part of the study.<br>
\-\-\- 16S.R: contains the code to reproduce the results and the plots present at the publication.<br>
\-\-\- otu_table.tsv: contains the OTU table, outputted by the Qiime pipeline mentioned in the publication.<br>
\-\-\- table_mc33942_sorted_L6.tsv: contains the proportion of each genera per sample, outputted by the Qiime pipeline mentioned in the publication.<br>
\-\-\- rep_set.tre: contains the phylogenetic tree, outputted by the Qiime pipeline mentioned in the publication.<br>
\-\- Part 2 - WGS based<br>
\-\-\- functions.R: contains code for the functions exclusive to this part of the study.<br>
\-\-\- pre-WGS.R: contains code to retrieve the BIOM (JSON-like) files containing count information for all metagenomes addressed in this part of the study, from the MG-RAST webserver.<br>
\-\-\- WGS.R: contains the code to reproduce the results and the plots present at the publication.<br>
\-\-\- metadata.tsv: contains metadata information.<br>

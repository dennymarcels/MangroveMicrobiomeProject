# This workflow is used to retrieve the BIOM (JSON-like) files, containing taxonomic and functional count tables for all metagenomes to be analyzed, from the MG-RAST webserver

library(rjson)

categories <- list(c("organism", "domain", "RefSeq"), c("organism", "phylum", "RefSeq"), c("organism", "class", "RefSeq"), c("organism", "order", "RefSeq"), c("organism", "family", "RefSeq"), c("organism", "genus", "RefSeq"), c("function", "level1", "Subsystems"), c("function", "level2", "Subsystems"), c("function", "level3", "Subsystems"), c("function", "function", "Subsystems"))
metagenomes.id <- read.delim('Part 2 - WGS based/metadata.tsv', header = T, stringsAsFactors = F)$sample[1:17]

## Prepare a list to retrieve the files from

n <- 1; submit <- list()
for (cat in categories){
    for (i in metagenomes.id){
        url <- paste0('http://api.metagenomics.anl.gov/1/matrix/', cat[[1]], '?id=mgm', i, '&group_level=', cat[[2]], '&source=', cat[[3]], '&result_type=abundance&evalue=5')
        cat(n, "/", (length(metagenomes.id)+length(mangroves.id))*length(categories), "\n", sep="")
        submit[[n]] <- c(n, i, cat[[1]], cat[[2]], cat[[3]], fromJSON(file=url)$url)
        n <- n+1
    }
}
rm(n, cat, i, url)

## Download the metagenomes contained in the list

war <- list(); m <- 1
if(!dir.exists('Part 2 - WGS based/biom')) dir.create('Part 2 - WGS based/biom')
for (sub in submit){
    tryCatch({
        status <- "processing"
        cat(sub[1], "/", length(submit), "\n", sep="")
        while(status != "done"){
            status <- fromJSON(file = sub[6])$status
            if(status != "done"){Sys.sleep(5)
            }else{
                download.file(sub[6], destfile=paste0('Part 2 - WGS based/biom/', sub[4], '_', sub[2],'.biom', sep=""), method="internal")
                closeAllConnections()
            }
        }
    }, warning = function(w){    
        cat("WARNING :", conditionMessage(w), "\n")
        war[[m]] <<- sub
        m <<- m+1
    }, error=function(e){
        cat("ERROR :", conditionMessage(e), "\n")
        war[[m]] <<- sub
        m <<- m+1
    }
    )
    warning <- war
}
rm(war, m, sub, status)
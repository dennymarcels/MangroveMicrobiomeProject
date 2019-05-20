# FUNCTIONS EXCLUSIVE TO THE WGS BASED WORKFLOW

# Function to retrieve the taxonomic/functional information from the BIOM files

retrieve <- function(file){
    json <- fromJSON(file = file)
    len_names <- numeric()
    for(n in 1:length(json$data$rows)){
        len_names <- c(len_names, length(json$data$rows[[n]]$metadata$hierarchy))    
    }
    len_names <- max(len_names)
    n <- 1
    while(length(json$data$rows[[n]]$metadata$hierarchy) < len_names){
        n <- n+1
    }
    names <- rev(sort(names(json$data$rows[[n]]$metadata$hierarchy)))
    biom <- data.frame()
    for(n in 1:length(json$data$rows)){
        this_biom <- c('id' = json$data$rows[[n]]$id)
        if(length(json$data$rows[[n]]$metadata) == 0){
            this_biom <- c(this_biom, rep('', length(names)))
            this_biom <- setNames(this_biom, c('id', names))
            this_biom <- t(this_biom)
        }
        else{
            for(name in names){
                this_this_biom <- c(json$data$rows[[n]]$metadata$hierarchy[name])
                if(this_this_biom == 'NULL'){
                    which = which(names == name)
                    this_this_biom = paste0('undefined (from ', json$data$rows[[n]]$metadata$hierarchy[names[which+1]], ')')
                }
                this_biom <- c(this_biom, this_this_biom)
            }
        }
        names(this_biom) <- c('id', names)
        biom <- rbind(biom, data.frame(this_biom), stringsAsFactors = T)
    }
    data <- as.numeric(json$data$data)
    table <- cbind.data.frame(biom, data, row.names=as.character(biom$id))
    mg_name <- substr(json$data$columns[[1]]$id, 4, 12)
    colnames(table)[length(colnames(table))] <- mg_name
    return(table)
}

# Function to merge all information retrieved from the BIOM files in a single data frame

multmerge <- function(pat, metagenomes){
    filenames <- lapply(metagenomes, function(x){x <- paste0('Part 2 - WGS based/biom/', pat, '_', x, '.biom')})
    datalist <- lapply(filenames, function(x){x <- retrieve(x)})
    table <- Reduce(function(x, y) {
        table <- merge(x, y, all=T)
    }, datalist)
    return(table)
}

# Function to consolidate the genus count table at all levels explicited in 'taxonomy', then plot boxplots of all these levels, separating sequences by domain

boxplots <- function(taxotable, functable, taxonomy, functional){
    
    # Taxonomic boxplots
    for(tax in taxonomy){
        tableAll <- consolidate(taxotable, tax)
        domains <- unique(tableAll[, 1])
        domains <- factor(domains, levels = c('unclassified', 'Viruses', 'Archaea', 'Bacteria', 'Eukaryota'))
        domains <- sort(domains)
        ord_domains <- unclass(factor(domains, levels = unique(tableAll[, 1])))
        if(tax == 'domain'){
            rownames(tableAll) <- tableAll[, 1]
            tableAll <- tableAll[ord_domains, ]
        }else{
            unclassified <- grepl('^unclassified', tableAll[, 2])
            tableAll[, 1][unclassified] <- 'unclassified'
            tables <- list()
            for(domain in domains){
                table <- tableAll[tableAll$domain == domain, , drop = F]
                n <- ifelse(domain %in% c('Bacteria', 'Eukaryota'), 10, 5)
                table <- table[order(rowMeans(table[, 3:ncol(table)]), decreasing = T), ][n:1, ]
                table <- table[!grepl('NA', rownames(table)), ]
                tables[[domain]] <- table
            }
            tableAll <- bind_rows(tables)
            rownames(tableAll) <- tableAll[, 2]
            tableAll <- tableAll[, -2]
        }
        num_each <- table(tableAll$domain)[ord_domains]
        num_each <- num_each[!is.na(num_each)]
        tableAll <- t(tableAll[, -1])
        
        png(paste0('Part 2 - WGS based/boxplot_mang_', tax, '.png'), height = 450, width = 600)
        par(mar = c(4, 2, 1, 17.5))
        boxplotting(tableAll, tax, num_each)
        if(tax != 'domain'){
            domains <- domains[domains %in% names(num_each)]
            mtext(sapply(domains, function(x) paste0(x, ' |')), at = cumsum(num_each)+0.5, side = 2, adj = 1, line = 0.5)
        }
        dev.off()
    }
    
    # Functional boxplots
    for(fun in functional){
        tableAll <- consolidate(functable, fun)
        tableAll <- tableAll[order(rowMeans(tableAll[, 2:ncol(tableAll)]), decreasing = T), ][20:1, ]
        rownames(tableAll) <- tableAll[, 1]
        tableAll <- t(tableAll[, -1])
        
        png(paste0('Part 2 - WGS based/boxplot_', fun, '_mang.png'), height = 350, width = 750)
        par(mar = c(4, 1, 1, 35))
        boxplotting(tableAll, fun)
        dev.off()
    }
}

# Function to consolidate the count table by environmental factor (taken from the chosen metadata column). While consolidating, it will rarefy among the corresponding metagenomes, complete empty counts with 0, remove empty categories and convert counts to proportions. It will return a list with filtered metadata and consolidated tables by each hierarchical level

consolidate_by_factor <- function(table_taxo_WGS, table_funct_WGS, table_16S = NULL, metadata, factor){
    metagenomes <- as.character(metadata[, 'sample'][!is.na(metadata[, factor])])
    
    if(factor == 'sequencing method'){
        taxo1 <- table_16S[, c(taxonomy, metagenomes[metagenomes %in% colnames(table_16S)])]
        taxo1[, taxonomy] <- sapply(taxo1[, taxonomy][, 1:ncol(taxo1[, taxonomy])], 
                                    function(x) {
                                        x <- gsub('\\[', '', x)
                                        x <- gsub('\\]', '', x)
                                        return(x)
                                    }
        )
        taxo2 <- table_taxo_WGS[, c(taxonomy, metagenomes[metagenomes %in% colnames(table_taxo_WGS)])]
        domains <- unique(taxo2[, 'domain'])
        domains <- as.character(domains[unique(domains) %in% unique(taxo1[, 'domain'])])
        taxo2 <- taxo2[taxo2[, 'domain'] %in% domains, ]
        taxo <- merge(taxo1, taxo2, all = T)
        taxo[is.na(taxo)] <- 0
    }
    else{
        taxo <- table_taxo_WGS[, c(taxonomy, metagenomes)]
    }
    
    matrix <- taxo[, metagenomes]
    set.seed(1)
    matrix <- t(Rarefy(t(matrix))$otu.tab.rff)
    taxo[, metagenomes] <- matrix
    
    if(factor != 'sequencing method'){
        funct <- table_funct_WGS[, c(functional, metagenomes)]
        matrix <- funct[, metagenomes]
        set.seed(1)
        matrix <- t(Rarefy(t(matrix))$otu.tab.rff)
        funct[, metagenomes] <- matrix
    }
    
    consolidated <- list()
    
    temp <- metadata[(metadata[, 'sample'] %in% metagenomes), c('sample', 'biome', factor)]
    colnames(temp)[3] <- 'groups'
    consolidated[['metadata']] <- temp
    
    for(tax in taxonomy){
        temp <- taxo %>% group_by_at(tax) %>% summarise_if(is.numeric, sum) %>% as.data.frame
        temp <- temp[rowSums(temp[, 2:ncol(temp)]) > 0, ]
        for(col in 2:ncol(temp)){
            temp[, col] <- temp[, col]/sum(temp[, col]) 
        }
        consolidated[[tax]] <- temp
    }
    
    if(factor != 'sequencing method'){
        for(fun in functional){
            temp <- funct %>% group_by_at(vars(fun)) %>% summarise_if(is.numeric, sum) %>% as.data.frame
            temp <- temp[rowSums(temp[, 2:ncol(temp)]) > 0, ]
            for(col in 2:ncol(temp)){
                temp[, col] <- temp[, col]/sum(temp[, col]) 
            }
            consolidated[[fun]] <- temp
        }
    }
    
    return(consolidated)
}

# Function to perform statistical analysis based on factor and hierarchical level. Performs PERMANOVA, pairwise PERMANOVA (when more then two features) and t-tests (when two features)

statistical_analysis <- function(factor, level){
    metadata <- comp[[factor]][['metadata']]
    table <- comp[[factor]][[level]]
    rownames(table) <- table[, 1]
    table <- t(table[, -1])
    
    set.seed(1)
    adonis <- adonis(table ~ groups, metadata)
    
    if(length(unique(metadata[, 'groups'])) > 2){
        factors <- as.character(unique(metadata[, 'groups']))
        mang_factor <- metadata[metadata[, 'biome'] == 'mangrove', 'groups']
        mang_factor <- names(rev(sort(table(mang_factor)))[1])
        factors <- factors[-which(factors == mang_factor)]
        
        pairwise <- data.frame()
        for(factor in factors){
            mdt <- metadata[metadata[, 'groups'] %in% c(mang_factor, factor), ]
            tbl <- table[mdt[, 'sample'], ]
            
            set.seed(1)
            adonis2 <- adonis(tbl ~ mdt[, 'groups'], mdt)
            df <- data.frame('factor' = factor,
                             'F.model' = adonis2$aov.tab$F.Model[1],
                             'R2' = adonis2$aov.tab$R2[1],
                             'p-value' = adonis2$aov.tab$`Pr(>F)`[1])
            pairwise <- rbind(pairwise, df)
        }
        
        pairwise <- list('compared_to' = mang_factor, 'pairwise' = pairwise)
        comparisons <- list('adonis' = as.data.frame(adonis$aov.tab), 'pairwise' = pairwise)
        
        return(comparisons)
    }
    
    # Groups to compare
    
    factors <- as.character(unique(metadata[, 'groups']))
    if(factors[2] == 'mangrove'){factors <- rev(factors)}
    
    ttests <- perform_ttests(table, metadata, factors)
    
    comparisons <- list('adonis' = as.data.frame(adonis$aov.tab), 'ttests' = ttests)
    
    return(comparisons)
}

# Function to generate a summary of the PERMANOVA results contained in a list of results

summarise_permanova <- function(results){
    summaries <- data.frame()
    for(factor in names(results)){
        for(level in names(results[[factor]])){
            F.model <- tryCatch(
                {results[[factor]][[level]][['adonis']]$F.Model[1]
                }, error = function(cond){NA})
            p.value <- tryCatch(
                {results[[factor]][[level]][['adonis']]$`Pr(>F)`[1]
                }, error = function(cond){NA})
            summary <- data.frame(factor = factor,
                                  level = level, 
                                  F.model = F.model,
                                  p.value = p.value
            )
            summaries <- rbind(summaries, summary)
        }
    }
    summaries <- summaries %>% arrange(factor)
    print(summaries)
    return(summaries)
}

# Function to truncate features names to 35 characters max (helps plotting)

truncateNames <- function(name){
    if(nchar(name) > 40){
        name <- unlist(strsplit(name, ''))
        name <- paste0(c(name[1:35], '[...]'), collapse = '')
    }
    return(name)
}

# Function to append a discriminator for those truncated features names that, due to truncating, end up being equal

dealRepeated <- function(names){
    repeated <- names(table(names)[table(names) > 1])
    for(r in repeated){
        which <- which(names == r)
        i <- 1
        for(w in which){
            names[w] <- sub('$', paste0('(', i, ')'), names[w])
            i <- i+1
        }
    }
    return(names)
}

# Function to calculate similarity indexes between groups, from distance metrics between individual metagenomes

similarity <- function(category.dist){
    category.sim <- as.matrix(1-category.dist)
    diag(category.sim) <- 1
    
    biomes <- factors.for.col$biome
    uniq.biomes <- unique(biomes)
    
    sim <- data.frame()
    for(biome in uniq.biomes){
        row <- biome
        which.row <- which(biomes == biome)
        for(biome in uniq.biomes){
            col <- biome
            which.col <- which(biomes == biome)
            category.sum <- 0
            for(i in which.row){
                for(j in which.col){
                    category.sum <- category.sum + category.sim[i, j]
                }
            }
            sim[row, col] <- category.sum/(length(which.row)*length(which.col))
        }
    }
    return(sim)
}
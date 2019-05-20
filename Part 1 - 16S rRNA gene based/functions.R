# FUNCTIONS EXCLUSIVE TO THE 16S rRNA GENE BASED WORKFLOW

# Function to adjust the name in the genus abundance table

adjustname <- function(names){
    names <- gsub('^Unassigned', 'unclassified', names) # 'Unassigned' at the beginning of a string will change to 'unclassified'
    names <- gsub('Other', 'unclassified', names) # All 'Other' strings will change to 'unclassified'
    names <- gsub('__$', '__unclassified', names) # All empty information at the end of a string will change to 'unclassified'
    names <- gsub('__;', '__unclassified;', names) # All empty information overall will change to 'unclassified'
    names <- gsub('(.__)(*)', '\\2', names) # Removing all prefix denoting hierarchical level, and the __ separator
    return(names)
}

# Function to split the column containing the whole taxonomic tree labels into separate columns, one for each taxonomic level

split_taxa_in_columns <- function(table){
    taxanames <- rownames(table)
    taxanames <- adjustname(taxanames)
    table <- as.data.frame(table)
    rownames(table) <- NULL
    
    splits <- strsplit(taxanames, ';')
    
    table$domain <- unlist(lapply(splits, '[[', 1))
    table$phylum <- unlist(lapply(splits, '[[', 2))
    table$class <- unlist(lapply(splits, '[[', 3))
    table$order <- unlist(lapply(splits, '[[', 4))
    table$family <- unlist(lapply(splits, '[[', 5))
    table$genus <- unlist(lapply(splits, '[[', 6))
    
    table <- table[, c(31:ncol(table), 1:30)]
    return(table)
}

# Function to append the closest taxonomic information to those categories labeled as 'unclassified' 

preserve_tax_info <- function(table){
    for(j in 6:2){
        for(i in 1:length(table[, j])){
            if(table[i, j] == 'unclassified'){
                row <- table[i, 1:j]
                where_not_uncl <- which(row != 'unclassified')
                if(length(where_not_uncl) > 0){
                    higher_order <- row[max(where_not_uncl)]
                    table[i, 2:6] <- gsub('unclassified', paste0('unclassified (from ', higher_order, ')'), table[i, 2:6]) 
                }
            }
        }
    }
    return(table)
}

# Function to process the genus count table, consolidating it at all levels explicited in 'taxonomy', then plotting boxplots of all these levels, separating Bacteria, Archaea and unclassified sequences

boxplots <- function(table, taxonomy){
    for(tax in taxonomy){
        tableAll <- consolidate(table, tax)
        domain <- tableAll[, 1]
        if(tax == 'domain'){
            rownames(tableAll) <- tableAll[, 1]
        }else{
            rownames(tableAll) <- tableAll[, 2]
            tableAll <- tableAll[, -2]
        }
        
        tableBacteria <- tableAll[domain != 'Archaea', , drop = F]
        tableBactC <- tableBacteria[!(grepl('unclassified', rownames(tableBacteria))), , drop = F]
        tableBactU <-tableBacteria[grepl('unclassified', rownames(tableBacteria)), , drop = F]
        tableArchaea <- tableAll[domain == 'Archaea', , drop = F]
        
        tableAll <- rbind(tableArchaea[order(rowMeans(tableArchaea[, 2:ncol(tableArchaea)]), decreasing=T), ][5:1, ], # limit Archaeal taxa to 5
                          tableBactU[order(rowMeans(tableBactU[, 2:ncol(tableBactU)]), decreasing=T), ][5:1, ], # limit unclassified taxa to 5
                          tableBactC[order(rowMeans(tableBactC[, 2:ncol(tableBactC)]), decreasing=T), ][10:1, ] # limit classified Bacterial taxa to 10
        )
        tableAll <- tableAll[!grepl('NA', rownames(tableAll)), ]
        num_each <- c(min(nrow(tableArchaea), 5), min(nrow(tableBactU), 5), min(nrow(tableBactC),10))
        tableAll <- t(tableAll[, -1])
        
        png(paste0('Part 1 - 16S rRNA gene based/boxplot_', tax, '.png'), height = 350, width = 600)
        par(mar = c(4, 2, 1, 17.5))
        boxplotting(tableAll, tax, num_each)
        if(tax != 'domain'){
            mtext(c('Archaea |', 'unc |', 'Bacteria'), at = c(num_each[1], num_each[1]+num_each[2], sum(num_each))+0.5, side = 2, adj = 1, line = 0.5)
        }
        dev.off()
    }
}

# Function to consolidate the counts table at the desired taxonomic level, then formats the table to be used for hierarchical clustering

prepare_for_hclust <- function(table, level){
    table <- consolidate(table, level)
    rownames(table) <- table[, 2]
    table <- t(table[, -c(1, 2)])
    table <- table[, order(colMeans(table), decreasing=T)]
    return(as.data.frame(table))
}

# Function to plot the first two dimensions of a RDA object with the chosen scaling

plot_rda <- function(rda_object, scaling){
    
    ### Getting coordinates
    
    plot <- plot(rda_object, scaling = scaling)
    
    species <- plot$species
    sites <- plot$sites
    centroids <- plot$centroids
    biplot <- plot$biplot * attr(plot$biplot, 'arrow.mul')
    rownames(biplot) <- gsub('`', '', rownames(biplot))
    
    ### Setting plot parameters
    
    xlim <- c(1.2*min(c(species[, 1], sites[, 1], centroids[ ,1], biplot[, 1])), 
              1.2*max(c(species[, 1], sites[, 1], centroids[ ,1], biplot[, 1])))
    ylim <- c(1.1*min(c(species[, 2], sites[, 2], centroids[ ,2], biplot[, 2])), 
              1.1*max(c(species[, 2], sites[, 2], centroids[ ,2], biplot[, 2])))
    
    comp <- c(summary(dbrda)$cont$importance[2]*100, summary(dbrda)$cont$importance[5]*100)
    cons <- c(summary(dbrda)$concont$importance[2]*100, summary(dbrda)$concont$importance[5]*100)
    
    xlab <- paste0('CAP1 (', round(comp[1],1), '% / ', round(cons[1],1), '%)')
    ylab <- paste0('CAP2 (', round(comp[2],1), '% / ', round(cons[2],1), '%)')
    
    ### Setting file
    
    png(paste0('Part 1 - 16S rRNA gene based/rda_scaling_', i, '.png'), width = 1000, height = 1000, res = 200)
    
    ### Ploting empty plot
    
    par(mgp = c(1.5, 0.5, 0), mar = c(3, 3, 1, 1))
    plot(0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', cex.lab = 0.7, type = 'none')
    x = round(xlim); s = seq(x[1], x[2], 1); axis(side = 1, at = s, cex.axis = 0.7)
    y = round(ylim); s = seq(y[1], y[2], 1); axis(side = 2, at = s, cex.axis = 0.7)
    abline(h = 0, v = 0, lty = 3)
    
    ### Including centroids for each sample
    
    ordispider(sites, groups = env$sample, col = 'lightgray')
    text_coords <- as.data.frame(sites) %>% group_by(env$sample) %>% summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2))
    text_coords <- as.data.frame(text_coords)
    rownames(text_coords) <- text_coords$`env$sample`; text_coords <- text_coords[, -1]
    text(text_coords, labels = rownames(text_coords), cex = 0.6)
    
    ### Highlighting sample L2L, and sampling spots 1 & 2 for H
    
    text(sites['L2L', ][1], sites['L2L', ][2], label = 'L2L', col = 'darkgray', cex = 0.6)
    text_coords <- as.data.frame(sites) %>% group_by(env$salinity, env$sample, env$replicate) %>% filter(`env$salinity` == 'H', `env$replicate` %in% c(1, 2))
    text(text_coords$CAP1, text_coords$CAP2, labels = text_coords$`env$replicate`, col = 'darkgray', cex = 0.6)
    
    ### Including centrois for each salinity
    
    ordihull(sites, group = fact$salinity, lty = 2, col = adjustcolor('purple', alpha.f = 0.5))
    text(centroids, labels = rownames(centroids), cex = 0.8, col = adjustcolor('purple', alpha.f = 0.5))
    text(centroids, labels = rownames(centroids), cex = 0.8, col = adjustcolor('purple', alpha.f = 0.5))
    
    ### Including species
    
    text(species, labels = rownames(species), cex = 0.6, col = 'darkgreen')
    
    ### Including environmental factors
    
    arrows(0, 0, biplot[, 1], biplot[, 2], length = 0.1, angle = 15, col = 'red')
    text(biplot + 0.02, labels = rownames(biplot), cex = 0.6, col = 'red', adj = 0)
    
    ### Saving file
    
    dev.off()
}
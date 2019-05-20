# 1 Setting environment ----

# 1.1 Loading dependencies ----

require(dplyr)
require(GUniFrac)
require(rjson)
require(mratios)
require(scales)
require(ggplot2)
require(gridExtra)
require(grid)
require(heatmap.plus)
require(randomForest)
require(candisc)
source('functions.R')
source('Part 2 - WGS based/functions.R')

# 1.2 Defining global variables ----

study <- 'WGS'
taxonomy <- c('domain', 'phylum', 'class', 'order', 'family', 'genus')
functional <- c('level1', 'level2', 'level3', 'level4')

# 2 Getting count ('BIG') tables ----

# 2.1 Getting BIG TABLE for 16S samples ----

taxo16S <- read.delim('Part 1 - 16S rRNA gene based/taxonomic_counts_16S.tsv', header = T, check.names = F)

# 2.2 Getting BIG TABLE for WGS samples (taxonomy) ----

metadata <- read.delim('Part 2 - WGS based/metadata.tsv', sep = '\t', header = T, check.names = F)

# Metadata contains information for all comparisons to be addressed; for now, I will extract information from the first 17 rows, which correspond to the WGS data
mg_ids <- metadata$sample[1:17]
mg_names <- metadata$name[1:17]

# I will load in only the genus tables, as it contains the whole taxonomic information
taxoWGS <- multmerge('genus', mg_ids)

# Which ids are there where id != genus? These are entries with no taxonomic information besides what is in id. I will move this information to the 'genus' column
mis_ids_n <- which((as.character(taxoWGS$id) == as.character(taxoWGS$genus)) == F)
mis_ids <- as.character(taxoWGS$id[mis_ids_n])
taxoWGS$genus <- taxoWGS$id

# Reordering the columns in taxonomic order
taxoWGS <- taxoWGS[, c(6, 2, 7, 3, 5, 4, 8:ncol(taxoWGS))]

# Now I will standardize some string entries in the rows with no taxonomic information identified in the previous step
for(mis_id_n in mis_ids_n){
    label <- as.character(taxoWGS$genus[mis_id_n])
    
    # If the entry contains the text 'viruses', the domain will be identified as 'Viruses'
    if(grepl('viruses', label)){
        which_tax <- 1
        label <- 'Viruses'
    }
    
    # If not, I will extracted the '*' text contained in the 'derived from *' string, and try to match it to its corresponding taxonomic level
    else{
        label <- sub('^unclassified \\(derived from (.*)\\)', '\\1', label)
        for(i in 1:6){
            tax <- taxonomy[i]
            which_tax <- grep(paste0('^', label, '$'), taxoWGS[, tax])
            if(length(which_tax) > 0){
                which_tax <- i
                break
            }
        }
        
        # Some of the texts extracted have more info than the necessary after the first word and won't find a match, I will update the text to contain only the first word
        if(length(which_tax) == 0){
            label <- sub('^([^ ]+) .+', '\\1', label)
            for(i in 1:6){
                tax <- taxonomy[i]
                which_tax <- grep(paste0('^', label, '$'), taxoWGS[, tax])
                if(length(which_tax) > 0){
                    which_tax <- i
                    break
                }
            }
        }
        
        # If, nonetheless, no match is found, then I will attribute the label 'unclassifed'
        if(length(which_tax) == 0){
            which_tax <- 1
            label <- 'unclassified'
        }
    }
    
    # This string vector will contain the whole taxonomic info available for the match 
    info <- taxoWGS[taxoWGS[, which_tax] == label, 1:which_tax]
    if(which_tax == 1){
        info <- info[1]
    }
    else{
        info <- info[1, ]
    }
    
    # Finally, the taxonomic info is inserted where it was missing
    taxoWGS[mis_id_n, 1:which_tax] <- info
}

remaining <- taxoWGS[is.na(taxoWGS$domain), ]$genus
remaining_n <- rownames(taxoWGS[is.na(taxoWGS$domain), ])
print(remaining_n); print(remaining)

# There is no taxonomic information for these 4 genera (plus, obviously, the totally unclassified sequences) in the table, but I found it online. Let me complete these manually.
levels(taxoWGS$domain) <- c(levels(taxoWGS$domain), 'unclassified')
levels(taxoWGS$phylum) <- c(levels(taxoWGS$phylum), 'Negarnaviricota', 'unclassified')
levels(taxoWGS$class) <- c(levels(taxoWGS$class), 'Monjiricites', 'unclassified')
levels(taxoWGS$order) <- c(levels(taxoWGS$order), 'Psocoptera', 'Chondrosida', 'Mononegavirales', 'unclassified')
levels(taxoWGS$family) <- c(levels(taxoWGS$family), 'Marseilleviridae', 'Lepidopsocidae', 'Chondrillidae', 'unclassified')

taxoWGS[984, 1:5] <- c('Viruses', '', '', '', 'Marseilleviridae')
taxoWGS[1141, 1:5] <- c('Eukaryota', 'Arthropoda', 'Insecta', 'Psocoptera', 'Lepidopsocidae')
taxoWGS[1177, 1:5] <- c('Eukaryota', 'Porifera', 'Demospongiae', 'Chondrosida', 'Chondrillidae') 
taxoWGS[1180, 1:5] <- c('Viruses', 'Negarnaviricota', 'Monjiricites', 'Mononegavirales', '')

# Now I will change all taxonomic strings containing 'unclassified' to 'unclassified' alone
for(i in 1:6){
    taxoWGS[, i] <- sub('.*unclassified.*', 'unclassified', taxoWGS[, i])
}

# Empty strings will also change to 'unclassified'
taxoWGS[, 1:6][is.na(taxoWGS[, 1:6])] <- 'unclassified'
taxoWGS[, 1:6][taxoWGS[, 1:6] == ''] <- 'unclassified'

# Empty count values will change to 0
taxoWGS[is.na(taxoWGS)] <- 0

# taxoWGS contains unclassified taxa indicates simply as 'unclassified'. I would also want to preserve some taxonomic information for these unclassified taxa, appending to its label the first higher order identified for it, as in '(from higher_order)'
taxoWGS <- preserve_tax_info(taxoWGS)

write.table(taxoWGS, 'Part 2 - WGS based/taxonomic_counts_WGS.tsv', sep = '\t', row.names = F)

# 2.3 Getting BIG TABLE for WGS samples (functional) ----

# I will load in only the function (level4) tables, as it contains the whole functional hierarchy information
functWGS <- multmerge('function', mg_ids)

# Reordering the columns in hierarchical order
functWGS <- functWGS[, c(5:1, 6:ncol(functWGS))]

# Empty count values will change to 0
functWGS[is.na(functWGS)] <- 0

# functWGS contains some categories that have more than one higher order. I would want to preserve this information for these categories, appending to its label the higher order related to it, as in '(from higher_order)'
for(i in 5:2){
    levs <- colnames(functWGS)[1:5]
    
    # This table will contain the cases where one category has more than one immediately higher order
    table <- functWGS %>% select_at(c(i, i-1)) %>% unique %>% group_by_at(1) %>% summarise(n = n()) %>% filter(n != 1) %>% data.frame
    if(nrow(table) > 0){
 
        prev_lev <- as.character(functWGS[, levs[i-1]])
        lev <- functWGS[, levs[i]]
        
        # Cycling throught the 'repeated' categories and appending its corresponding higher level information
        for(j in 1:length(lev)){
            if(functWGS[j, i] %in% table[, 1]){
                sub <- gsub('(.*)', paste0('\\1 (from ', prev_lev[j], ')'), functWGS[j, i])
                levels(functWGS[, levs[i]]) <- c(levels(functWGS[, levs[i]]), sub)
                functWGS[j, i] <- sub
            }
        }
    }
}

write.table(functWGS, 'Part 2 - WGS based/functional_counts_WGS.tsv', sep = '\t', quote = F, row.names = F)

# 3 Taxonomic and functional mangrove profiles based on WGS data ----

# Now I will separate the mangrove WGS data
taxoMangWGS <- taxoWGS[, c(1:6, grep('mangrove', metadata$biome)+6)]
functMangWGS <- functWGS[, c(1:5, grep('mangrove', metadata$biome)+5)]

# 3.1 Generating boxplots for all taxonomic and functional levels ----

boxplots(taxoMangWGS, functMangWGS, taxonomy, functional)

# 4 Taxonomic and functional profiles for all comparisons desired ----

# Generating a comparison list for all environmental factors to be analyzed
comp <- list()
factors <- colnames(metadata[, 4:ncol(metadata)])
for(factor in factors){
    comp[[factor]] <- consolidate_by_factor(taxoWGS, functWGS, taxo16S, metadata, factor)
}

# 5 Statistical analysis ----

# Generating a list with results for PERMANOVA, pairwise PERMANOVA (where there are more than two groups) and pairwise t-tests (when factor has only two levels)
results = list()
for(factor in factors){
    for(level in taxonomy){
        results[[factor]][[level]] <- statistical_analysis(factor, level)
    }
    for(level in functional){
        if(factor == 'sequencing method'){next}
        results[[factor]][[level]] <- statistical_analysis(factor, level)
    }
}

# 5.1 Generating a summary with the PERMANOVA results ----
summaries <- summarise_permanova(results)

# 6 Pairwise comparisons between sequencing methods ----

# Plotting pairwise comparisons at all hierarchical levels possible
remove_unclassified <- F
factor <- 'sequencing method'

level <- 'domain'
tweak_sup <- -2.9
tweak_pval <- 1.6
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'phylum'
tweak_sup <- -1.5
tweak_pval <- 1.5
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'class'
tweak_sup <- 0
tweak_pval <- 1.56
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'order'
tweak_sup <- 1
tweak_pval <- 1.57
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'family'
tweak_sup <- 1
tweak_pval <- 1.57
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'genus'
tweak_sup <- 1
tweak_pval <- 1.58
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

# 7 Investigating mangrove data ----

# 7.1 Domain / WGS data ----

mang_domain <- comp$salinity$domain
rownames(mang_domain) <- mang_domain$domain
mang_domain <- mang_domain[, -1]
mang_domain$mean <- rowMeans(mang_domain)*100
mang_domain$sd <- apply(mang_domain[, 1:4], 1, sd)*100
print(mang_domain[, c('mean', 'sd')])

# 7.2 Other taxonomic levels / WGS data x 16S data ----

enrich <- list()
factor <- 'sequencing method'
for(level in taxonomy){
    enrich[[level]] <- enriched(results[[factor]][[level]][['ttests']])
}

# Taxa enriched for WGS compared to 16S
for(name in names(enrich)){
    g <- enrich[[name]][[1]]
    print(name)
    g <- g[abs(g[, 'diff']) >= 1, c('factor', 'diff', 'p')]
    print(names(enrich[[name]]))
    print(g)
}

# 7.3 Level1 / WGS data ----

mang_level1 <- comp$salinity$level1
rownames(mang_level1) <- mang_level1$level1
mang_level1 <- mang_level1[, -1]
mang_level1$mean <- rowMeans(mang_level1)*100
mang_level1$sd <- apply(mang_level1[, 1:4], 1, sd)*100
mang_level1 <- mang_level1[order(mang_level1$mean, decreasing = T), ]
mang_level1$cumsum <- cumsum(mang_level1$mean)
print(mang_level1[, c('mean', 'cumsum', 'sd')])

# 7.4 Salinity sectors / WGS data ----

# Max difference between salinity sectors at all taxonomic/functional levels

factor <- 'salinity'

enrich <- list()
for(level in taxonomy){
    enrich[[level]] <- enriched(results[[factor]][[level]][['ttests']])
}
print(names(results[[factor]][[level]][['ttests']]))
print(sort(sapply(enrich, function(x) x[[1]][order(abs(x[[1]]$diff), decreasing = T), ][1, 2])))

enrich <- list()
for(level in functional){
    enrich[[level]] <- enriched(results[[factor]][[level]][['ttests']])
}
print(sort(sapply(enrich, function(x) x[[1]][order(abs(x[[1]]$diff), decreasing = T), ][1, 2])))

# 8 Hierarchical clustering/heatmap of all levels ----

# Generating a color list
factors.for.col <- metadata[1:17, c('biome', 'salinity status', 'matrix type', 'stable organization')]
colours.biome <- rainbow(length(unique(factors.for.col$biome)))[factors.for.col$biome]
colours.salinity <- colorRampPalette(c('gray', 'black'))(length(unique(factors.for.col$`salinity status`)))[factors.for.col$`salinity status`]
colours.matrix <- colorRampPalette(c('brown', 'lightblue'))(length(unique(factors.for.col$`matrix type`)))[factors.for.col$`matrix type`]
colours.org <- gray.colors(length(unique(factors.for.col$`stable organization`)), start = 0)[factors.for.col$`stable organization`]

# To be used in heatmaps
ann.colors <- list('matrix type' = colours.matrix, 'salinity status' = colours.salinity, 'biome' = colours.biome)

# 8.1 Plotting ----
for(level in c(taxonomy, functional)){
    table <- comp[['biome']][[level]]
    rownames(table) <- table[, 1]
    table <- t(table[, -1])
    rownames(table) <- mg_names
    table <- table[, order(colMeans(table), decreasing=T)]
    dist <- vegdist(table, 'bray')
    clust <- hclust(dist, method = 'average')
    table <- simplify(table)

    plot_heatmap(table, level)
}

# 9 Internal/external similarity of mangroves to other biomes ----

# 9.1 Taxonomic similarities ----

# First I will calculate mean taxonomic similarities between biomes
tax.sim <- data.frame()
for(level in rev(taxonomy)){    
    table <- consolidate(taxoWGS, level)
    rownames(table) <- table[, level]
    table <- t(table[, -c(1:which(colnames(table) == level))])
    dist <- vegdist(table) # dissimilarities between metagenomes
    sim <- similarity(dist) # this derives group (mean) similarities from individual dissimilarities
    sim <- t(sim[ ,'mangrove', drop = F]) # keep only comparison to mangrove
    rownames(sim) <- level
    tax.sim <- rbind(tax.sim, sim)
}
tax.sim <- t(tax.sim)
tax.sim <- cbind(tax.sim[order(rowMeans(tax.sim), decreasing = T), ], c(0, 1, 0, 0, 0)) # I am attaching a column with 0 and 1 for automatic plot coloring (0 = 0% similarity, 1 = 100% similarity)

# 9.2 Functional similarities ----

# Repeat for functional hierarchy
fun.sim <- data.frame()
for(level in rev(functional)){    
    table <- consolidate(functWGS, level)
    rownames(table) <- table[, level]
    table <- t(table[, -c(1:which(colnames(table) == level))])
    dist <- vegdist(table)
    sim <- similarity(dist)
    sim <- t(sim[ , 'mangrove', drop = F])
    rownames(sim) <- level
    fun.sim <- rbind(fun.sim, sim)
}
fun.sim <- t(fun.sim)
fun.sim <- cbind(fun.sim[order(rowMeans(fun.sim), decreasing = T), ], c(0, 1, 0, 0, 0))

# 9.3 Plotting similarities ----

png("Part 2 - WGS based/tax_sim.png", width = 170, height = 230, res = 250)
par(xpd = TRUE)
heatmap.plus(t(as.matrix(tax.sim)), Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "black"))(100), margins = c(3, 1.75), cexRow = 0.35, cexCol = 0.35, scale = "none")
legend(x = 1.1, y = 0.86, legend = "", cex = 0.5, lty = 0, box.lwd = 0, box.col = "white", bg = "white", text.width = 2) # this line covers the extra column I added with 0/1 values
dev.off()

png("Part 2 - WGS based/fun_sim.png", width = 170, height = 230, res = 250)
par(xpd = TRUE)
heatmap.plus(t(as.matrix(fun.sim)), Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "black"))(100), margins = c(3, 1.75), cexRow = 0.35, cexCol = 0.35, scale = "none")
legend(x = 1.1, y = 0.8, legend = "", cex = 0.5, lty = 0, box.lwd = 0, box.col = "white", bg = "white", text.width = 2)
dev.off()

# Common legend for both plots
png("Part 2 - WGS based/legend_sim.png", width = 68, heigh = 60, res = 250)
par(xpd = TRUE)
par(mar = c(0, 0, 0, 0))
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("topleft", title = "Similarity", legend=c(rep("", 20), "0%", rep("",20), "50%", rep("",20), "100%"), col = c(rep("white", 20), colorRampPalette(c("white", "black"))(43)), lty = 1, lwd = 0.1, box.lwd = 0.5, y.intersp=0.06, cex=0.25)
dev.off()

# 10 PCOA ----

# 10.1 Calculating PCOA ----

factor <- 'biome'
levels <- c('phylum', 'order', 'level1', 'level3')
pcoa_list <- list()
for(level in levels){
    table <- comp[[factor]][[level]]
    rownames(table) <- table[, 1]
    table <- t(table[, -1])
    rownames(table) <- mg_names
    
    dist <- vegdist(table)
    pcoa_list[[level]] <- pcoa(dist)
}

# 10.2 Ploting ----

factors <- factors.for.col[, c('stable organization', 'biome')]
pcoa.colors <- list('stable organization' = colours.org, 'biome' = colours.biome)

for(i in 1:2){
    plot_pcoa(pcoa_list, pcoa.colors, factors, i)
}

# 11 Random Forest ----

factors <- list('biome' = c('phylum', 'level1'), 'stable organization' = c('order', 'level3'))

# 11.1 Calculating random forests ----

tables_for_cda <- list() # The formatted profiles generated below will be also used in the sequence
groups <- list() # And also the vectors containing the groups names
rf_list <- list()
for(factor in names(factors)){
    tables <- list()
    grps <- list()
    rf <- list()
    for(level in factors[[factor]]){
        table <- comp[[factor]][[level]]
        rownames(table) <- table[, 1]
        table <- t(table[, -1])
        rownames(table) <- mg_names
        colnames <- sapply(colnames(table), truncateNames)
        colnames <- dealRepeated(colnames)
        colnames(table) <- colnames
        tables[[level]] <- table
        
        meta <- comp[[factor]][['metadata']]
        group <- meta[, 'groups']
        grps[[level]] <- group
        
        set.seed(1)
        rf[[level]] <- randomForest(table, group, importance=TRUE, proximity=TRUE, ntree=5000)
    }
    tables_for_cda[[factor]] <- tables
    groups[[factor]] <- grps
    rf_list[[factor]] <- rf
}

# 11.2 Plotting ----

imp_names <- plot_rf_importances(rf_list)

# 11.3 Important features ----

# Indexes of the most important features for each level (read from the plot) 
m_imp <- list('biome' = list('phylum' = c(32, 3), 'level1' = c(7, 2)), 'stable organization' = list('order' = c(5, 7), 'level3' = c(10, 7)))

# Do all important features in Gini measure occur within most important features in Accuracy?
for(factor in names(m_imp)){
    for(level in names(m_imp[[factor]])){
        cat(c(factor, '|', level, '\n'))
        cat(sapply(imp_names[[2]][[factor]][[level]][1:m_imp[[factor]][[level]][2]], 
            function(x) x %in% imp_names[[1]][[factor]][[level]][1:m_imp[[factor]][[level]][1]]))
        cat('\n')
        }
}

# Names of the most important features for each level, after pulling accuracy and Gini together
most_important <- list()
for(factor in names(m_imp)){
    most_important_level <- list()
    for(level in names(m_imp[[factor]])){
        most_important_level[[level]] <- unique(c(imp_names[['accs']][[factor]][[level]][1:m_imp[[factor]][[level]][1]], imp_names[['ginis']][[factor]][[level]][1:m_imp[[factor]][[level]][2]]))
        # most_important[[lev]] <- most_important[[lev]][1:min(length(most_important[[lev]]), nrow(df)-1)]
        most_important_level[[level]] <- most_important_level[[level]][1:min(length(most_important_level[[level]]), 14)]
    }
    most_important[[factor]] <- most_important_level
}

# 12 CDA ----

# 12.1 Calculating CDA ----

tables_for_oob_error <- list() # tables_for_cda will be formatted for CDA, I will keep the formatted versions to be used later on
cda_list <- list()

for(factor in names(m_imp)){
    tables <- list()
    cda <- list()
    for(level in names(m_imp[[factor]])){
        table <- as.data.frame(tables_for_cda[[factor]][[level]])
        table$groups <- groups[[factor]][[level]]
        table <- table[, c('groups', most_important[[factor]][[level]])]
        
        table <- formatforcda(table)
        tables[[level]] <- table
        
        model <- lm(as.matrix(table[, 2:ncol(table)]) ~ groups, data = table)
        
        set.seed(1)
        can <- candisc(model, data = df)
        cda[[level]] <- can
        
        print(factor)
        print(level)
        print(can)
        cat('\n')
    }
    tables_for_oob_error[[factor]] <- tables
    cda_list[[factor]] <- cda
}

# 12.2 Ploting ----
    
plot_cda(cda_list)

# 12.3 Estimating OOB error ----

for(factor in names(tables_for_oob_error)){
    for(level in names(tables_for_oob_error[[factor]])){
        table <- tables_for_oob_error[[factor]][[level]]
        # the functions have several limitations to characters in features names. Since I do not need the features names, the easiest way to circumvent this is simply changing features names to a simple mask (I am using LETTERS)
        colnames(table)[2:ncol(table)] <- LETTERS[2:ncol(table)] 
        set.seed(1)
        print(level)
        tryCatch({print(cdaErrorEst.fun(table, 1, 2:ncol(table)))}, error = function(condition){NA})
    }
}

# 13 Pairwise comparisons between mangrove and sea ----

# Plotting pairwise comparisons at all hierarchical levels
factor <- 'mangrove vs sea'
remove_unclassified <- F

level <- 'domain'
tweak_sup <- -2.9
tweak_pval <- 1.45
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'phylum'
tweak_sup <- 0.15
tweak_pval <- 1.55
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'class'
tweak_sup <- 0.2
tweak_pval <- 1.55
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'order'
tweak_sup <- 1
tweak_pval <- 1.6
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'family'
tweak_sup <- 1
tweak_pval <- 1.6
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'genus'
tweak_sup <- 1
tweak_pval <- 1.6
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'level1'
tweak_sup <- 2.9
tweak_pval <- 1.75
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'level2'
tweak_sup <- 2.1
tweak_pval <- 1.65
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'level3'
tweak_sup <- 2.75
tweak_pval <- 1.7
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

level <- 'level4'
tweak_sup <- 2.65
tweak_pval <- 1.7
plot_ttest(results, factor, level, remove_unclassified, tweak_sup, tweak_pval)

# 13.1 Checking enrichment ----

enrich <- list()
factor <- 'mangrove vs sea'
for(level in c('phylum', 'order', 'level1', 'level3')){
    enrich[[level]] <- enriched(results[[factor]][[level]][['ttests']])
}

# Taxa enriched for mangrove compared to sea (max 2%)
for(name in c('phylum', 'order')){
    g <- enrich[[name]][[1]]
    print(name)
    g <- g[abs(g[, 'diff']) >= 2, c('factor', 'diff', 'p')]
    print(names(enrich[[name]]))
    print(g)
    cat('\n')
}

# Functions enriched for mangrove compared to sea (max 0.5%)
for(name in c('level1', 'level3')){
    g <- enrich[[name]][[1]]
    print(name)
    g <- g[abs(g[, 'diff']) >= 0.5, c('factor', 'diff', 'p')]
    print(names(enrich[[name]]))
    print(g)
    cat('\n')
}

save.image('Part 2 - WGS based/WGS.RData')

# FUNCTIONS COMMON TO BOTH WORKFLOWS

# Function to append taxonomic information to those levels described as 'unclassified', using information from higher orders

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

# Function to consolidate the genus count table by the desired taxonomic level, and convert counts to proportions

consolidate <- function(table, level){
    if(level %in% taxonomy){
        table <- table %>% group_by_at(c('domain', level)) %>% summarise_if(is.numeric, sum)
        first_num <- 3
        start <- ifelse(level == 'domain', 2, 3)
    } else{
        table <- table %>% group_by_at(level) %>% summarise_if(is.numeric, sum)
        first_num <- 2
        start <- 2
    }
    table <- table[rowSums(table[, first_num:ncol(table)]) > 0, ]
    table[, start:ncol(table)] <- sapply(table[, start:ncol(table)], function(x){x/sum(x)})
    return(as.data.frame(table))
}

boxplotting <- function(table, level, num_each = NULL){
    title <- unlist(strsplit(level, ''))
    title <- paste0(c(toupper(title[1]), title[2:length(title)]), collapse = '')
    boxplot(table*100, main = title, cex.axis = 1, horizontal = T, las = 2, axes = F, cex.lab = 1)
    axis(1); axis(2, at = 1:ncol(table), side = 4, labels = NA); box()
    if(!is.null(num_each)){
        abline(h = cumsum(num_each)[1:(length(num_each)-1)]+0.5, lty = 3)
    }
    title(xlab = '% sequences', line = 2.5, cex.lab = 1)
    mtext(colnames(table), at = 1:ncol(table), side = 4, line = 1, las = 2)
}

# Function that accepts a proportions table, metadata, and a vector containing the factors to perform T-tests

perform_ttests <- function(table, metadata, factors){
    ttests <- list()
    for(i in 1:(length(factors)-1)){
        for(j in (i+1):length(factors)){
            
            group1 <- metadata[, 'sample'][metadata['groups'] == factors[i]]
            group2 <- metadata[, 'sample'][metadata['groups'] == factors[j]]
            
            table1 <- table[rownames(table) %in% group1, , drop = F]
            table2 <- table[rownames(table) %in% group2, , drop = F]
            
            keepcols <- which(colnames(table) %in% colnames(rbind(table1, table2))[colSums(rbind(table1, table2)) > 0])
            
            table1 <- table1[, keepcols, drop = F]
            table2 <- table2[, keepcols, drop = F]
            
            ttest <- data.frame()
            for(col in colnames(table1)){
                this_ttest <- tryCatch(
                    {t.test(table1[, col], table2[, col], var.equal = F)
                    }, error = function(cond) {list('estimate' = c('mean of x' = mean(table1[, col]),
                                                                   'mean of y' = mean(table2[, col])),
                                                    'conf.int' = NA,
                                                    'p.value' = NA)})
                this_ratio <- tryCatch(
                    {
                        ttestratio(table1[, col], table2[, col], var.equal = F)
                    }, error = function(cond) {list('conf.int' = NA)})
                this_ttest <- data.frame(mean_of_x = this_ttest$estimate['mean of x'], 
                                         mean_of_y = this_ttest$estimate['mean of y'],
                                         diff = this_ttest$estimate['mean of x']-this_ttest$estimate['mean of y'],
                                         diff_lci = this_ttest$conf.int[1],
                                         diff_uci = this_ttest$conf.int[2],
                                         ratio = this_ttest$estimate['mean of x']/this_ttest$estimate['mean of y'],
                                         ratio_lci = this_ratio$conf.int[1],
                                         ratio_uci = this_ratio$conf.int[2],
                                         p = this_ttest$p.value)
                rownames(this_ttest) <- col
                ttest <- rbind(ttest, as.data.frame(this_ttest))
            }
            # ttest$p_adj <- p.adjust(ttest$p, 'bonferroni')
            ttest$sig <- sapply(ttest$p, function(x){if(is.na(x)){NA} else if(x <= 0.001){'***'} else if(x <= 0.01){'**'}  else if(x <= 0.05){'*'} else if(x <= 0.10){'.'} else{''}})
            ttest$sig <- factor(ttest$sig, levels = c('***', '**', '*', '.', ''), ordered = T)
            ttest$factor <- rownames(ttest)
            ttest <- ttest %>% arrange(sig, -abs(diff), -ratio) %>% select_at(c(11, 1:10))    
            
            ttests[[paste0(factors[i], ' vs ', factors[j])]] <- ttest
        }
    }
    return(ttests)
}

# Function to format the T-test results in a plot-friendly format

formatttest <- function(ttest, forPlot = T){
    table <- ttest
    table <- table[!is.na(table$p), ]
    table$sig_dummy <- sapply(table$sig, function(x) ifelse(as.numeric(x) <= 3, '*', ' '))
    table$sig_dummy <- factor(table$sig_dummy, levels = c('*', ' ', ordered = T))
    table <- table %>% arrange(sig_dummy, -abs(diff), -ratio)
    if(forPlot){
        table <- table[1:min(nrow(table), 30), ]
        
        table$factor <- sapply(as.character(table$factor), function(x){
            if(nchar(x) > 50){
                truncated <- paste(c(unlist(strsplit(x, ''))[1:45], '[...]'), collapse = '')
                levels(table$factor) <- c(levels(table$factor), truncated)
                truncated
            }else{
                x
            }})
        repeated <- names(table(table$factor)[table(table$factor) > 1])
        for(rep in repeated){
            which <- which(table$factor == rep)
            i <- 1
            for(w in which){
                table$factor[w] <- sub('(\\[\\.\\.\\.\\])', paste0('\\1\\(', i, '\\)'), table$factor[w])
                i <- i+1
            }
        }
    }
    table$factor <- factor(table$factor, levels = rev(table$factor))
    
    table$p <- sapply(X = table$p, FUN = function(x){
        if(x >= 0.001){sprintf('%1.3f', x)
        } else{scientific(x, digits = 2)}
    })
    table$p <- gsub('([\\+\\-])0(.)', '\\1\\2', table$p)
    
    table$sig <- gsub('^$', ' ', table$sig)
    table$sig <- factor(table$sig, levels = c('***', '**', '*', '.', ' '))
    return(table)
}

# Function to plot the T-test comparisons

plot_ttest <- function(results, factor, level, remove_unclassified = T, tweak_sup = 0, tweak_pval = 1, std = study){
    ttest <- results[[factor]][[level]][['ttests']]
    comparison <- names(ttest)
    ttest <- ttest[[comparison]]
    table <- formatttest(ttest, forPlot = T)
    last <- nrow(table)
    
    colors <- gray.colors(5, start = 0, end = 1)
    pal <- c(
        '***' = colors[1],
        '**' = colors[2],
        '*' = colors[3],
        '.' = colors[4],
        ' ' = colors[5]
    )
    
    ploting <- unlist(strsplit(comparison, ' '))[c(1, 3)]
    
    limits <- c(table$ratio, table$ratio_lci, table$ratio_uci)
    limits <- limits[(abs(limits) != Inf) & (limits != 0)]
    limits <- abs(limits[!is.na(limits)])
    if(length(limits) == 0){
        limits <- 10
    }
    limit <- max(c(limits, 1/min(limits)))*1.01
    decimals <- ceiling(log10(limit))
    decimals <- ceiling(decimals/2)
    brk <- c(1*10^-decimals, 1, 1*10^decimals)
    if(brk[3] > limit){limit <- brk[3]}
    
    g1 <- ggplot(table, aes(x = factor, y = ratio, fill = sig)) +
        scale_y_log10(limits = c(1/limit, limit), breaks = brk, labels = function(x){format(x, scientific = F, drop0trailing = T)}) +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
        geom_errorbar(width = 1, aes(ymin = ratio_lci, ymax = ratio_uci)) +
        geom_point(shape = 21, size = 3, show.legend = F) +
        scale_fill_manual(values = pal) +
        theme(plot.title = element_text(hjust = 0.5, size = 9, margin = ggplot2::margin(2, 0, 2, 0)),
              plot.subtitle = element_text(hjust = 0.5, size = 8, face = 'italic'),
              plot.margin = unit(c(0, (-3.15-tweak_sup), 0, 0), 'lines'),
              plot.background = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()
        ) +
        ggtitle('Ratio (log scale)', subtitle = paste0(ploting[2], ' | ', ploting[1])) +
        coord_flip()
    
    whiches <- unique(c(which(table$ratio == Inf), which(is.na(table$ratio_lci)), which(table$ratio_lci < 0), which(table$ratio_uci < 0)))
    for(which in whiches){
        g1 <- g1 +
            geom_point(data = table[which, ], aes(x = factor, y = ratio), shape = 21, size = 3, color = 'red', show.legend = F)
    }
    
    limits <- c(table$diff, table$diff_lci, table$diff_uci)
    limits <- limits[abs(limits) != Inf]
    limits <- abs(limits[!is.na(limits)])
    limit <- max(c(limits, -min(limits)))
    decimals <- ceiling(-log10(limit))
    brk0 <- ceiling(limit*10^decimals)/10^decimals
    brk <- c(-brk0/2, 0, brk0/2)
    
    g2 <- ggplot(table, aes(x = factor, y = diff, fill = sig)) +
        scale_x_discrete(position = 'top') +
        scale_y_continuous(breaks = brk, labels = scales::percent_format(accuracy = 0.001, drop0trailing = T)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
        geom_errorbar(width = 1, aes(ymin = diff_lci, ymax = diff_uci)) +
        geom_point(shape = 21, size = 3) +
        scale_fill_manual(values = pal) +
        geom_text(aes(x = last+1.1, label = 'p-value', y = limit*tweak_pval), hjust = 1, size = 2.5, color = 'darkgray') +
        geom_text(aes(x = factor, label = p, y = limit*tweak_pval), hjust = 1, size = 3, color = 'darkgray') +
        theme(plot.title = element_text(hjust = 0.5, size = 9, margin = ggplot2::margin(2, 0, 2, 0)),
              plot.subtitle = element_text(hjust = 0.5, size = 8, face = 'italic'),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              plot.margin = unit(c(0, 2, 0, (-5.65-tweak_sup)), 'lines'),
              plot.background = element_blank(),
              legend.title = element_text(size = 6),
              legend.text = element_text(size = 6),
              legend.key.size = unit(0.5, 'lines'),
              legend.background = element_rect(fill = alpha('white', 0.5)),
              legend.position = c(1, 0),
              legend.justification = c(1, 0)
        ) +
        ggtitle('Difference', subtitle = paste0(ploting[2], ' | ', ploting[1])) +
        coord_flip(xlim = 1:last, ylim = c(-limit, limit), clip = 'off')
    
    gt1 <- ggplotGrob(g1)
    gt2 <- ggplotGrob(g2)
    
    newWidth <- unit.pmax(gt1$widths[2:4], gt2$widths[2:4])
    
    gt1$widths[2:4] = as.list(newWidth)
    gt2$widths[2:4] = as.list(newWidth)
    
    folder <- ifelse(std == '16S', 'Part 1 - 16S rRNA gene based/', 'Part 2 - WGS based/')
    
    title <- unlist(strsplit(level, ''))
    title <- paste0(c(toupper(title[1]), title[2:length(title)]), collapse = '')
    png(paste0(folder, 'ttest_', level, '_', gsub(' ', '_', comparison), '.png'), width = 900, height = 900, res = 200)
    grid.arrange(gt1, gt2, ncol = 2, top = title)
    dev.off()
}

# Function to check enrichment between two groups

enriched <- function(ttests){
    enriched <- list()
    names <- names(ttests)
    for(name in names){
        ttest <- ttests[[name]]
        table <- formatttest(ttest, F)
        table <- table[as.numeric(table$sig) <= 3, c('factor', 'diff', 'ratio', 'p')]
        table$diff <- table$diff*100
        enriched[[name]] <- table
    }
    return(enriched)
}

# Function to simplify a 'prepared for hclust' table, making it to contain only the 20 most abundant taxa, collapsing all others in a __minor__ category

simplify <- function(table){
    if(ncol(table) >= 20){
        table <- cbind(rowSums(table[, 21:dim(table)[1]]), table[, 20:1])
        colnames <- colnames(table)
        colnames[1] <- '__minor__'
    }else{
        table <- table[, ncol(table):1]
        colnames <- colnames(table)
    }
    for(i in 1:length(colnames)){
        if(nchar(colnames[i]) > 30){
            colnames[i] <- paste0(substr(colnames[i], 1, 24), ' [...]')
        }
    }
    repeated <- names(table(colnames)[table(colnames) > 1])
    for(rep in repeated){
        which <- which(colnames == rep)
        i <- 1
        for(w in which){
            colnames[w] <- sub('...(\\[\\.\\.\\.\\])', paste0('\\1\\(', i, '\\)'), colnames[w])
            i <- i+1
        }
    }
    colnames(table) <- colnames
    return(table)
}

# Function to plot a heatmap of the desired simplified table, from the desired hierarchical level (which needs to be discriminated for naming output file)

plot_heatmap <- function(table, level, std = study){
    folder <- ifelse(std == '16S', 'Part 1 - 16S rRNA gene based/', 'Part 2 - WGS based/')
    png(paste0(folder, 'heatmap_', level, '.png'), width = 850, height = 650, res = 200)
    par(xpd = TRUE)
    
    title <- unlist(strsplit(level, ''))
    title <- paste0(c(toupper(title[1]), title[2:length(title)]), collapse = '')
    
    heatmap.plus(t(table), main = title, margins = c(3.3, 10), Rowv = NA, Colv = as.dendrogram(clust), ColSideColors = as.matrix(as.data.frame(ann.colors)), cexRow = 0.8, cexCol = 0.8, col=colorRampPalette(c('blue', 'yellow', 'red'))(100), scale= 'none')
    
    abundance_legend <- c('0%',
                          rep('',20),
                          paste0(signif(mean(c(min(table), max(table))), 0)*100, '%'),
                          rep('', 20),
                          paste0(signif(max(table), 0)*100, '%'))
    legend(x = 0.73, y = -0.45, title = 'Abundance', title.adj = 0.95, legend = abundance_legend, col = colorRampPalette(c('blue', 'yellow', 'red'))(43), lty = 1, lwd = 0.5, y.intersp = 0.05, cex = 0.45, text.width = 0.17)
    
    text(x = 0.72, y = 1.05, cex = 0.53, expression(italic('Environmental levels')), adj = 0)
    
    envir_names = names(ann.colors)
    if(study == '16S'){
        envir_x_pos = c(0.72, 0.82)
        envir_y_pos = c(1.35, 1.35)
    }else{
        envir_x_pos = c(0.72, 0.72, 0.92)
        envir_y_pos = c(1.35, 1.225, 1.35)
    }
    for(i in 1:length(envir_names)){
        legend(x = envir_x_pos[i], y = envir_y_pos[i], title = envir_names[i], title.adj = 0, legend = levels(factors.for.col[, tolower(envir_names[i])]), col = unique(ann.colors[[envir_names[i]]]), lty = 1, lwd = 5, seg.len = 1, bty = 'n', y.intersp = 0.72, cex = 0.35)
    }
    
    dev.off()
}

# Function to plot two chosen dimensions of a PCOA object

plot_pcoa <- function(pcoa_list, colors_list, metadata, first_dimension, second_dimension = first_dimension + 1, std = study){
    folder <- ifelse(std == '16S', 'Part 1 - 16S rRNA gene based/', 'Part 2 - WGS based/')
    
    for(n in 1:length(pcoa_list)){
        level <- names(pcoa_list)[n]
        pcoa <- pcoa_list[[n]]
        
        png(paste0(folder, 'pcoa_', level, '_', first_dimension, 'vs', second_dimension, '.png'), width = 960, height = 960, res = 200)
        
        var_exp <- round(pcoa$values$Relative_eig*100, 1)
        points <- pcoa$vectors
        labels <- dimnames(pcoa$vectors)[[1]]
        
        limits <- c(min(points[, first_dimension]), max(points[, first_dimension]), min(points[, second_dimension]), max(points[, second_dimension]))
        ranges <- c(limits[2]-limits[1], limits[4]-limits[3])
        xlim <- c(limits[1]-ranges[1]*0.3, limits[2]+ranges[1]*0.1)
        ylim <- c(limits[3]-ranges[2]*0.3, limits[4]+ranges[2]*0.1)
        
        names = colnames(metadata)
        
        par(mar = c(4, 4, 0.5, 0.5))
        plot(points[, c(first_dimension, second_dimension)], 
             xlim = xlim,
             ylim = ylim,
             xlab = paste0('Axis ', first_dimension, ' (', var_exp[first_dimension], '%)'),
             ylab = paste0('Axis ', second_dimension, ' (', var_exp[second_dimension], '%)'),
             col = colors_list[[names[1]]],
             bg = colors_list[[names[2]]],
             pch = 22,
             lwd = 1.5,
             cex = 1
        )
        text(points[, i]-0.03*(xlim[2]-xlim[1]), points[, i+1]-0.03*(ylim[2]-ylim[1]), labels, cex = 0.7, col = adjustcolor('black', alpha.f = 0.5))
        
        
        legend('bottomleft',
               title = names[1],
               legend = c(unique(as.character(metadata[, names[1]]))), 
               pt.bg = 'white',
               col = c(unique(colors_list[[names[1]]])),
               pch = 22,
               cex = 0.7
        )
        legend('bottomright',
               title = names[2],
               legend = c(unique(as.character(metadata[, names[2]]))), 
               pt.bg = c(unique(colors_list[[names[2]]])),
               col = 'white',
               pch = 22,
               cex = 0.7
        )
        dev.off()
    }
}

# Function to plot Random Forest importances (Accuracy and Gini index). It will also return a list containing the most important features by Accuracy and by Gini.

plot_rf_importances <- function(rf_list, std = study){
    folder <- ifelse(std == '16S', 'Part 1 - 16S rRNA gene based/', 'Part 2 - WGS based/')
    n_features <- ifelse(std == '16S', 30, 50)
    cex <- ifelse(std == '16S', 0.7, 0.5)
    factors <- names(rf_list)
    accs <- list()
    ginis <- list()
    for(factor in factors){
        levels <- names(rf_list[[factor]])
        acc0 <- list()
        gini0 <- list()
        for(level in levels){
            rf <- rf_list[[factor]][[level]]
            png(paste0(folder, 'rf_', factor, '_', level, '.png'), height = 800, width = 1500, res = 200)
            par(mfrow = c(1, 2), mar = c(4.1, 0, 0.5, 1))
            
            # Ploting MeanDecreaseAccuracy
            acc <- varImpPlot(rf, n.var = min(n_features, nrow(rf$importance)), type = 1, pch = NA, cex = cex, main='')
            
            # Defining pch according to the group to which a given phylum is more important accuracy-wise
            min <- min(n_features, nrow(acc))
            acc_names <- rownames(acc)[rev(order(acc))][1:min]
            imp <- rf$importance[acc_names, ]
            imp <- imp[rev(order(imp[, 'MeanDecreaseAccuracy'])), c(1:length(levels(rf$y)))]
            imp <- factor(apply(imp, 1, function(x) names(which(x == max(x)))))
            pch <- c(15, 16, 17, 18, 4)[imp]
            
            # Adding points
            points(rev(sort(acc))[1:min], min:1, pch = pch)
            
            legend('bottomright', title = 'Groups', legend = unique(imp), pch = unique(pch), cex = 0.55)
            
            # Ploting MeanDecreaseGini
            gini <- varImpPlot(rf, n.var = min(n_features, nrow(rf$importance)), type = 2, pch = NA, cex = cex, main='')
            
            # Defining pch according to the group to which a given phylum is more important Gini-wise
            min <- min(n_features, nrow(gini))
            gini_names <- rownames(gini)[rev(order(gini))][1:min]
            imp <- rf$importance[gini_names, ]
            imp <- imp[rev(order(imp[, 'MeanDecreaseGini'])), c(1:length(levels(rf$y)))]
            imp <- factor(apply(imp, 1, function(x) names(which(x == max(x)))))
            pch <- c(15, 16, 17, 18, 4)[imp]
            
            # Adding points
            points(rev(sort(gini))[1:min], min:1, pch = pch)
            
            legend('bottomright', title = 'Groups', legend = unique(imp), pch = unique(pch), cex = 0.55)
            
            dev.off()
            
            acc0[[level]] <- acc_names
            gini0[[level]] <- gini_names
        }
    accs[[factor]] <- acc0
    ginis[[factor]] <- gini0
    }
    return(list('accs' = accs, 'ginis' = ginis))
}

# Function to remove features with no intragroup standard deviation

formatforcda <- function(table){
    sds <- table %>% group_by(groups) %>% summarise_if(is.numeric, sd)
    sds <- sds[, 2:ncol(sds)]
    zerosd <- apply(sds, 2, function(x) any(round(x, 5) == 0))
    zerosd <- colnames(sds)[zerosd]
    if(length(zerosd) > 0){
        table <- table[, -sapply(zerosd, function(x) which(colnames(table) == x))]    
    }
    return(as.data.frame(table))
}

# This function is a slight modification from the standard plot function applied to a candisc object, since the original function was not plotting colors correctly

plot_candisc <- function(x, which = 1:2, conf = 0.95, col, pch, scale, asp = 1, 
                         var.col = "blue", var.lwd = par("lwd"), var.labels, var.cex = 1, 
                         var.pos, rev.axes = c(FALSE, FALSE), ellipse = FALSE, ellipse.prob = 0.68, 
                         fill.alpha = 0.1, prefix = "Can", suffix = TRUE, titles.1d = c("Canonical scores", 
                                                                                        "Structure"), ...) 
{
    term <- x$term
    factors <- x$factors
    rev.axes <- rep(rev.axes, length.out = 2)
    nlev <- if (is.matrix(x$means)) 
        nrow(x$means)
    else length(x$means)
    if (missing(col)) 
        col <- rep(palette()[-1], length.out = nlev)
    fill.col <- heplots::trans.colors(col, fill.alpha)
    if (x$ndim < 2 || length(which) == 1) {
        which <- which[1]
        op <- par(no.readonly = TRUE)
        ng <- length(x$means)
        structure <- as.vector(x$structure[, which])
        canvar <- paste("Can", which, sep = "")
        scores <- x$scores[, canvar]
        if (isTRUE(rev.axes[1])) {
            scores <- -scores
            structure <- -structure
        }
        ns <- length(structure)
        wid <- if (ns < 2 * ng) 
            c(2, 1)
        else c(1.2, 1)
        layout(matrix(c(1, 2), 1, 2), widths = wid)
        par(mar = c(5, 4, 4, 0) + 0.1)
        if (is.logical(suffix) & suffix) 
            suffix <- paste(" (", round(x$pct[which], 1), "%)", 
                            sep = "")
        else suffix <- NULL
        canlab <- paste(prefix, which, suffix, sep = "")
        formule <- formula(paste(canvar, " ~", term, sep = ""))
        boxplot(formule, data = x$scores, ylab = canlab, xlab = term, 
                col = fill.col, main = titles.1d[1], ...)
        xx <- 1:ns
        par(mar = c(5, 0, 4, 1) + 0.1)
        ylim <- range(c(structure, if (all(structure > 0)) +1 else 0, 
                        if (all(structure < 0)) -1 else 0))
        xlim <- c(0.5, ns + 0.5)
        plot(xx, structure, type = "n", ylab = "", xlab = "", 
             xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
             main = titles.1d[2])
        arrows(xx, 0, xx, structure, length = 0.1, angle = 15, 
               col = var.col, lwd = var.lwd)
        abline(h = 0, lty = 2, col = "gray")
        vars <- if (missing(var.labels)) 
            rownames(x$structure)
        else var.labels
        adj1 <- as.vector(ifelse(structure > 0, 1, 0))
        adj2 <- rep(-0.3, ns)
        adj2[1] <- 1.1
        for (i in 1:ns) text(xx[i], structure[i], paste("  ", 
                                                        vars[i], "  "), adj = c(adj1[i], adj2[i]), col = var.col, 
                             srt = 90, cex = var.cex, xpd = TRUE)
        par(op)
        return(invisible())
    }
    canvar <- paste("Can", which, sep = "")
    if (is.logical(suffix) & suffix) 
        suffix <- paste(" (", round(x$pct[which], 1), "%)", 
                        sep = "")
    else suffix <- NULL
    canlab <- paste(prefix, which, suffix, sep = "")
    nlev <- nrow(x$means)
    if (missing(col)) 
        col <- rep(palette(), length.out = nlev)
    if (missing(pch)) 
        pch <- rep(1:18, length.out = nlev)
    scores <- x$scores[, canvar]
    means <- x$means[, which]
    labels <- rownames(x$means)
    structure <- x$structure[, which]
    if (isTRUE(rev.axes[1])) {
        scores[, 1] <- -scores[, 1]
        means[, 1] <- -means[, 1]
        structure[, 1] <- -structure[, 1]
    }
    if (isTRUE(rev.axes[2])) {
        scores[, 2] <- -scores[, 2]
        means[, 2] <- -means[, 2]
        structure[, 2] <- -structure[, 2]
    }
    Ind <- dataIndex(x$scores, term)
    plot(scores, asp = asp, xlab = canlab[1], ylab = canlab[2], 
         col = col[Ind], pch = pch[Ind], ...)
    abline(h = 0, v = 0, lty = 2, col = "grey")
    if (ellipse) {
        fill.col <- heplots::trans.colors(col, fill.alpha)
        radius <- sqrt(qchisq(ellipse.prob, df = 2))
        angles <- (0:60) * 2 * pi/60
        circle <- radius * cbind(cos(angles), sin(angles))
        for (i in 1:nlev) {
            sigma <- var(scores[Ind == i, ])
            mu <- as.numeric(means[i, ])
            Q <- chol(sigma, pivot = TRUE)
            order <- order(attr(Q, "pivot"))
            ell <- sweep(circle %*% Q[, order], 2, mu, FUN = "+")
            polygon(ell, col = fill.col[i], border = col[i], 
                    lty = 1)
        }
    }
    ord_factors <- unique(factors[, 1])
    means1 <- means[ord_factors, ]
    points(means1[, 1], means1[, 2], col = col, pch = "+", cex = 2)
    pos <- ifelse(means[, 2] > 0, 3, 1)
    text(means[, 1], means[, 2], labels = labels, pos = pos)
    if (missing(scale)) {
        scale <- vecscale(structure)
        message("Vector scale factor set to ", round(scale, 
                                                     3))
    }
    cs <- scale * structure
    if (!missing(var.labels)) 
        rownames(cs) <- var.labels
    vectors(cs, col = var.col, cex = var.cex, lwd = var.lwd, 
            pos = var.pos, xpd = TRUE)
    circle <- function(center, radius, segments = 41, ...) {
        angles <- (0:segments) * 2 * pi/segments
        unit.circle <- cbind(cos(angles), sin(angles))
        circle <- t(as.numeric(center) + radius * t(unit.circle))
        lines(circle, col = col, ...)
    }
    if (conf > 0) {
        n <- as.vector(table(factors))
        radii <- sqrt(qchisq(conf, 2)/n)
        symbols(means, circles = radii, inches = FALSE, add = TRUE, 
                fg = col)
    }
}

# Function to plot CDA objects stored in a list

plot_cda <- function(cda_list, std = study){
    folder <- ifelse(std == '16S', 'Part 1 - 16S rRNA gene based/', 'Part 2 - WGS based/')
    factors <- names(cda_list)
    for(factor in factors){
        levels <- names(cda_list[[factor]])
        for(level in levels){
            can <- cda_list[[factor]][[level]]
            
            png(paste0(folder, 'cda_', factor, '_', level, '.png'), height = 1750, width = 1750, res = 200)
    
            par(mar = c(4, 4, 1, 1))
            
            col <- if(factor == 'biome') {unique(colours.biome)} else if(factor == 'stable organization') {unique(colours.org)} else {rainbow(can$dfh+1)}
            
            xlims = c(min(can$scores$Can1), max(can$scores$Can1))
            ylims = c(min(can$scores$Can2), max(can$scores$Can2))
            lims = c(min(xlims[1], ylims[1]), max(xlims[2], ylims[2]))
            
            plot_candisc(can, xlim = xlims*1.1, ylim = ylims*1.1, pch = rep(21, can$dfh+1), col = col, scale = 0, conf = 0, cex = 1.2, cex.lab = 1, cex.axis = 1, var.labels = NULL, var.col = 'darkgreen', prefix = 'Canonical Discriminant ')
            
            arrows <- can$structure
            
            mult <- c()
            for(row in 1:nrow(arrows)){
                arrow <- arrows[row, , drop = F]
                mult_x <- abs(arrow[, 1] - lims)
                mult_x <- ifelse(arrow[, 1] < 0, mult_x[1], mult_x[2])
                mult_x <- mult_x/abs(arrow[, 1])
                mult_y <- abs(arrow[, 2] - lims)
                mult_y <- ifelse(arrow[, 2] < 0, mult_y[1], mult_y[2])
                mult <- min(mult, mult_x, mult_y)
            }
            
            for(row in 1:nrow(arrows)){
                arrows(0, 0, arrows[row, 1]*mult, arrows[row, 2]*mult, length = 0.1, col = 'darkgreen')
            }
            
            mult <- mult*1.1
            
            if((factor == 'stable organization') & (level == 'order')){
                condition <- arrows[, 2] < 0
                y <- -30
            }else if((factor == 'stable organization') & (level == 'level3')){
                condition <- arrows[, 1] > 0
                y <- -38
            }else if((factor == 'mangroups') & (level == 'class')){
                condition <- arrows[, 1] >= arrows[8, 1]
                y = -10
            }else{
                condition <- rep(F, nrow(arrows))
            }
            
            dislocate <- data.frame()
            for(row in 1:nrow(arrows)){
                arrow <- arrows[row, , drop = F]
                if(condition[row]){
                    dislocate <- rbind(dislocate, arrow)
                    arrow <- NULL
                }else{
                    if(arrow[, 1] <= 0){
                        text(arrow[, 1]*mult, arrow[, 2]*mult, rownames(arrow), cex = 0.9, col = 'darkgreen', adj = 0.5)
                    }else{
                        text(arrow[, 1]*mult, arrow[, 2]*mult, rownames(arrow), cex = 0.9, col = 'darkgreen', adj = 0.5)
                    }
                }
                
            }
            if(nrow(dislocate > 0)){
                dislocate <- dislocate[rev(order(dislocate$Can2)), ]
                text(x = xlims[2]*1.1, y = y, labels = paste0(rownames(dislocate), '\n', collapse = ''), cex = 0.9, col = 'darkgreen', adj = c(1, 0))
            }
            
            dev.off()
        }
    }
}

# The next four functions were obtained from Dinsdale et al. (2013). They estimate the OOB error based on a CDA model.

# define a function to get the out of bag error
bagData.fun <- function(data, outPercent, response=c(1)){
    # Get information about the classes
    envColumn <- getSingleResponse.fun(data,response)
    vCategories <- levels(data[,envColumn])
    num_levels <- length(vCategories)
    
    # For each class, take out 20% as out of bag.
    for (i in 1:num_levels){
        # Pick out just the data for that class
        subData <- data[grep(vCategories[i],data[,envColumn]),] 
        
        # Check if there are more rows than sample size
        sampleSize <- round(outPercent * length(subData[,1]))
        randValues <- sample(1:length(subData[,1]), sampleSize) # Randomly select
        bData <- subData[-1 * randValues,] 
        outData <- subData[randValues,]
        if (i == 1){ # First time through, create new data frame
            bagData <- bData
            oobData <- outData
        }
        else{ # Subsequent times through, add rows to new data frame
            bagData <- rbind(bagData, bData)
            oobData <- rbind(oobData, outData)
        }
    }
    return(list(bag=bagData, oob=oobData))
}

createFormula.fun <- function(y, response, predictors){
    ifelse(is.character(predictors[1]),for(i in 1:length(predictors)){(1:length(names(y)))[names(y)==paste(predictors[i],sep="")]->predictors[i]},predictors<-predictors);
    predictors.names<-paste(names(y)[predictors], sep="");
    
    ifelse(is.character(response),(1:length(names(y)))[names(y)==paste(response,sep="")]->response,response<-response);
    response.name<-paste(names(y)[response], sep="");
    fmla<-as.formula(paste("cbind(",paste(predictors.names,collapse=","),")"," ~ ",paste(response.name,collapse=""),collapse=""));
    return(fmla)
}

getSingleResponse.fun <- function(y, response){
    if (is.character(response)){
        for (i in 1:length(names(y))){
            if (names(y)[i] == response) response = i
        }
    }
    return(response)
}

# estimate the CDA error
cdaErrorEst.fun <- function(metagenome, response, predictors, outPercent=.2, trials=length(metagenome[,1])){
    totalTest = 0
    totalError = 0
    fmla <- createFormula.fun(metagenome,response,predictors)
    
    # Initialize confusion matrix
    envColumn <- getSingleResponse.fun(metagenome,response)
    envFactor <- levels(metagenome[,envColumn])
    facLength <- length(envFactor)
    confuMatrix <-  matrix(rep(0, facLength * (facLength+1)),nrow=facLength,ncol=(facLength+1),dimnames=list(envFactor,c(envFactor,"class.error")))
    for (i in 1:trials){
        # Get in-bag and out-of-bag sets
        metagenomeBags <- bagData.fun(metagenome, outPercent, response)
        bagMG <- metagenomeBags$bag
        oobMG <- metagenomeBags$oob
        
        # Do CDA on bag data
        metagenome.mod <- lm(fmla, data=bagMG)
        metagenome.can <- candisc(metagenome.mod, data=bagMG)
        can.scores <- metagenome.can$scores
        can.env <- can.scores[,1]
        
        # Do LDA on that CDA
        metagenome.lda <- lda(can.scores[,-1], can.env) 
        
        # Get CDA scores on out-of-bag data
        bag.means <- apply(bagMG[,predictors],2,mean)
        oobMG.norm <- oobMG[,predictors]
        for (j in 1:length(oobMG[,1])){  # Recenter each row
            oobMG.norm[j,] <- oobMG.norm[j,] - bag.means
        }
        oobMG.scores <- as.matrix(oobMG.norm) %*% as.matrix(metagenome.can$coeffs.raw)
        
        # Predict using those scores and compare to actual environment
        metagenome.predict <- predict(metagenome.lda, oobMG.scores)
        for (k in 1:length(oobMG[,1])){
            totalTest = totalTest + 1
            currentRow <- oobMG[k,]
            trueClass <- as.character(oobMG$groups[k])
            predictClass <- as.character(metagenome.predict$class[k])
            confuMatrix[trueClass,predictClass] <- confuMatrix[trueClass,predictClass] + 1
            if (predictClass != trueClass){
                totalError = totalError + 1
            }
        }
    }
    confuMatrix <- confuMatrix / trials
    for (m in 1:facLength) { confuMatrix[m,facLength+1] = 1 - (confuMatrix[m,m] / sum(confuMatrix[m,])) }
    return (list(error=(totalError / totalTest), confusion=confuMatrix))
}
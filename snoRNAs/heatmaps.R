#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)
library(pheatmap)
library(foreach)

o <- foreach(filename=list.files(pattern='*.csv'), .combine='rbind') %do% {
    dat <- fread(filename)
    filename.split <- unlist(strsplit(filename, split='_|\\.'))
    grp <- filename.split[2]
    treatment <- filename.split[3]
    dat <- melt(dat, measure.vars=colnames(dat)[-c(1)], variable.name='sample')
    dat[, 'grp' := grp]
    dat[, 'treatment' := treatment]
    setnames(dat, 'V1', 'gene')
    return(dat[])
}


o[sample %like% 'Neurite', zone := 'Neurite']
o[sample %like% 'Soma', zone := 'Soma']
o[grp == 'tdp' & sample %like% paste0('_S',1:5, '$', collapse='|'), zone := 'Neurite']
o[grp == 'tdp' & sample %like% paste0('_S',6:10, '$', collapse='|'), zone := 'Soma']
o[grp == 'tdp' & sample %like% paste0('_S',11:15, '$', collapse='|'), zone := 'Neurite']
o[grp == 'tdp' & sample %like% paste0('_S',16:20, '$', collapse='|'), zone := 'Soma']

o[, sample := tstrsplit(sample, split='_S')[2]]
o[, sample := as.numeric(sample)]
o[, sample := factor(sample)]

# Ensure all rows have Neurite or Soma labeled zone
stopifnot(nrow(o[is.na(zone)]) == 0)

o[, 'log10Value' := log10(value+1)]

gene_superset <- rev(sort(unique(o$gene)))
o[, gene := factor(gene, levels=gene_superset)]



# Prettify labels
o[grp=='fus', grp := 'FUS']
o[grp=='tdp', grp := 'TDP-43']
o[grp=='hnrnpa1', grp := 'hnRNPA1']
o[treatment=='ctrl', treatment := 'Non-targeting control']
o[treatment=='kd', treatment := 'Knockdown']

o[, facet_grp := paste0(zone, ' ', treatment)]

# Add in combinations that don't exist, give NA values, so levels aren't dropped
all_levels <- CJ('gene'=sort(unique(o$gene)),
    'treatment'=sort(unique(o$treatment)),
    'zone'=sort(unique(o$zone)),
    'grp'=sort(unique(o$grp)))

setkey(all_levels, gene, treatment, zone, grp)
setkey(o, gene, treatment, zone, grp)

o <- merge(o, all_levels, all=T)



# Plot log scale
for(i in unique(o$grp)) {
    g <- ggplot(o[grp == i], aes(x=sample, y=gene, fill=log10Value)) +
            geom_tile() +
            scale_fill_viridis(limits=c(0,max(o$log10Value)), option='magma') +
            facet_grid(~facet_grp, scales='free_x') +
            theme_few() +
            labs(y='snoRNA', title=i)
    ggsave(g, file=paste0(i, '_snoRNAs_log10.pdf'), width=35, height=45, units='cm')
}

# Plot linear scale
for(i in unique(o$grp)) {
    g <- ggplot(o[grp == i], aes(x=sample, y=gene, fill=value)) +
            geom_tile() +
            scale_fill_viridis(limits=c(0,max(o$value)), option='magma') +
            facet_grid(~facet_grp, scales='free_x') +
            theme_few() +
            labs(y='snoRNA', title=i)
    ggsave(g, file=paste0(i, '_snoRNAs_linear.pdf'), width=35, height=45, units='cm')
}



## Nick's WT Analysis
library(data.table)
fn <- 'WT_data.tsv'
dat <- fread(fn)
setnames(dat, 'hgnc_symbol', 'gene')
dat[, treatment := 'WT']

dat[sample %like% 'WT_N', zone := 'Neurite']
dat[sample %like% 'WT_S', zone := 'Soma']
dat <- dat[gene %in% levels(o$gene)]
dat[, gene := factor(gene, levels=levels(o$gene))]
dat[, log10Value := log10(value+1)]

# Plot log scale
g <- ggplot(dat, aes(x=sample, y=gene, fill=log10Value)) +
        geom_tile() +
        scale_fill_viridis(option='magma') +
        facet_grid(~zone, scales='free_x') +
        theme_few() +
        labs(y='snoRNA')

ggsave(g, file='WT_snoRNAs_log10.pdf', width=20, height=45, units='cm')


# Plot linear scale
g <- ggplot(dat, aes(x=sample, y=gene, fill=value)) +
        geom_tile() +
        scale_fill_viridis(limits=c(0,15), option='magma') +
        facet_grid(~zone, scales='free_x') +
        theme_few() +
        labs(y='snoRNA')

ggsave(g, file='WT_snoRNAs_linear.pdf', width=20, height=45, units='cm')
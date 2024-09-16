#!/usr/bin/env Rscript

# R 4.3.2

library(data.table)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(viridis)
library(ggthemes)

####################################################################################################
enrich_pvalue <- 0.05

getGO <- function(entrez_list, all_gene_list, enrich_pvalue=0.05) {
    desired_cols <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count','group')
    d1 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d1)) {
        d1 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d1, desired_cols)
    } else {
        d1 <- as.data.table(d1@result)
    }
    d1[, 'group' := 'go_cc']

    d2 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d2)) {
        d2 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d2, desired_cols)
    } else {
        d2 <- as.data.table(d2@result)
    }
    d2[, 'group' := 'go_bp']

    d3 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "MF",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d3)) {
        d3 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d3, desired_cols)
    } else {
        d3 <- as.data.table(d3@result)
    }
    d3[, 'group' := 'go_mf']

    d4 <- enrichKEGG(gene         = entrez_list,
                                    organism     = 'hsa',
                                    universe      = all_gene_list,
                                    pAdjustMethod = "BH",
                                    qvalueCutoff  = enrich_pvalue,
                                    pvalueCutoff = enrich_pvalue)
    if(is.null(d4)) {
        d4 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d4, desired_cols)
    } else {
        d4 <- as.data.table(d4@result)
    }
    d4[, 'group' := 'kegg']
    d4 <- d4[, .SD, .SDcols=desired_cols]
    rbindlist(list(d1,d2,d3,d4))
}

plot_enrich <- function(dat, filestem, lims=NULL, x.lbl=NULL) {
        # Concatenate into one table
        dat <-dat[qvalue <= 0.05]
        dat[group=='go_mf', 'Domain' := 'GO Molecular Function']
        dat[group=='go_cc', 'Domain' := 'GO Cellular Component']
        dat[group=='go_bp', 'Domain' := 'GO Biological Process']
        dat[group=='kegg', 'Domain' := 'KEGG Pathway']

        dat.hold <- copy(dat)

        for (grp in c('go_mf','go_cc','go_bp','kegg')) {
        dat <- copy(dat.hold)
        dat <- dat[group==grp]
        # continue with next if no signif
        if (nrow(dat) == 0) { next }

        # Create qsort column to sort within 'counts' ties
        dat[, 'qsort' := -1*qvalue]

        # Set order for plotting
        setkey(dat, 'Domain', Count, qsort)
        #dat[, Description := stringr::str_wrap(Description, width=35)]
        desc_order <- dat$Description
        dat[Count > 0, 'hj' := 1]
        dat[Count > 0, Description := paste0(Description, ' ')]
        dat[Count < 0, 'hj' := 0]
        dat[Count < 0, Description := paste0(' ', Description)]
        dat[, Description := factor(Description, levels=desc_order)]


        g <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
                theme_few() +
                scale_fill_viridis(limits=c(0,0.05), direction=-1) +
                geom_vline(xintercept=0) +
                labs(x=x.lbl) +
                labs(y='GO') +
                #ggplot2::facet_grid(`GO domain` ~ ., drop=T, scales='free', space='free', switch='y') +
                theme(strip.text.y.left = element_text(angle = 0)) +
                geom_text(data=dat, aes(x=0, label=Description, hjust=hj)) +
                theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
                scale_x_continuous(expand = expansion(mult = 0.1)) +
                scale_x_continuous(breaks=scales::pretty_breaks(), limits=lims) +
                facet_grid(.~Domain)

        height_cm <- 5 + 0.5 * nrow(dat)
        outfile <- paste0(filestem, '-', grp, '.pdf')
        ggsave(g, file=outfile, width=40, height=height_cm, units='cm', limitsize = FALSE)
    }
}


symbol_to_entrez <- function(symbols) {
    entrez <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID",keytype="SYMBOL")
    return(unique(unlist(entrez)))
}

plot_MA <- function(DT, color_values) {
        g <- ggplot(DT, aes(x=AveExpr, y=logFC, color=Group)) +
                geom_point(data=DT[Group=='NS']) +    # Plot insigificant first
                geom_point(data=DT[Group !='NS']) +   # Plot remaining data on layer ABOVE NS
                scale_color_manual(values=color_vals) +     # Force manual colors in `color_values`
                geom_hline(yintercept=c(-1,1), linetype='dashed') + # foldchange threshold
                theme_minimal(12)
        return(g)
}

plot_Volcano <- function(DT, color_values) {
        # Volcano plot
        g <- ggplot(dat, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
                geom_point(data=dat[Group=='NS']) +     # Plot insigificant first
                geom_point(data=dat[Group !='NS']) +    # Plot remaining data on layer ABOVE NS
                scale_color_manual(values=color_vals) +     # Force manual colors in `color_values`
                geom_hline(yintercept=2, linetype='dashed') +   # P-value threshold
                geom_vline(xintercept=c(-1,1), linetype='dashed') + # foldchange threshold
                theme_minimal(12)
        return(g)
}

####################################################################################################
# Plotting aesthetics

color_vals <- c('Neurite Heavy' =   '#67a9cfff',    # Blue-ish
                'NS'            =   '#969696ff',    # Grayish
                'Soma Heavy'    =   '#ef8a62ff'     # Orange-ish
)    
####################################################################################################
## Added 20240720

dat <- fread('mean_value_HN_H_HS_L.csv')
dat[, gene := tstrsplit(PG.Genes, split=';')[1]]



# test set is all the non-NA genes in the row_mean.HN_H column,
test_set <- dat[!is.na(row_mean.HN_H), gene]

# background set is all non-NA in row_mean.HS_L column
# (plus the test set)
background <- unique(c(test_set, dat[!is.na(row_mean.HS_L), gene]))

test_set <- symbol_to_entrez(test_set)
background <- symbol_to_entrez(background)

HN_H_vs_HS_L <- getGO(entrez_list = test_set, all_gene_list=background)

plot_enrich(HN_H_vs_HS_L, filestem='HN_H_vs_HS_L', lims=c(-150,210), x.lbl='HN_H_vs_HS_L')




# Format SILAC data

silac_filename <- 'Neutri_H_vs_Soma_H_limma.csv'
dat <- fread(silac_filename)

dat[Group=='DOWN', Group := 'Soma Heavy']
dat[Group=='UP', Group := 'Neurite Heavy']
dat[Group=='Others', Group := 'NS']

dat[, Group := factor(Group, levels=c('Soma Heavy','NS','Neurite Heavy'))]
dat[, logFC := -1*logFC]

dat.silac <- dat
rm(dat)

dat.silac[, gene := tstrsplit(Genes, split=';')[1]]
heavy.neurite.symbol <- dat.silac[Group=='Neurite Heavy', gene]
heavy.soma.symbol <- dat.silac[Group=='Soma Heavy', gene]

heavy.neurite.entrez <- symbol_to_entrez(heavy.neurite.symbol)

heavy.soma.entrez <- symbol_to_entrez(heavy.soma.symbol)

silac_genes.entrez <- symbol_to_entrez(dat.silac$gene)


heavy_neurite_GO <- getGO(entrez_list = heavy.neurite.entrez, all_gene_list=silac_genes.entrez)
heavy_soma_GO <- getGO(entrez_list = heavy.soma.entrez, all_gene_list=silac_genes.entrez)
heavy_soma_GO[, Count := -1*Count]

heavy.combined <- rbindlist(list(heavy_neurite_GO, heavy_soma_GO))

plot_enrich(heavy.combined, filestem='heavy_Neurite_vs_soma', lims=c(-125,100), x.lbl='Heavy peptide enrichment in Neurites')

####################################################################################################
# Format QUANTCAT data

quantcat_filename <- 'QH_L_VS_QH_H_limma.csv'
dat <- fread(quantcat_filename)

dat <- dat[, c('Protein_Group','Genes','QH1.Channel3','QH2.Channel3','QH3.Channel3')]
dat <- melt(dat, measure.vars=grep('QH[123]',colnames(dat), value=T), variable.name='channel')
dat[is.na(value), value := 0]
dat[, rnk := frank(value)]

# Take leading representative from each protein group
dat[, lead_gene := tstrsplit(Genes, split=';')[1]]


lead_gene_order <- dat[, mean(rnk), by=lead_gene][order(V1), lead_gene]

dat[, lead_gene := factor(lead_gene, levels=lead_gene_order)]

in_all_3.symbol <- as.character(dat[value > 0, .N, by=lead_gene][N==3]$lead_gene)        # 136 genes
in_2_or_more.symbol <- as.character(dat[value > 0, .N, by=lead_gene][N>=2]$lead_gene)    # 257 genes


# dat <- merge(dat,  dat[value>0, .N, by=lead_gene], all=TRUE)
# dat[is.na(N), N := 0]


# Get universe of all detected genes
dat.mean_values <- fread('mean_value.csv')
dat.mean_values[, lead_gene := tstrsplit(Genes, split=';')[1]]
all_gene_vector <- unique(as.character(dat.mean_values$lead_gene))
gene_id <- symbol_to_entrez(all_gene_vector)
all_gene_vector <- unique(as.character(na.omit(gene_id)))


# Get test sets


in_all_3.entrez <- symbol_to_entrez(in_all_3.symbol)

in_2_or_more.entrez <- symbol_to_entrez(in_2_or_more.symbol)













# Genes found in all 3 heavy samples
d1 <- as.data.table(enrichGO(gene          = in_all_3.entrez,
                                universe      = all_gene_vector,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = enrich_pvalue,
                                qvalueCutoff  = enrich_pvalue,
                                readable      = TRUE)@result)
d1[, 'group' := 'go_cc']

d2 <- as.data.table(enrichGO(gene          = in_all_3.entrez,
                                universe      = all_gene_vector,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = enrich_pvalue,
                                qvalueCutoff  = enrich_pvalue,
                                readable      = TRUE)@result)
d2[, 'group' := 'go_bp']

d3 <- as.data.table(enrichGO(gene          = in_all_3.entrez,
                                universe      = all_gene_vector,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = enrich_pvalue,
                                qvalueCutoff  = enrich_pvalue,
                                readable      = TRUE)@result)
d3[, 'group' := 'go_mf']

d4 <- as.data.table(enrichKEGG(gene         = in_all_3.entrez,
                                organism     = 'hsa',
                                universe      = all_gene_vector,
                                pAdjustMethod = "BH",
                                qvalueCutoff  = enrich_pvalue,
                                pvalueCutoff = enrich_pvalue)@result)
d4[, 'group' := 'kegg']
desired_cols <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count','group')
d4 <- d4[, .SD, .SDcols=desired_cols]
in_all_3.GO <- rbindlist(list(d1,d2,d3,d4))



# Genes found in 2+ heavy samples
d1 <- as.data.table(enrichGO(gene          = in_2_or_more.entrez,
                                universe      = all_gene_vector,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = enrich_pvalue,
                                qvalueCutoff  = enrich_pvalue,
                                readable      = TRUE)@result)
d1[, 'group' := 'go_cc']

d2 <- as.data.table(enrichGO(gene          = in_2_or_more.entrez,
                                universe      = all_gene_vector,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = enrich_pvalue,
                                qvalueCutoff  = enrich_pvalue,
                                readable      = TRUE)@result)
d2[, 'group' := 'go_bp']

d3 <- as.data.table(enrichGO(gene          = in_2_or_more.entrez,
                                universe      = all_gene_vector,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = enrich_pvalue,
                                qvalueCutoff  = enrich_pvalue,
                                readable      = TRUE)@result)
d3[, 'group' := 'go_mf']

d4 <- as.data.table(enrichKEGG(gene         = in_2_or_more.entrez,
                                organism     = 'hsa',
                                universe      = all_gene_vector,
                                pAdjustMethod = "BH",
                                qvalueCutoff  = enrich_pvalue,
                                pvalueCutoff = enrich_pvalue)@result)
d4[, 'group' := 'kegg']
desired_cols <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count','group')
d4 <- d4[, .SD, .SDcols=desired_cols]
in_2_or_more.GO <- rbindlist(list(d1,d2,d3,d4))



plot_enrich(in_2_or_more.GO, filestem='2_or_more', lims=c(-70,105), x.lbl='Enrichment in heavy peptide MS')
plot_enrich(in_all_3.GO, filestem='all_3', lims=c(-70,65), x.lbl='Enrichment in heavy peptide MS')







quit()





ggsave(plot_MA(dat, color_vals), file='MA-plot.pdf', width=15, height=15, units='cm')

ggsave(plot_Volcano(dat), file='volcano-plot.pdf', width=15, height=15, units='cm')

####################################################################################################
# Prepare QUANCAT data

# things that have incorporated the heavy label in the QH H sample

####################################################################################################
# Plotting functions

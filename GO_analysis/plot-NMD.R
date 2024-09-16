#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(gprofiler2)
library(viridis)
library(ggthemes)
library(org.Hs.eg.db)
library(clusterProfiler)

####################################################################################################



plot_gprofiler_enrich <- function(DT, filestem, lims=NULL, x.lbl=NULL) {
        dat <- copy(DT)
        # Concatenate into one table
        setnames(dat, 'p_value', 'qvalue')
        setnames(dat, 'intersection_size', 'Count')
        setnames(dat, 'term_name', 'Description')
        dat <-dat[qvalue <= 0.05]
        dat[source=='GO:MF', 'Domain' := 'GO Molecular Function']
        dat[source=='GO:CC', 'Domain' := 'GO Cellular Component']
        dat[source=='GO:BP', 'Domain' := 'GO Biological Process']
        dat[source=='KEGG', 'Domain' := 'KEGG Pathway']
        dat[, source := gsub(':', '_', source)]
        dat <- dat[!is.na(Domain)]

        dat.hold <- copy(dat)

        for (grp in c('GO_MF','GO_CC','GO_BP','KEGG')) {
        dat <- copy(dat.hold)
        dat <- dat[source==grp]
        # continue with next if no signif
        if (nrow(dat) == 0) { next }

        # Create qsort column to sort within 'counts' ties
        dat[, 'qsort' := -1*qvalue]

        # Set order for plotting
        setkey(dat, 'Domain', Count, qsort)
        #dat[, Description := stringr::str_wrap(Description, width=35)]
        dat[, Description := paste0(Description, ' ')]
        desc_order <- dat$Description
        dat[, Description := factor(Description, levels=desc_order)]


        g <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
                theme_few() +
                scale_fill_viridis(limits=c(0,0.05), direction=-1) +
                geom_vline(xintercept=0) +
                labs(x=x.lbl) +
                labs(y='GO') +
                #ggplot2::facet_grid(`GO domain` ~ ., drop=T, scales='free', space='free', switch='y') +
                theme(strip.text.y.left = element_text(angle = 0)) +
                geom_text(data=dat, aes(x=0, label=Description, hjust=1)) +
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


sensitive_transcripts <- symbol_to_entrez(readLines('nmd_sensitive_transcripts.txt'))
insensitive_transcripts <- symbol_to_entrez(readLines('nmd_insensitive_transcripts.txt'))

sensitive <- gost(
    sensitive_transcripts,
    organism = "hsapiens",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = NULL,
    numeric_ns = "",
    sources = NULL,
    as_short_link = FALSE,
    highlight = FALSE
)
sensitive.dt <- as.data.table(sensitive$result)
plot_gprofiler_enrich(sensitive.dt, filestem='sensitive_gprofiler_v_annotated', x.lbl='Count', lims=c(-300,300))

insensitive <- gost(
    insensitive_transcripts,
    organism = "hsapiens",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = "fdr",
    domain_scope = "annotated",
    custom_bg = NULL,
    numeric_ns = "",
    sources = NULL,
    as_short_link = FALSE,
    highlight = FALSE
)
insensitive.dt <- as.data.table(insensitive$result)
plot_gprofiler_enrich(insensitive.dt, filestem='insensitive_gprofiler_v_annotated', x.lbl='Count', lims=c(-75,75))


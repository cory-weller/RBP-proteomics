library(data.table)
library(ggplot2)
library(viridis)
library(ggthemes)

upfile <- '/home/wellerca/VR-proteomics/ProtPipe/RNABP_KD/Enrichiment_Analysis/FUS-Neurite-KD_vs_FUS-Neurite-NT/up.enrich_res.tsv'
downfile <- '/home/wellerca/VR-proteomics/ProtPipe/RNABP_KD/Enrichiment_Analysis/FUS-Neurite-KD_vs_FUS-Neurite-NT/down.enrich_res.tsv'


plot_GO <- function(upTable=NULL, downTable=NULL, plotname=NULL, downsample=NULL) {
    upgenes <- NULL
    downgenes <- NULL
    if(!is.null(upTable)) { upgenes <- fread(upTable) }
    if(!is.null(downTable)) {
        downgenes <- fread(downTable)

        # Invert Count for depleted genes
        downgenes[, Count := -1*Count]
    }

    # Concatenate into one table
    dat <- rbindlist(list(upgenes, downgenes))[qvalue <= 0.05]

    # Create qsort column to sort within 'counts' ties
    dat[, 'qsort' := -1*qvalue]

    # Set order for plotting
    setkey(dat, Count, qsort)
    dat[, Description := stringr::str_wrap(Description, width=35)]
    GO_order <- dat$Description
    dat[, Description := factor(Description, levels=GO_order)]

    # Take <downsample> random rows
    dat <- dat[sample(.N, size=downsample, replace=FALSE)]
    max_char <- max(nchar(as.character(dat$Description)))


    g <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
        theme_few() +
        scale_fill_viridis(limits=c(0,0.05), direction=-1) +
        geom_vline(xintercept=0) +
        labs(x='Enrichment in Knockdown relative to non-targeting control') +
        labs(y='GO')

    height_cm <- 5 + 0.5 * nrow(dat)

    ggsave(g, file=paste0(plotname, '.png'), width=25, height=height_cm, units='cm')

}

#
plot_GO(upfile, downfile, '15', 15)


# 10
plot_GO(upfile, downfile, '10', 10)

# 7
plot_GO(upfile, downfile, '7', 7)

# 3
plot_GO(upfile, downfile, '3', 3)

# 2
plot_GO(upfile, downfile, '2', 2)

# 1
plot_GO(upfile, downfile, '1', 1)


n_y <- 
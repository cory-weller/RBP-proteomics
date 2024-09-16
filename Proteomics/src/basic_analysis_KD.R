#!/usr/bin/env Rscript
# R/4
#proteomics analysis for DIA-NN and Spectronaut quantity estimates

#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
    make_option(
        "--pgfile",
        default=NULL,
        help=paste(
            'Input file of Protein Group Intensity (from DIA-NN or Spectronaut)',
            'Required.',
            sep=optparse_indent
        )
    ),
    make_option(
        "--out",
        dest="outdir",
        default='output', 
        help=paste(
            'Directory to direct all output. Directory will be created if does not exist).',
            'Defaults to the current working directory:',
            pwd,
            sep=optparse_indent
        )
    ),
    make_option(
        "--labelgene",
        dest="labelgene",
        default=NULL, 
        help='Gene to always label in output plots'
    ),
    make_option(
      "--heatmap",
      dest="heatmap",
      default=NULL, 
      help='Gene to plot heatmap'
    ),
    make_option(
        "--base",
        dest="log_base",
        default=10, 
        help='Base for log transformation of intensity data. Default: 10'
    ),
    make_option(
        "--normalize",
        default='none',
        type='character',
        help=paste(
            'shift: adjust sample intensities to match global median by adding a constant',
            'scale: adjust sample intensities to match global median by multiplicative scaling',
            'none: do not normalize',
            sep=optparse_indent
        )
    ),
    make_option(
        "--exclude",
        default=NULL,
        type='character',
        help=paste(
            'semicolon-separated string of files to exclude from analysis'
        )
    ),
    make_option(
        "--sds",
        dest = 'sds',
        default=3,
        type='numeric',
        help=paste(
            'Filter out samples with protein group counts > N standard deviations from the mean.',
            'Increase to higher values for greater tolerance of variance in protein group counts.',
            'Default: 3',
            sep=optparse_indent
        )
    ),
    make_option(
        "--minintensity",
        dest = 'minintensity',
        default=0,
        type='numeric',
        help='Minimum LINEAR (not log) intensity. Default: 0'
    ),
    make_option(
        "--fdr",
        dest = 'fdr_threshold',
        default=0.01,
        type='numeric',
        help=paste(
            'False Discovery Rate threshold for differential abundance analysis.',
            'Default: 0.01',
            sep=optparse_indent
        )
    ),
    make_option(
        "--foldchange",
        dest = 'foldchange',
        default=2,
        type='numeric',
        help=paste(
            'Minimum LINEAR fold change [NOT log, as log base can be modified] for labeling',
            'protein groups in differential abundance analysis. Default: 2 (equivalent to',
            'log2 fold-change threshold of 1)',
            sep=optparse_indent
        )
    ),
    make_option(
        "--imputation",
        action = 'store_true',
        default=FALSE, 
        type='logical',
        help=paste(
            'Applies data imputation. Not yet implimented.',
            sep=optparse_indent
        )
    ),
    make_option(
        "--design",
        default=NULL,  
        help=paste(
            'Comma- or tab-delimited, three-column text file specifying the experimental design.',
            'File should contain headers. Header names do not matter; column order DOES matter.',
            'Columns order: <sample_name> <condition> <control>',
            sep=optparse_indent
        )
    ),
    make_option(
        "--neighbors",
        default=15,
        type='numeric',
        help=paste(
            'N Neighbors to use for UMAP. Default: 15',
            sep=optparse_indent
        )
    ),
    make_option(
        "--dry",
        action = 'store_true',
        default=FALSE, 
        type='logical',
        help=paste(
            'Applies data imputation. Not yet implimented.',
            sep=optparse_indent
        )
    ),
    make_option(
      "--enrich",
      dest = 'enrich_pvalue',
      default=0.01,
      type='numeric',
      help=paste(
        'The cutoff of p-value for gene enrichment analysis.',
        'Default: 0.01',
        sep=optparse_indent
      )
    ),
    make_option(
      "--gsea",
      dest = 'gsea_fdr_cutoff',
      default=0.01,
      type='numeric',
      help=paste(
        'The cutoff False Discovery Rate of gsea analysis.',
        'Default: 0.01',
        sep=optparse_indent
      )
    )
)

usage_string <- "Rscript %prog --pgfile [filename] --design [filename] [other options] "
opt <- parse_args(OptionParser(usage = usage_string, option_list))

# opt$

source('src/functions.R')

# Set to TRUE when running interactively for debugging, to set test opts
if(TRUE) {
    opt <- list()
    opt$pgfile <-  'renamed.pg_matrix.tsv'
    opt$design <-  'design_matrix.csv'
    opt$outdir <-  'Proteomics_07022024/'
    opt$sds <-  3
    opt$normalize <-  'shift'
    opt$log_base <-  2
    opt$exclude <-  c('HNRNPA1-Neurite-NT_2','FUS-Soma-KD_5','TDP43-Soma-NT_6','TDP43-Soma-KD_6')
    opt$fdr_threshold <- 0.05
    opt$foldchange <- 2
    opt$minintensity <- 0
    opt$neighbors <- 15
    opt$sds <- 2
    opt$enrich_pvalue <- 0.01
    opt$dry <- FALSE
}

badargs <- FALSE

if(opt$dry) {
    cat("INFO: Quitting due to --dry run\n")
    quit(status=0)
}

if(! opt$normalize %in% c('shift','scale','none')) {
    cat("ERROR: --normalize must be 'shift' 'scale' or 'none'\n")
    badargs <- TRUE
}

if (is.null(opt$pgfile)) {
    cat("ERROR: --pgfile <file> must be provided\n")
    badargs <- TRUE
}

if (badargs == TRUE) {
    quit(status=1)
}

#### PACKAGES ######################################################################################
package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist','ggbeeswarm', 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db','clusterProfiler','pheatmap')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
    cat("INFO: All packages successfully loaded\n")
} else {
    cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on

opt$lfc_threshold <- log(opt$foldchange, base=2)
cat(paste0('INFO: LFC threshold of log2(Intensity) > ', opt$lfc_threshold, '\n'))
cat(paste0('INFO: FDR threshold of ', opt$fdr_threshold, '\n'))

#### MAKE DIRS #####################################################################################

QC_dir <- paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
}

cluster_dir <- paste0(opt$outdir, '/Clustering/')
if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
}

if (!is.null(opt$design)) {
  DI_dir <- paste0(opt$outdir, '/Differential_Intensity/')
  if(! dir.exists(DI_dir)){
    dir.create(DI_dir, recursive = T)
  }
  
  EA_dir <- paste0(opt$outdir, '/Enrichiment_Analysis/')
  if(! dir.exists(EA_dir)){
    dir.create(EA_dir, recursive = T)
  }
}
  
#### IMPORT AND FORMAT DATA#########################################################################

tryTo(paste0('INFO: Reading input file ', opt$pgfile),{
    dat <- fread(opt$pgfile)
}, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))

tryTo(paste0('INFO: Massaging data from ', opt$pgfile, ' into a common style format for processing'), {
    dat <- standardize_format(dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')

tryTo(paste0('INFO: Trimming extraneous column name info'), {
    setnames(dat, trim_colnames(dat))
    #set column order
    col_order=c(colnames(dat)[1:2],sort(colnames(dat)[3:ncol(dat)]))
    setcolorder(dat,col_order)
}, 'ERROR: failed! Check for missing/corrupt headers?')


#Converting to long format
tryTo(paste0('INFO: Converting to long format'), {
  dat.long <- melt_intensity_table(dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')

tryTo('INFO: Excluding all unquantified or zero intensities', {
  dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]
}, 'ERROR: failed!')

tryTo(paste0('INFO: Applying Filter Intensity > ',opt$minintensity),{
  dat.long <- dat.long[Intensity > opt$minintensity]
}, 'ERROR: failed!')

#### QC ############################################################################################

## Plotting intensity distribution

tryTo('INFO: Plotting intensity distribution',{
    plot_pg_intensities(dat.long, QC_dir, 'intensities.pdf', plot_title='Intensity Distribution Pre-normalization')
    plot_pg_intensities(dat.long[! Sample %like% 'FUS'], QC_dir, 'intensities-no-FUS.pdf', plot_title='Intensity Distribution Pre-normalization')
}, 'ERROR: failed!')

## Normalization takes place by default, and can be modified with the --normalize flag. See opts.
if (opt$normalize == 'none') {
    cat('Skipping median-normalization due to --normalize none\n')
} else {
    if (opt$normalize == 'shift') {
        tryTo('INFO: Calculating median-normalized intensities by shifting sample intensities',{
            dat.long <- shift_normalize_intensity(dat.long)
        }, 'ERROR: failed!')

        tryTo('INFO: Plotting shift-normalized intensity distributions',{
          plot_pg_intensities(dat.long, QC_dir, 'intensities_shift_normalized.pdf', plot_title='Intensity Distribution Post-normalization')
          plot_pg_intensities(dat.long[! Sample %like% 'FUS'], QC_dir, 'intensities_shift_normalized-no-FUS.pdf', plot_title='Intensity Distribution Post-normalization')
        }, 'ERROR: failed!')
    } else if (opt$normalize == 'scale') {

        tryTo('INFO: Calculating median-normalized intensities by scaling sample intensities',{
            dat.long <- scale_normalize_intensity(dat.long)
        }, 'ERROR: failed!')

        tryTo('INFO: Plotting scale-normalized intensity distributions',{
            plot_pg_intensities(dat.long, QC_dir, 'intensities_scale_normalized.pdf', plot_title='scale-normalized intensities')
        }, 'ERROR: failed!')
    }

    tryTo('INFO: Re-generating wide table with normalized intensities',{
        original_colorder <- colnames(dat)
        dat <- data.table::dcast(dat.long, Protein_Group+Genes~Sample, value.var='Intensity')
        setcolorder(dat, original_colorder)
    }, 'ERROR: failed!')
}

# pgcounts represents the distribution of Protein Groups with Intensity > 0
# Visually, it is represented as a bar plot with x=sample, y=N, ordered by descending N
# Get counts of [N=unique gene groups with `Intensity` > 0]
tryTo('INFO: Tabulating protein group counts',{
    pgcounts <- dat.long[, .N, by=Sample]
    # Order samples by ascending counts
    ezwrite(pgcounts, QC_dir, 'protein_group_nonzero_counts.tsv')
    plot_pg_counts(pgcounts, QC_dir, 'protein_group_nonzero_counts.pdf')
    plot_pg_counts(pgcounts[! Sample %like% 'FUS'], QC_dir, 'protein_group_nonzero_counts-no-FUS.pdf')
}, 'ERROR: failed!')

tryTo('INFO: Plotting sample intensity correlations',{
    dat.correlations <- get_spearman(dat)
    ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
    plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.pdf')
    dat.correlations.NOFUS <- dat.correlations[! SampleA %like% 'FUS'][! SampleB %like% 'FUS']
    plot_correlation_heatmap(dat.correlations.NOFUS, QC_dir, 'sample_correlation-no-FUS.pdf')
}, 'ERROR: failed!')


if (!is.null(opt$design)) {
  tryTo('INFO: Importing experimental design',{
    design <- fread(opt$design, header=TRUE)
    setnames(design, c('sample_name', 'condition', 'control'))
  }, 'ERROR: failed!')
  tryTo('INFO: Validating experimental design',{
    print(design[])
    cat('\n')
    conditions <- unique(design$condition)
    for (condition.i in conditions) {
      samples <- design[condition == condition.i, sample_name]
      control <- unique(design[condition == condition.i, control])
      if (length(control) != 1) {
        cat(paste0('ERROR: condition ', condition.i, ' maps to multiple controls: ', control, '\n'))
        cat(paste0('       Check the design matrix and esure no more than one control label per condition\n'))
        quit(exit=1)
      } else {
        cat(paste0('INFO: condition ', condition.i, ' maps to control ', control, '\n'))
      }
    }
    cat(paste0('INFO: all conditions pass check (i.e. map to one control condition)\n'))
  }, 'ERROR: failed!')
}




## Exclude samples with N protein groups < opt$sds away from mean
## Default value: 3 standard deviations, modifiable with --sds [N]
tryTo('INFO: Identifying samples with protein group count outliers',{
    cat(paste0('INFO: defining outliers as samples with [N protein groups] > ', opt$sds, ' standard deviations from the mean\n'))
    stdev <- sd(pgcounts[,N])
    mean_count <- mean(pgcounts[,N])
    min_protein_groups <- floor(mean_count - (opt$sds * stdev))
    max_protein_groups <- ceiling(mean_count + (opt$sds * stdev))
    cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
    low_count_samples <- as.character(pgcounts[N < min_protein_groups, Sample])
    high_count_samples <- as.character(pgcounts[N > max_protein_groups, Sample])
    if(length(low_count_samples)==0) {
        cat('\nINFO: No low group count samples to remove\n')
    } else {
        cat(paste0('\nINFO: Pruning low-count outlier ', low_count_samples))
        cat('\n\n')
        print(pgcounts[Sample %in% low_count_samples])
        cat('\n')
        dat[, c(low_count_samples) := NULL]    # remove sample columns from wide table
        dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples)==0) {
        cat('INFO: No high group count samples to remove\n')
    } else {
        cat(paste0('\nINFO: Pruning high-count outlier ', high_count_samples))
        cat('\n')
        print(pgcounts[Sample %in% high_count_samples])
        dat[, c(high_count_samples) := NULL]    # remove sample columns from wide table
        dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    }
}, 'ERROR: failed!')

if (!is.null(opt$exclude)) {
  tryTo(paste('INFO: excluding samples', opt$exclude),{
    design <- design[! sample_name %in% opt$exclude]
    dat[, opt$exclude := NULL]
    dat.long <- dat.long[! Sample %in% opt$exclude]
  }, 'ERROR: failed!')
}

# #### CLUSTERING ####################################################################################
# # PCA
# tryTo('INFO: running PCA and plotting first two components',{
#     pca <- get_PCs(dat)
#     ezwrite(pca$components, cluster_dir, 'PCA-all.tsv')
#     ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
# }, 'ERROR: failed!')

get_pca_subset <- function(DT, col_pattern) {
    DT <- DT[, .SD, .SDcols=c('Protein_Group','Genes', grep(col_pattern, colnames(DT), value=T))]
    out <- get_PCs(DT)
    return(out)
}

plot_RNABP_PCs <- function(PCAobj, prefix) {
    DT <- as.data.table(PCAobj$components)
    PCAsummary <- as.data.table(PCAobj$summary)
    pc1 <- PCAobj$summary[component=='PC1', percent]
    pc2 <- PCAobj$summary[component=='PC2', percent]
    DT[, c('gene','location','treatment','replicate') := tstrsplit(sample, split='-|_')]
    DT[treatment=='KD' & gene == 'FUS', fillcolor := '#FF89D7']
    DT[treatment=='KD' & gene == 'HNRNPA1', fillcolor := '#0095FF']
    DT[treatment=='KD' & gene == 'TDP43', fillcolor := '#D682FF']
    DT[treatment=='NT', fillcolor := NA]

    DT[, fillcolor := factor(fillcolor)]

    g <- ggplot(DT, aes(x=PC1, y=PC2, color=gene, shape=location, fill=fillcolor)) + 
    geom_point(size=2) +
    scale_color_manual(values=c('FUS'='#FF89D7','HNRNPA1'='#0095FF','TDP43'='#D682FF')) +
    scale_shape_manual(values=c('Neurite'=21, 'Soma'=22)) +
    scale_fill_manual(values=levels(DT$fillcolor), na.value='transparent') +
    theme_minimal() +
    guides(fill = 'none',
          shape = guide_legend(override.aes = list(fill = "black"))) +
    labs(title='Proteomics PCA -- empty shapes are non-targeting controls') +
    labs(x=paste0('PC1 (', pc1, '%)')) +
    labs(y=paste0('PC2 (', pc2, '%)'))

    ggsave(g, file=paste0(prefix,'.pca.pdf'), width=15, height=15, units='cm')

    DT[, fillcolor := NULL]

    fwrite(DT, paste0(prefix, '.pca.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
    fwrite(PCAsummary, paste0(prefix, '.var-explained.tsv'), quote=F, row.names=F, col.names=T, sep='\t')
    
}

# Plot ALL proteomics samples
plot_RNABP_PCs(get_pca_subset(dat, '[0-9]$'), paste0(cluster_dir,'/Proteomics-ALL'))

# Plot FUS
plot_RNABP_PCs(get_pca_subset(dat, 'FUS-'), paste0(cluster_dir,'/FUS-Both'))
plot_RNABP_PCs(get_pca_subset(dat, 'FUS-Soma-'), paste0(cluster_dir,'/FUS-Soma'))
plot_RNABP_PCs(get_pca_subset(dat, 'FUS-Neurite-'), paste0(cluster_dir,'FUS-Neurite'))

# Plot TDP43
plot_RNABP_PCs(get_pca_subset(dat, 'TDP43-'), paste0(cluster_dir,'/TDP43-Both'))
plot_RNABP_PCs(get_pca_subset(dat, 'TDP43-Soma-'), paste0(cluster_dir,'TDP43-Soma'))
plot_RNABP_PCs(get_pca_subset(dat, 'TDP43-Neurite-'), paste0(cluster_dir,'TDP43-Neurite'))

# Plot HNRNPA1
plot_RNABP_PCs(get_pca_subset(dat, 'HNRNPA1-'), paste0(cluster_dir,'/HNRNPA1-Both'))
plot_RNABP_PCs(get_pca_subset(dat, 'HNRNPA1-Soma-'), paste0(cluster_dir,'HNRNPA1-Soma'))
plot_RNABP_PCs(get_pca_subset(dat, 'HNRNPA1-Neurite-'), paste0(cluster_dir,'HNRNPA1-Neurite'))

genes_to_label <- c('TARDBP','FUS','HNRNPA1','HNRNPA3','TAF15','STMN2')
opt$labelgene <- paste0('^', genes_to_label, '$')
label_gene_pattern <- paste0(opt$labelgene, collapse='|')
### DIFFERENTIAL INTENSITY########################################################################
if (!is.null(opt$design)) {
  tryTo('INFO: Running differential intensity t-tests and pathway analysis',{
    for (treatment in conditions) {
      print(paste0(treatment, ' vs ', control, ' DE analysis'))
      treatment_sample_names <- intersect(colnames(dat), design[condition == treatment, sample_name])
      if (length(treatment_sample_names)==0) {next}
      control <- unique(design[condition == treatment, control])
      control_sample_names <- colnames(dat)[colnames(dat) %like% control]
      if (length(control_sample_names)==0) {next}
      if(treatment != control) {
        t_test <- do_t_test(dat, treatment_sample_names, control_sample_names)
        ezwrite(t_test[order(p.adj)], DI_dir, paste0(treatment, '_vs_', control, '.tsv'))
        plot_volcano(t_test, treatment, control, opt$log_base, opt$lfc_threshold, opt$fdr_threshold, DI_dir, label_gene_pattern)
        enrich_pathway(t_test, treatment, control, EA_dir, opt$lfc_threshold, opt$fdr_threshold,opt$enrich_pvalue)
      }
    }
  }, 'ERROR: failed!')
}

library(viridis)
library(ggthemes)


plot_enrich <- function(dir_path=NULL, lims=NULL) {
    upfile <- list.files(dir_path, 'up.*.tsv', full.names=T)
    downfile <- list.files(dir_path, 'down.*.tsv', full.names=T)
    upgenes <- NULL
    downgenes <- NULL
    if(length(upfile) > 0) { upgenes <- fread(upfile) }
    if(length(downfile) > 0) {
        downgenes <- fread(downfile) 
        # Invert Count for depleted genes
        downgenes[, Count := -1*Count]
    }

    # Concatenate into one table
    dat <- rbindlist(list(upgenes, downgenes))[qvalue <= 0.05]
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

        # Take <downsample> random rows
        # dat <- dat[sample(.N, size=downsample, replace=FALSE)]
        # max_char <- max(nchar(as.character(dat$Description)))


        g <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
            theme_few() +
            scale_fill_viridis(limits=c(0,0.05), direction=-1) +
            geom_vline(xintercept=0) +
            labs(x='Enrichment in Knockdown relative to non-targeting control') +
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
        plotname <- basename(dir_path)
        ggsave(g, file=paste0(dir_path, '/', plotname, '.', grp, '.pdf'), width=40, height=height_cm, units='cm', limitsize = FALSE)

    }

}

plot_enrich(dir_path='Proteomics_07022024/Enrichiment_Analysis/FUS-Neurite-KD_vs_FUS-Neurite-NT', lims=c(-3.5,2.5))
plot_enrich(dir_path='Proteomics_07022024/Enrichiment_Analysis/FUS-Soma-KD_vs_FUS-Soma-NT', lims=c(-2.5,4))
plot_enrich(dir_path='Proteomics_07022024/Enrichiment_Analysis/HNRNPA1-Neurite-KD_vs_HNRNPA1-Neurite-NT', lims=c(-150,110))
plot_enrich(dir_path='Proteomics_07022024/Enrichiment_Analysis/HNRNPA1-Soma-KD_vs_HNRNPA1-Soma-NT', lims=c(-1,3))
plot_enrich(dir_path='Proteomics_07022024/Enrichiment_Analysis/TDP43-Neurite-KD_vs_TDP43-Neurite-NT', lims=c(-12,4))
plot_enrich(dir_path='Proteomics_07022024/Enrichiment_Analysis/TDP43-Soma-KD_vs_TDP43-Soma-NT', lims=c(-8,3))








quit()
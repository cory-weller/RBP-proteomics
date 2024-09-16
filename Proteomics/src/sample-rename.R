#!/usr/bin/env Rscript

library(data.table)

dat <- data.table(filename =list.files('mzML', pattern='*.mzML$'))
dat[, 'sample' := tstrsplit(filename, '.mzml', split='\\.mzML')[1]]
dat[, 'replicate' := substr(sample, nchar(sample), nchar(sample))]
dat[, 'replicate' := as.numeric(replicate)]

dat[sample %like% '011023', gene := 'FUS']
dat[sample %like% '011323', gene := 'HNRNPA1']
dat[sample %like% '100822', gene := 'TDP43']
dat[sample %like% 'N_{0,1}[1-6]', location := 'Neurite']
dat[sample %like% 'S_{0,1}[1-6]', location := 'Soma']

# Fix one differently-labeled sample
dat[sample %like% '011323_VHR_1S_5', gene := 'FUS']


dat[gene=='FUS' & filename %like% '_1[NS]_[1-6]\\.', treatment := 'NT']
dat[gene=='FUS' & filename %like% '_52[NS]_[1-6]\\.', treatment := 'KD']

dat[gene=='HNRNPA1' & filename %like% '_6[NS]_[1-6]\\.', treatment := 'NT']
dat[gene=='HNRNPA1' & filename %like% '_8[NS]_[1-6]\\.', treatment := 'KD']

dat[gene=='TDP43' & filename %like% '_1[NS][1-6]\\.', treatment := 'NT']
dat[gene=='TDP43' & filename %like% '_2[NS][1-6]\\.', treatment := 'KD']

dat <- dat[, .SD, .SDcols=c('sample','gene','location','treatment','replicate')]
setkey(dat, gene, location, treatment, replicate)


# Export conditions table
fwrite(dat, file='conditions.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# Build design matrix table for ProtPipe
design_matrix <- dat[treatment == 'KD']
design_matrix[, sample_name := paste0(gene, '-', location, '-', treatment, '_', replicate)]
design_matrix[, condition := paste0(gene, '-', location, '-', treatment)]
design_matrix[, control := paste0(gene, '-', location, '-', 'NT')]
fwrite(design_matrix[, .SD, .SDcols=c('sample_name','condition','control')], file='design_matrix.csv')





# Rename columns in DIA-NN report to be parsable by ProtPipe
dat[, oldname := paste0('mzML/', sample, '.mzML')]
dat[, newname := paste0('mzML/', gene, '-', location, '-', treatment, '_', replicate)]
abundance <- fread('DIA/report.pg_matrix.tsv')
setnames(abundance, dat$oldname, dat$newname)
fwrite(abundance, 'renamed.pg_matrix.tsv', quote=F, row.names=F, col.names=T, sep='\t')

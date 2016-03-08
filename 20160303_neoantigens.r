
library(ggplot2)
setwd('/Volumes/yjiao/prepost/reneoantigendetection/out/')
options(jupyter.plot_mimetypes = 'image/png')

files <- list.files(pattern='.tsv')
patients <- sapply(files, function(x) {
    strsplit(x, '\\.')[[1]][1]
})
names(files) <- patients

cutoff <- 500
#sapply(patients, function(pID){
pID <- '208T'
neoag <- read.table(files[pID], header=TRUE, sep='\t', stringsAsFactors = FALSE)

options(repr.plot.width=4.5, repr.plot.height=3)
neoag$is_neoantigen <- neoag$min_aff_mt <= cutoff # McGranahan definition
neoag$is_neoantigen2 <- neoag$min_aff_mt <= cutoff & neoag$min_aff_wt > cutoff & neoag$min_aff_wt >= 3*neoag$min_aff_mt # old definition

c <- ggplot(neoag, aes(x=log(min_aff_wt), y=log(min_aff_mt)))
c + geom_point(alpha=.1, aes(color=is_neoantigen), size=1) +
    labs(title='Mininum Affinity: McGranahan', y='Log Mutant', x='Log Wildtype') +
    theme_classic() +
    scale_color_manual(values=c('#000000','#00BFFF'),
                     name='Classification',
                     breaks=c(TRUE, FALSE),
                     labels=c('Neoantigen', 'Non-neoantigen'))

c <- ggplot(neoag, aes(x=log(min_aff_wt), y=log(min_aff_mt)))
c + geom_point(alpha=.1, aes(color=is_neoantigen2), size=1) +
    labs(title='Mininum Affinity: Original', y='Log Mutant', x='Log Wildtype') +
    theme_classic() +
    scale_color_manual(values=c('#000000','#00BFFF'),
                     name='Classification',
                     breaks=c(TRUE, FALSE),
                     labels=c('Neoantigen', 'Non-neoantigen'))
#})

groups <- tapply(1:nrow(neoag), neoag$Hugo_Symbol, function(idx){
    temp <- neoag[idx[1],]
    temp$min_aff_mt <- min(neoag[idx,]$min_aff_mt)
    temp$min_aff_wt <- min(neoag[idx,]$min_aff_mt)
    temp
})
groups <- do.call(rbind, groups)

options(repr.plot.width=3, repr.plot.height=3)
c <- ggplot(groups, aes(x=min_aff_wt, y=log(expression_allele_tpm)))
c + geom_point(size=1, alpha=.1) + scale_alpha(range = c(0,1))
c <- ggplot(groups, aes(x=min_aff_wt <= 500, y=log(expression_allele_tpm)))
c + geom_boxplot()

c <- ggplot(groups, aes(x=min_aff_mt, y=log(expression_allele_tpm)))
c + geom_point(size=1, alpha=.1) + scale_alpha(range = c(0,1))
c <- ggplot(groups, aes(x=is_neoantigen, y=log(expression_allele_tpm)))
c + geom_boxplot()

options(repr.plot.width=8, repr.plot.height=3)
c <- ggplot(neoag[neoag$ccf_hat > 0,], aes(x=is_neoantigen, y=ccf_hat)) #filter out all the ccf_hat = 0's
c + geom_boxplot() + facet_wrap(~order, nrow=1)


options(repr.plot.width=8, repr.plot.height=3)
neorate <- lapply(unique(neoag$order), function(i){
    ineoag <- neoag[neoag$order == i,]
    neorate <- lapply(seq(0.1,1,.1), function(thresh){
        subclonal <- ineoag$ccf_hat < thresh
        dat1 <- sum(ineoag[subclonal, 'is_neoantigen'])/sum(subclonal)
        dat2 <- sum(ineoag[!subclonal, 'is_neoantigen'])/sum(!subclonal)
        out <- data.frame(CCF_threshold=c(thresh, thresh))
        out$rate <- c(dat1, dat2)
        out$order <- c(i, i)
        out$class <- c('subclonal', 'clonal')
        out
    })
    neorate <- do.call(rbind, neorate)
})
neorate <- do.call(rbind, neorate)

c <- ggplot(neorate, aes(x=CCF_threshold, color=class, y=rate))
c + 
    geom_line() + 
    geom_point() +
    expand_limits(y=0) +
    facet_wrap(~order)

#neoag$is_neoantigen <- neoag$min_aff_mt <= cutoff & neoag$min_aff_wt > cutoff & neoag$min_aff_wt >= 3*neoag$min_aff_mt

options(repr.plot.width=8, repr.plot.height=3)
neorate <- lapply(unique(neoag$order), function(i){
    ineoag <- neoag[neoag$order == i,]
    neorate <- lapply(seq(0.1,1,.1), function(thresh){
        subclonal <- ineoag$ccf_hat < thresh
        dat1 <- sum(ineoag[subclonal, 'is_neoantigen2'])/sum(subclonal)
        dat2 <- sum(ineoag[!subclonal, 'is_neoantigen2'])/sum(!subclonal)
        out <- data.frame(CCF_threshold=c(thresh, thresh))
        out$rate <- c(dat1, dat2)
        out$order <- c(i, i)
        out$class <- c('subclonal', 'clonal')
        out
    })
    neorate <- do.call(rbind, neorate)
})
neorate <- do.call(rbind, neorate)

c <- ggplot(neorate, aes(x=CCF_threshold, color=class, y=rate))
c + 
    geom_line() + 
    geom_point() +
    expand_limits(y=0) +
    facet_wrap(~order)


options(repr.plot.width=2, repr.plot.height=2)
nsamples <- length(unique(neoag$order))
ags <- neoag[neoag$is_neoantigen,]
grouped <- tapply(1:nrow(ags), ags$Hugo_Symbol, function(idx){ ags[idx,]})
nNeoPerGene <- data.frame(nNeoAg = sapply(grouped, nrow)/nsamples)

ggplot(nNeoPerGene, aes(x=nNeoAg)) + geom_histogram(bins=10)



#!/usr/bin/env Rscript
# csaw.r bam_files num_reps_cond1 output_prefix narrow/broad/TF trended/composition

# This script is a wrapper of csaw that follows the workflows demonstrated in the
# chipseqDB package (https://bioconductor.org/packages/3.11/workflows/html/chipseqDB.html).

# <bam_files> should be separated by a comma.
# Each BAM file must be sorted (by the genomic coordinates of reads) and indexed.

# The first <num_reps_cond1> samples are taken as the replicates for condition 1.
# The remaining samples are considered to belong to condition 2.
# At least 1 condition should be associated with at least two replicates.
# Resulting Mval = log2(condition_2 / condition_1)

# <narrow/broad/TF> specifies the type of ChIP-seq data.
# The three options are for histone modifications associated with relatively narrow peaks,
# histone modifications with broad peaks, and transcription factors respectively.

# <trended/composition> specifies the normalization strategy to be applied.
# The two options are for normalizing away trended and composition biases respectively.

library(rtracklayer)
library(csaw)
library(edgeR)

argv <- commandArgs(TRUE)
if (length(argv) < 5) {
    cat("Usage:", sep = "\n")
    cat("csaw.r bam_files num_reps_cond1 output_prefix narrow/broad/TF trended/composition",
        sep = "\n")
    cat("", sep = "\n")
    stop("Missing mandatory arguments!")
}


## Parameter settings.

# Modify blacklist and frag.len as needed.
blacklist <- "/home/tushiqi/tushiqi/genome/hg19/hg19-blacklist.v2.bed"
minq <- 20
frag.len <- 200
trended_binWidth <- 2000
composition_binWidth <- 10000
merge.tol <- 100

if (argv[4] == "narrow") {
    win.width <- 150
    win.spacing <- 50
    filter.min.fc <- 3
    merge.max.width <- 5000
} else if (argv[4] == "broad") {
    win.width <- 2000
    win.spacing <- 500
    filter.min.fc <- 2
    merge.max.width <- 30000
    consolidation.win.width <- 500
    consolidation.win.spacing <- 100
} else if (argv[4] == "TF") {
    win.width <- 10
    win.spacing <- 50
    filter.min.fc <- 3
    merge.max.width <- 5000
} else {
    stop("The 4th argument is invalid!")
}


## Counting reads into windows.
blacklist <- import(blacklist)
param <- readParam(minq = minq, discard = blacklist)
bams <- strsplit(argv[1], split = ",", fixed = TRUE)[[1]]
win.data <- windowCounts(bams, param = param, width = win.width, ext = frag.len,
                         spacing = win.spacing)


## Filtering and normalization.
if (argv[5] == "trended") {
    # Filtering windows by abundance.
    bins <- windowCounts(bams, bin = TRUE, width = trended_binWidth, param = param)
    filter.stat <- filterWindowsGlobal(win.data, bins)
    keep <- filter.stat$filter > log2(filter.min.fc)
    filtered.data <- win.data[keep, ]

    # Normalizing for sample-specific trended biases.
    filtered.data <- normOffsets(filtered.data)

} else if (argv[5] == "composition") {
    # Normalization for composition biases.
    bins <- windowCounts(bams, bin = TRUE, width = composition_binWidth, param = param)
    win.data <- normFactors(bins, se.out = win.data)

    # Filtering out low-abundance windows.
    filter.stat <- filterWindowsGlobal(win.data, bins)
    keep <- filter.stat$filter > log2(filter.min.fc)
    filtered.data <- win.data[keep, ]

} else {
    stop("The 5th argument is invalid!")
}


## Statistical modelling of biological variability.

# Construct the design matrix.
x <- as.integer(argv[2])
celltype <- factor(c(rep("c1", x), rep("c2", length(bams) - x)))
design <- model.matrix(~ 0 + celltype)
colnames(design) <- levels(celltype)

# Model fitting.
y <- asDGEList(filtered.data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)


## Testing for differential binding and controlling FDR.

contrast <- makeContrasts(c2-c1, levels = design)
res <- glmQLFTest(fit, contrast = contrast)

if (argv[4] == "broad") {
    win.data2 <- windowCounts(bams, param = param, width = consolidation.win.width,
                              ext = frag.len, spacing = consolidation.win.spacing)
    if (argv[5] == "trended") {
        # Filtering windows by abundance.
        filter.stat2 <- filterWindowsGlobal(win.data2, bins)
        keep2 <- filter.stat2$filter > log2(filter.min.fc)
        filtered.data2 <- win.data2[keep2, ]

        # Normalizing for sample-specific trended biases.
        filtered.data2 <- normOffsets(filtered.data2)

    } else if (argv[5] == "composition") {
        # Re-using the same normalization factors.
        win.data2$norm.factors <- win.data$norm.factors

        # Filtering out low-abundance windows.
        filter.stat2 <- filterWindowsGlobal(win.data2, bins)
        keep2 <- filter.stat2$filter > log2(filter.min.fc)
        filtered.data2 <- win.data2[keep2, ]

    }
    y2 <- asDGEList(filtered.data2)
    y2 <- estimateDisp(y2, design)
    fit2 <- glmQLFit(y2, design, robust = TRUE)
    res2 <- glmQLFTest(fit2, contrast = contrast)
    merged <- mergeResultsList(list(filtered.data, filtered.data2),
                               tab.list = list(res$table, res2$table),
                               equiweight = TRUE, tol = merge.tol,
                               merge.args = list(max.width = merge.max.width))
} else {
    merged <- mergeResults(filtered.data, res$table, tol = merge.tol,
                           merge.args = list(max.width = merge.max.width))
}
tabcom <- merged$combined
tabbest <- merged$best


## Writing results.
export(merged$regions, paste0(argv[3], ".csaw_", argv[5], ".bed"))
out <- data.frame(Mval = tabbest$rep.logFC, pval = tabcom$PValue, fdr = tabcom$FDR)
write.table(out, file = paste0(argv[3], ".csaw_", argv[5], ".res"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




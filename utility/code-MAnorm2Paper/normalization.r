## Using the H3K4me3 ChIP-seq data for GM12891 and GM12892 LCLs (lymphoblastoid cell lines)
## to illustrate how to use size factor, TMM and quantile normalization methods to derive
## normalized log2 read counts.

# Read data.
profile_bins <- read.table("data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
head(profile_bins)
n <- 6
count <- profile_bins[4:(3 + n)]
rownames(count) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                         sep = "_")


## For size factor.

# Derive size factors based on the median-ratio strategy.
library(DESeq2)
condition <- rep("x", n)
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = data.frame(condition = condition),
                              design = ~1)
dds <- estimateSizeFactors(dds)
sizeFactor <- dds@colData@listData$sizeFactor

# Perform normalization.
res1 <- t(apply(count, 1, function(y) log2(y / sizeFactor + 0.5)))
head(res1)


## For TMM.

# Derive effective library sizes.
library(limma)
library(edgeR)
dge <- DGEList(counts = count)
dge <- calcNormFactors(dge)

# Perform normalization.
design <- matrix(1, nrow = ncol(count), ncol = 1)
rownames(design) <- names(count)
colnames(design) <- "mu"
res2 <- voom(dge, design)$E
head(res2)


## For quantile normalization.
res3 <- normalizeBetweenArrays(log2(as.matrix(count) + 0.5), method = "quantile")
head(res3)




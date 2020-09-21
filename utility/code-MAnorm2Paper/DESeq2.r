## Using the differential analysis of H3K4me3 ChIP-seq data between GM12891 and GM12892
## LCLs (lymphoblastoid cell lines) to illustrate the workflow of DESeq2.
library(DESeq2)


## Apply DESeq2 to raw read counts.

# Read data.
profile_bins <- read.table("data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
head(profile_bins)
n1 <- 3
n <- 6
count <- profile_bins[4:(3 + n)]
rownames(count) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                         sep = "_")
condition <- rep("cond2", n)
condition[1:n1] <- "cond1"

# Instantiate a DESeqDataSet.
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = data.frame(condition = condition),
                              design = ~condition)

# Perform normalization and estimation of dispersion.
dds <- DESeq(dds)

# Identify differential signals.
res <- results(dds, contrast = c("condition", "cond2", "cond1"))
head(res)
summary(res)


## Apply DESeq2 to normalized log2 read counts.

# Read data.
normalized <- read.table("data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
                         header = TRUE, stringsAsFactors = FALSE)
head(normalized)
rownames(normalized) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                              sep = "_")

# Instantiate a DESeqDataSet (still using raw read counts).
dds2 <- DESeqDataSetFromMatrix(countData = count,
                               colData = data.frame(condition = condition),
                               design = ~condition)

# Design normalization factors.
nf <- as.matrix(count) / (2 ^ as.matrix(normalized))
nf <- apply(nf, 2, function(x) {
    y <- x[x > 0]
    a <- exp(mean(log(y)))
    x[x <= 0] <- a
    x
})
normalizationFactors(dds2) <- nf

# Follow the remaining original workflow for differential analysis.
dds2 <- DESeq(dds2)
res2 <- results(dds2, contrast = c("condition", "cond2", "cond1"))
head(res2)
summary(res2)




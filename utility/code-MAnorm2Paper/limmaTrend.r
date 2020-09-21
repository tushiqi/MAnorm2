## Using the differential analysis of H3K4me3 ChIP-seq data between GM12891 and GM12892
## LCLs (lymphoblastoid cell lines) to illustrate the workflow of limma-trend.
library(limma)
library(edgeR)


## Apply limma-trend to raw read counts.

# Read data.
profile_bins <- read.table("data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
head(profile_bins)
n1 <- 3
n <- 6
count <- profile_bins[4:(3 + n)]
rownames(count) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                         sep = "_")

# Construct the design matrix.
design <- cbind(cond1 = 1, cond2_minus_cond1 = c(rep(0, n1), rep(1, n - n1)))
rownames(design) <- names(count)

# Filter out genomic intervals with low read counts.
min_count <- 6
min_sample <- 1
flag <- apply(count, 1, function(x){ sum(x >= min_count) >= min_sample })
count <- count[flag, ]

# Instantiate a DGEList and perform normalization.
dge <- DGEList(counts = count)
dge <- calcNormFactors(dge)

# Identify differential signals.
logCPM <- cpm(dge, log = TRUE, prior.count = 3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend = TRUE)
topTable(fit, coef = 2, number = 10, adjust.method = "BH", sort.by = "P")


## Apply limma-trend to normalized log2 read counts.

# Read data.
normalized <- read.table("data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
                         header = TRUE, stringsAsFactors = FALSE)
head(normalized)
rownames(normalized) <- paste(profile_bins$chrom, profile_bins$start, profile_bins$end,
                              sep = "_")

# Follow the original workflow for differential analysis.
design2 <- design
rownames(design2) <- names(normalized)
fit2 <- lmFit(as.matrix(normalized), design2)
fit2 <- eBayes(fit2, trend = TRUE)
topTable(fit2, coef = 2, number = 10, adjust.method = "BH", sort.by = "P")




## Using the differential analysis of H3K4me3 ChIP-seq data between GM12891 and GM12892
## LCLs (lymphoblastoid cell lines) to illustrate the workflow of MAnorm2.
library(MAnorm2)


## Apply MAnorm2 to raw read counts.

# Read data.
profile_bins <- read.table("data/GM12891_cmp_GM12892.H3K4me3_profile_bins.xls",
                           header = TRUE, stringsAsFactors = FALSE)
head(profile_bins)
n1 <- 3
n <- 6

# Perform within-group normalization.
norm <- normalize(profile_bins, count = 4:(3 + n1), occupancy = 4:(3 + n1) + n)
norm <- normalize(norm, count = (4 + n1):(3 + n), occupancy = (4 + n1):(3 + n) + n)

# Construct a bioCond for each group of samples.
conds <- list(GM12891 = bioCond(norm[4:(3 + n1)], norm[4:(3 + n1) + n], name = "GM12891"),
              GM12892 = bioCond(norm[(4 + n1):(3 + n)], norm[(4 + n1):(3 + n) + n], name = "GM12892"))

# Perform between-group normalization.
# Restrict common peak regions to autosomes since the two LCLs are of different genders.
autosome <- !(profile_bins$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)

# Fit a mean-variance curve.
conds <- fitMeanVarCurve(conds, method = "parametric", init.coef = c(0.1, 10))

# Perform statistical tests for identifying differential signals.
res <- diffTest(conds[[1]], conds[[2]])
head(res)


## Apply MAnorm2 to normalized log2 read counts.

# Read data.
normalized <- read.table("data/GM12891_cmp_GM12892.H3K4me3_hieMA.xls",
                         header = TRUE, stringsAsFactors = FALSE)
head(normalized)
norm2 <- profile_bins
norm2[4:(3 + n)] <- normalized

# Follow the original workflow for differential analysis.
conds2 <- list(GM12891 = bioCond(norm2[4:(3 + n1)], norm2[4:(3 + n1) + n], name = "GM12891"),
               GM12892 = bioCond(norm2[(4 + n1):(3 + n)], norm2[(4 + n1):(3 + n) + n], name = "GM12892"))
conds2 <- fitMeanVarCurve(conds2, method = "parametric", init.coef = c(0.1, 10))
res2 <- diffTest(conds2[[1]], conds2[[2]])
head(res2)




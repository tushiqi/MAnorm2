# Testing R code for comparing H3K4me3 ChIP-seq signals between GM19238 and
# GM19239 cell lines.
#
# Last update: 2019-01-02


## Make a comparison of H3K4me3 ChIP-seq signals between GM19238 and GM19239
## cell lines.
library(MAnorm2)

# Prepare the dataset.
cnt <- read.table("../data/LCLs_Yoruban.H3K4me3/LCLs_Yoruban.H3K4me3_profile_bins.xls",
                  header=TRUE, stringsAsFactors=FALSE)
cnt <- cnt[c(1:3, 10:15, 22:27)]

# Perform MA normalization and construct bioConds.
norm <- normalize(cnt, 4:6, 10:12)
norm <- normalize(norm, 7:9, 13:15)
conds <- list(GM19238 = bioCond(norm[4:6], norm[10:12], name = "GM19238"),
              GM19239 = bioCond(norm[7:9], norm[13:15], name = "GM19239"))
autosome <- !(cnt$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
png("GM19238_cmp_GM19239.H3K4me3_NormMAplot.png")
MAplot(conds[[1]], conds[[2]], main = "GM19238 vs. GM19239")
abline(h = 0, lwd = 2, lty = 5)
dev.off()

# Fit a mean-variance curve.
conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE,
                         init.coef = c(0.1, 10))
png("GM19238_cmp_GM19239.H3K4me3_MVC.png")
plotMeanVarCurve(conds, subset = "occupied")
dev.off()

# Perform differential tests.
res <- diffTest(conds[[1]], conds[[2]])
png("GM19238_cmp_GM19239.H3K4me3_DiffMAplot.png")
MAplot(res, padj = 0.001, main = "GM19238 vs. GM19239")
abline(h = 0, lwd = 2, lty = 5, col = "green3")
dev.off()

# Write the test results.
write.table(res, file = "GM19238_cmp_GM19239.H3K4me3_diffTest.xls", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)



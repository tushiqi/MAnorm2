# Example R code for demonstrating the usage of MAnorm2 functions.
#
# Maximum column width: 77
#
# Last update: 2018-12-31


## Load the dataset.
data(H3K27Ac, package = "MAnorm2")
attr(H3K27Ac, "metaInfo")


## Estimate size factors
# Use all the genomic intervals.
estimateSizeFactors(H3K27Ac[4:8])
# Use only the genomic intervals occupied by all the ChIP-seq samples.
estimateSizeFactors(H3K27Ac[4:8], subset = apply(H3K27Ac[9:13], 1, all))


## Perform MA normalization on the whole set of ChIP-seq samples once for
## all.
# Exclude the genomic intervals in sex chromosomes from common peak regions,
# since the ChIP-seq samples to be normalized are associated with different
# genders.
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
norm <- normalize(H3K27Ac, 4:8, 9:13, common.peak.regions = autosome)
# Inspect the normalization effects.
attributes(norm)[5:8]
plot(attr(norm, "MA.cor"), symbreaks = TRUE, margins = c(8, 8))
MAplot(norm[[4]], norm[[5]], norm[[9]], norm[[10]],
       main = "GM12890_rep1 vs. GM12891_rep1")
abline(h = 0, lwd = 2, lty = 5)


## Alternatively, apply MA normalization first within each cell line, and
## then normalize across cell lines. In practice, this strategy is more
## recommended than the aforementioned one.
# Normalize samples separately for each cell line.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
# Construct separately a bioCond object for each cell line, and perform MA
# normalization on the resulting bioConds. Genomic intervals in sex
# chromosomes are not allowed to be common peak regions, since the cell
# lines are from different genders.
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Inspect the normalization effects.
attributes(conds)
plot(attr(conds, "MA.cor"), symbreaks = TRUE, margins = c(8, 8))
MAplot(conds[[1]], conds[[2]], main = "GM12890 vs. GM12891")
abline(h = 0, lwd = 2, lty = 5)


## Normalize directly the whole set of ChIP-seq samples by their size
## factors.
# Use only the genomic intervals that are occupied by all the ChIP-seq
# samples to be normalized to estimate the size factors.
norm <- normalizeBySizeFactors(H3K27Ac, 4:8,
                               subset = apply(H3K27Ac[9:13], 1, all))
# Inspect the normalization effects.
attr(norm, "size.factor")
MAplot(norm[[4]], norm[[5]], norm[[9]], norm[[10]],
       main = "GM12890_rep1 vs. GM12891_rep1")
abline(h = 0, lwd = 2, lty = 5)


## Alternatively, perform the normalization first within each cell line, and
## then normalize across cell lines. In practice, this strategy is more
## recommended than the aforementioned one.
# Normalize samples separately for each cell line.
norm <- normalizeBySizeFactors(H3K27Ac, 4)
norm <- normalizeBySizeFactors(norm, 5:6,
                               subset = apply(norm[10:11], 1, all))
norm <- normalizeBySizeFactors(norm, 7:8,
                               subset = apply(norm[12:13], 1, all))
# Construct separately a bioCond object for each cell line, and normalize
# the resulting bioConds by their size factors.
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
conds <- normBioCondBySizeFactors(conds)
# Inspect the normalization effects.
attr(conds, "size.factor")
MAplot(conds[[1]], conds[[2]], main = "GM12890 vs. GM12891")
abline(h = 0, lwd = 2, lty = 5)


## Construct a bioCond object for the GM12891 cell line.
# Apply MA normalization to the ChIP-seq samples of GM12891.
norm <- normalize(H3K27Ac, 5:6, 10:11)
# Call the constructor and optionally attach some meta information to the
# resulting bioCond, such as the coordinates of genomic intervals.
GM12891 <- bioCond(norm[5:6], norm[10:11], name = "GM12891",
                   meta.info = norm[1:3])
# Alternatively, you may assign different weights to the replicate samples
# for estimating the mean signal intensities of genomic intervals in this
# cell line. Here the weight of the 2nd sample is reduced to half the weight
# of the 1st one.
GM12891_2 <- bioCond(norm[5:6], norm[10:11], name = "GM12891",
                     weight = c(1, 0.5))
# Equivalently, you can achieve the same effect by setting the strMatrix
# parameter.
GM12891_3 <- bioCond(norm[5:6], norm[10:11], name = "GM12891",
                     strMatrix = list(diag(c(1, 2))))


## Set the weights of replicate ChIP-seq samples in a bioCond.
# Construct a bioCond object for the GM12891 cell line. By default, all the
# ChIP-seq samples belonging to the bioCond have the same weight for
# estimating the mean signal intensities of genomic intervals in the cell
# line.
norm <- normalize(H3K27Ac, 5:6, 10:11)
GM12891 <- bioCond(norm[5:6], norm[10:11], name = "GM12891")
# Now we set the weight of the 2nd sample to half of the 1st one.
GM12891_2 <- setWeight(GM12891, weight = c(1, 0.5))
# Equivalently, you can achieve the same effect by setting the strMatrix
# parameter.
GM12891_3 <- setWeight(GM12891, strMatrix = list(diag(c(1, 2))))


## Construct two bioConds comprised of the male and female individuals,
## respectively.
# First, normalize ChIP-seq samples separately for each individual (i.e.,
# cell line).
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
# Then, construct separately a bioCond for each individual, and perform MA
# normalization on the resulting bioConds. Genomic intervals in sex
# chromosomes are not allowed to be common peak regions, since the
# individuals are of different genders.
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Finally, group individuals into bioConds based on their genders.
female <- cmbBioCond(conds[c(1, 3)], name = "female")
male <- cmbBioCond(conds[2], name = "male")
summary(female)
summary(male)


## Fit a mean-variance curve treating each cell line (i.e., individual) as a
## biological condition.
# Perform the MA normalization and construct bioConds to represent cell
# lines.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Compare variance ratio factor between cell lines.
varRatio(conds$GM12891, conds$GM12892)
# The comparison is only possible when both of the bioConds to be compared
# have replicate samples.
varRatio(conds$GM12891, conds$GM12890)
# Deduce the relative variance ratio factors of a set of bioConds by making
# pairwise comparisons among them.
estimateVarRatio(conds)
estimateVarRatio(conds, base.cond = "GM12891")
# Variations in ChIP-seq signals across biological replicates of a cell line
# are generally of a low level, and their relationship with the mean signal
# intensities is expected to be well modeled by the presumed parametric
# form.
conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE)
summary(conds[[1]])
plotMeanVarCurve(conds, subset = "occupied")
# \dontrun
# Sometimes the parametric fitting algorithm cannot automatically deduce
# proper starting values for estimating the coefficients.
fitMeanVarCurve(conds[3], method = "parametric", occupy.only = TRUE)
# In such cases, explicitly specify the initial values of the coefficients.
fitMeanVarCurve(conds[3], method = "parametric", occupy.only = TRUE,
                init.coef = c(0.1, 10))


## Fit a mean-variance curve treating each gender as a biological condition,
## and each individual a replicate.
# First perform the MA normalization and construct bioConds to represent
# individuals.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Group individuals into bioConds based on their genders.
female <- cmbBioCond(conds[c(1, 3)], name = "female")
male <- cmbBioCond(conds[2], name = "male")
# The dependence of variance of ChIP-seq signal intensity across individuals
# on the mean signal intensity is not as regular as in the case for modeling
# biological replicates of cell lines. Better use the local regression to
# adaptively capture the mean-variance trend.
genders <- list(female = female, male = male)
genders1 <- fitMeanVarCurve(genders, method = "local", occupy.only = TRUE)
genders2 <- fitMeanVarCurve(genders, method = "local", occupy.only = FALSE)
# Suppose the local regression is performed using only the occupied genomic
# intervals as input. Good chances are that the extrapolation algorithm
# implemented in the regression method will produce over-estimated variances
# for the non-occupied intervals.
plotMeanVarCurve(genders1, subset = "all")
plotMeanVarCurve(genders2, subset = "all")
plotMeanVarCurve(genders1, subset = "non-occupied")
plotMeanVarCurve(genders2, subset = "non-occupied")
# On the other hand, applying the local regression on all genomic intervals
# may considerably reduce the estimated number of prior degrees of freedom
# associated with the fitted mean-variance curve, as ChIP-seq signals in the
# non-occupied intervals are generally of less data regularity compared with
# those in the occupied intervals.
summary(genders1$female)
summary(genders2$female)
# To split the difference, fit the mean-variance curve on all genomic
# intervals and re-estimate the number of prior degrees of freedom using
# only the occupied intervals, which is also the most recommended strategy
# in practice.
genders3 <- estimatePriorDf(genders2, occupy.only = TRUE)
plotMeanVarCurve(genders3, subset = "all")
plotMeanVarCurve(genders3, subset = "non-occupied")
summary(genders3$female)


## Fit a mean-variance curve for the GM12892 cell line (i.e., individual)
## and set the number of prior degrees of freedom of the curve to Inf.
# Perform the MA normalization and construct a bioCond to represent GM12892.
norm <- normalize(H3K27Ac, 7:8, 12:13)
GM12892 <- bioCond(norm[7:8], norm[12:13], name = "GM12892")
# Variations in ChIP-seq signals across biological replicates of a cell line
# are generally of a low level, and typically their relationship with the
# mean signal intensities could be well modeled by the presumed parametric
# form.
GM12892 <- fitMeanVarCurve(list(GM12892), method = "parametric",
                           occupy.only = TRUE, init.coef = c(0.1, 10))[[1]]
# In the vast majority of cases for modeling biological replicates of cell
# lines, the associated variance structure is so regular that variances of
# individual genomic intervals could be reliably estimated by fully
# depending on the mean-variance curve.
GM12892_2 <- setPriorDf(list(GM12892), Inf, occupy.only = TRUE)[[1]]
# The resulting model makes few differences from the original one, though.
# This is because MAnorm2 will adaptively deduce a large number of prior
# degrees of freedom for the mean-variance curve if the underlying variance
# structure is of high regularity. In practice, we recommend leaving the
# specification of prior df to the estimation method implemented in MAnorm2
# all the time.
summary(GM12892)
summary(GM12892_2)


## Fit a mean-variance curve based on the GM12891 cell line and associate
## the resulting curve with the other two cell lines.
# Perform the MA normalization and construct bioConds to represent cell
# lines.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Fit a mean-variance curve using only the GM12891 bioCond.
conds[2] <- fitMeanVarCurve(conds[2], method = "parametric",
                            occupy.only = TRUE)
summary(conds[[2]])
plotMeanVarCurve(conds[2], subset = "occupied")
# Associate the resulting curve with the other two bioConds.
conds[c(1, 3)] <- extendMeanVarCurve(conds[c(1, 3)], conds[[2]],
                                     occupy.only = TRUE)
summary(conds[[1]])
summary(conds[[3]])
plotMeanVarCurve(conds[3], subset = "occupied")
# Re-estimate number of prior degrees of freedom using all the bioConds,
# though the estimation result doesn't change in this example. But note the
# change of variance ratio factor of the bioCond without replicates (i.e.,
# GM12890).
conds2 <- estimatePriorDf(conds, occupy.only = TRUE)
summary(conds2[[1]])


## Make a comparison between GM12891 and GM12892 cell lines.
# Perform MA normalization and construct bioConds to represent the two cell
# lines.
norm <- normalize(H3K27Ac, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Variations in ChIP-seq signals across biological replicates of a cell line
# are generally of a low level, and their relationship with the mean signal
# intensities is expected to be well modeled by the presumed parametric
# form.
conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE)
summary(conds[[1]])
plotMeanVarCurve(conds, subset = "occupied")
# Perform differential tests between the two cell lines.
res <- diffTest(conds[[1]], conds[[2]])
head(res)
MAplot(res, padj = 0.001)
abline(h = 0, lwd = 2, lty = 5, col = "green3")


## Make a comparison between GM12891 and GM12892 cell lines using only their
## first replicates.
# Perform MA normalization and construct bioConds to represent the two cell
# lines.
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
norm <- normalize(H3K27Ac, c(5, 7), c(10, 12),
                  common.peak.regions = autosome)
conds <- list(GM12891 = bioCond(norm[5], norm[10], name = "GM12891"),
              GM12892 = bioCond(norm[7], norm[12], name = "GM12892"))
# Construct a "blind" bioCond that treats the two samples as replicates and
# fit a mean-variance curve accordingly. Only common peak regions of the two
# samples are considered to be occupied by the "blind" bioCond, and only
# these regions are used for fitting the mean-variance curve. This setting
# is for capturing underlying non-differential intervals as accurately as
# possible and avoiding over-estimation of prior variances (i.e., variances
# read from a mean-variance curve).
conds$blind <- bioCond(norm[c(5, 7)], norm[c(10, 12)], occupy.num = 2,
                       name = "blind")
conds[3] <- fitMeanVarCurve(conds[3], method = "parametric",
                            occupy.only = TRUE, init.coef = c(0.1, 10))
summary(conds[[3]])
plotMeanVarCurve(conds[3], subset = "occupied")
# Associate the resulting mean-variance curve with the two cell lines.
conds[1:2] <- extendMeanVarCurve(conds[1:2], conds[[3]])
summary(conds[[1]])
summary(conds[[2]])
# The following block of code is equivalent to the above two blocks.
conds$blind <- bioCond(norm[c(5, 7)], norm[c(10, 12)], occupy.num = 2,
                       name = "blind")
conds <- fitMeanVarCurve(conds, method = "parametric",
                         occupy.only = TRUE, init.coef = c(0.1, 10))
summary(conds[[1]])
summary(conds[[2]])
summary(conds[[3]])
plotMeanVarCurve(conds, subset = "occupied")
# Perform differential tests between the two cell lines.
res <- diffTest(conds[[1]], conds[[2]])
head(res)
MAplot(res, pval = 0.01)
abline(h = 0, lwd = 2, lty = 5, col = "green3")


## Make a comparison between the male and female genders by treating each
## individual (i.e., cell line) as a replicate.
# First perform the MA normalization and construct bioConds to represent
# individuals.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Group individuals into bioConds based on their genders.
female <- cmbBioCond(conds[c(1, 3)], name = "female")
male <- cmbBioCond(conds[2], name = "male")
# The dependence of variance of ChIP-seq signal intensity across individuals
# on the mean signal intensity is not as regular as in the case for modeling
# biological replicates of cell lines. Better use the local regression to
# adaptively capture the mean-variance trend.
genders <- list(female = female, male = male)
genders <- fitMeanVarCurve(genders, method = "local", occupy.only = FALSE)
genders <- estimatePriorDf(genders, occupy.only = TRUE)
summary(genders$female)
summary(genders$male)
plotMeanVarCurve(genders, subset = "all")
# Perform differential tests between the two genders.
res <- diffTest(genders[[1]], genders[[2]])
head(res)
MAplot(res, pval = 0.01)
abline(h = 0, lwd = 2, lty = 5, col = "green3")
# Examine the distribution of p-values in Y chromosome.
hist(res$pval[H3K27Ac$chrom == "chrY"], col = "red",
     main = "P-values in Y chromosome")


## Call differential genomic intervals among GM12890, GM12891 and GM12892
## cell lines.
# Perform MA normalization and construct bioConds to represent the cell
# lines.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Variations in ChIP-seq signals across biological replicates of a cell line
# are generally of a low level, and their relationship with the mean signal
# intensities is expected to be well modeled by the presumed parametric
# form.
conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE)
summary(conds[[1]])
plotMeanVarCurve(conds, subset = "occupied")
# Perform a moderated ANOVA on these cell lines.
res <- aovBioCond(conds)
head(res)
plot(res, padj = 1e-6)


## Call hypervariable and invariant genomic intervals across biological
## replicates of the GM12891 cell line.
## Call hypervariable and invariant genomic intervals across cell lines.
library(scales)
# Perform MA normalization and construct bioConds to represent cell lines.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
# Fit a mean-variance curve for the GM12891 cell line using the parametric
# method.
GM12891 <- fitMeanVarCurve(conds[2], method = "parametric",
                           occupy.only = TRUE)[[1]]
summary(GM12891)
plotMeanVarCurve(list(GM12891), subset = "occupied")
# Assess the observed variances of ChIP-seq signal intensities in GM12891.
res <- varTestBioCond(GM12891)
head(res)
# Inspect only the test results of occupied genomic intervals.
res <- res[GM12891$occupancy, ]
res$padj <- p.adjust(res$pval, method = "BH")
plot(res, padj = 0.2, col = alpha(c("black", "red"), c(0.04, 0.5)))
# Normalize the cell lines.
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Combine the cell lines into a single bioCond and use local regression to
# adaptively capture the mean-variance trend. Only genomic intervals that
# are occupied by each of the cell lines are considered to be occupied by
# the combined bioCond, which is for avoiding over-estimation of the prior
# variances.
LCLs <- cmbBioCond(conds, occupy.num = length(conds),
                   name = "lymphoblastoid_cell_lines")
LCLs <- fitMeanVarCurve(list(LCLs), method = "local",
                        occupy.only = FALSE)[[1]]
LCLs <- estimatePriorDf(list(LCLs), occupy.only = TRUE)[[1]]
summary(LCLs)
plotMeanVarCurve(list(LCLs), subset = "all")
# Assess the observed variances of ChIP-seq signal intensities across these
# cell lines.
res <- varTestBioCond(LCLs)
head(res)
plot(res, pval = 0.01, col = alpha(c("black", "red"), c(0.04, 0.5)))
# Non-occupied intervals are occupied by some of the cell lines but not all
# of them. These intervals tend to be more variable across the cell lines
# and more significant in the tests than occupied intervals.
f <- !(LCLs$occupancy)
wilcox.test(res$fold.change[f], res$fold.change[!f],
            alternative = "greater")
wilcox.test(res$pval[f], res$pval[!f], alternative = "less")
# Intervals in Y chromosome tend to be more variable across the cell lines
# and more significant in the tests than the other intervals, since the cell
# lines are of different genders.
f <- H3K27Ac$chrom == "chrY"
wilcox.test(res$fold.change[f], res$fold.change[!f],
            alternative = "greater")
wilcox.test(res$pval[f], res$pval[!f], alternative = "less")


## Cluster a set of ChIP-seq samples from different cell lines (i.e.,
## individuals).
# Perform MA normalization and construct a bioCond.
norm <- normalize(H3K27Ac, 4:8, 9:13)
cond <- bioCond(norm[4:8], norm[9:13], name = "all")
# Fit a mean-variance curve.
cond <- fitMeanVarCurve(list(cond), method = "local",
                        occupy.only = FALSE)[[1]]
plotMeanVarCurve(list(cond), subset = "all")
# Measure the distance between each pair of samples and accordingly perform
# a hierarchical clustering. Note that biological replicates of each cell
# line are clustered together.
d1 <- distBioCond(cond, method = "prior")
plot(hclust(d1, method = "average"), hang = -1)
# Measure the distances using only hypervariable genomic intervals. Note the
# change of scale of the distances.
res <- varTestBioCond(cond)
f <- res$fold.change > 1 & res$pval < 0.05
d2 <- distBioCond(cond, subset = f, method = "prior")
plot(hclust(d2, method = "average"), hang = -1)
# Apply a variance-stabilizing transformation and associate a constant
# function with the resulting bioCond as its mean-variance curve.
vst_cond <- vstBioCond(cond)
vst_cond <- setMeanVarCurve(list(vst_cond), function(x)
                            rep_len(1, length(x)), occupy.only = FALSE,
                            method = "constant prior")[[1]]
plotMeanVarCurve(list(vst_cond), subset = "all")
# Repeat the clustering analyses on the VSTed bioCond.
d3 <- distBioCond(vst_cond, method = "none")
plot(hclust(d3, method = "average"), hang = -1)
res <- varTestBioCond(vst_cond)
f <- res$fold.change > 1 & res$pval < 0.05
d4 <- distBioCond(vst_cond, subset = f, method = "none")
plot(hclust(d4, method = "average"), hang = -1)


## Cluster a set of individuals.
# Perform MA normalization and construct bioConds to represent individuals.
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
conds <- normBioCond(conds)
# Group the individuals into a single bioCond and fit a mean-variance curve
# for it.
cond <- cmbBioCond(conds, name = "all")
cond <- fitMeanVarCurve(list(cond), method = "local",
                        occupy.only = FALSE)[[1]]
plotMeanVarCurve(list(cond), subset = "all")
# Measure the distance between each pair of individuals and accordingly
# perform a hierarchical clustering. Note that GM12891 and GM12892 are
# actually a couple and they are clustered together.
d1 <- distBioCond(cond, method = "prior")
plot(hclust(d1, method = "average"), hang = -1)
# Measure the distances using only hypervariable genomic intervals. Note the
# change of scale of the distances.
res <- varTestBioCond(cond)
f <- res$fold.change > 1 & res$pval < 0.05
d2 <- distBioCond(cond, subset = f, method = "prior")
plot(hclust(d2, method = "average"), hang = -1)
# Apply a variance-stabilizing transformation and associate a constant
# function with the resulting bioCond as its mean-variance curve.
vst_cond <- vstBioCond(cond)
vst_cond <- setMeanVarCurve(list(vst_cond), function(x)
                            rep_len(1, length(x)), occupy.only = FALSE,
                            method = "constant prior")[[1]]
plotMeanVarCurve(list(vst_cond), subset = "all")
# Repeat the clustering analyses on the VSTed bioCond.
d3 <- distBioCond(vst_cond, method = "none")
plot(hclust(d3, method = "average"), hang = -1)
res <- varTestBioCond(vst_cond)
f <- res$fold.change > 1 & res$pval < 0.05
d4 <- distBioCond(vst_cond, subset = f, method = "none")
plot(hclust(d4, method = "average"), hang = -1)


## Perform differential analysis on bioConds that have gone through a
## variance-stabilizing transformation.
# Perform MA normalization and construct bioConds to represent cell lines
# (i.e., individuals).
norm <- normalize(H3K27Ac, 4, 9)
norm <- normalize(norm, 5:6, 10:11)
norm <- normalize(norm, 7:8, 12:13)
conds <- list(GM12890 = bioCond(norm[4], norm[9], name = "GM12890"),
              GM12891 = bioCond(norm[5:6], norm[10:11], name = "GM12891"),
              GM12892 = bioCond(norm[7:8], norm[12:13], name = "GM12892"))
autosome <- !(H3K27Ac$chrom %in% c("chrX", "chrY"))
conds <- normBioCond(conds, common.peak.regions = autosome)
# Fit a mean-variance curve.
conds <- fitMeanVarCurve(conds, method = "parametric", occupy.only = TRUE)
plotMeanVarCurve(conds, subset = "occupied")
# Apply a variance-stabilizing transformation.
vst_conds <- list(GM12890 = vstBioCond(conds$GM12890))
vst.func <- attr(vst_conds$GM12890, "vst.func")
temp <- matrix(vst.func(as.numeric(conds$GM12891$norm.signal)),
               nrow = nrow(norm))
vst_conds$GM12891 <- bioCond(temp, norm[10:11], name = "GM12891")
temp <- matrix(vst.func(as.numeric(conds$GM12892$norm.signal)),
               nrow = nrow(norm))
vst_conds$GM12892 <- bioCond(temp, norm[12:13], name = "GM12892")
# Associate a constant function with the resulting bioConds as their
# mean-variance curve.
vst_conds <- setMeanVarCurve(vst_conds, function(x) rep_len(1, length(x)),
                             occupy.only = TRUE, method = "constant prior")
plotMeanVarCurve(vst_conds, subset = "occupied")
# Make a comparison between GM12891 and GM12892.
res1 <- diffTest(conds$GM12891, conds$GM12892)
res2 <- diffTest(vst_conds$GM12891, vst_conds$GM12892)
# Examine the consistency of analysis results between using ordinary and
# VSTed signal intensities. Here we map p-values together with observed
# directions of signal changes to the standard normal distribution.
z1 <- qnorm(res1$pval / 2)
z1[res1$Mval > 0] <- -z1[res1$Mval > 0]
z2 <- qnorm(res2$pval / 2)
z2[res2$Mval > 0] <- -z2[res2$Mval > 0]
plot(z1, z2, xlab = "Ordinary", ylab = "VSTed")
abline(a = 0, b = 1, lwd = 2, lty = 5, col = "red")
cor(z1, z2)
cor(z1, z2, method = "sp")
# Simultaneously compare GM12890, GM12891 and GM12892 cell lines.
res1 <- aovBioCond(conds)
res2 <- aovBioCond(vst_conds)
# Examine the consistency of analysis results between using ordinary and
# VSTed signal intensities by mapping p-values to the standard normal
# distribution.
z1 <- qnorm(res1$pval, lower.tail = FALSE)
z1[z1 == Inf] <- 39
z2 <- qnorm(res2$pval, lower.tail = FALSE)
z2[z2 == Inf] <- 39
plot(z1, z2, xlab = "Ordinary", ylab = "VSTed")
abline(a = 0, b = 1, lwd = 2, lty = 5, col = "red")
cor(z1, z2)
cor(z1, z2, method = "sp")



## For performing statistical simulation to evaluate the p-value distribution
## under the null model of HyperChIP.
set.seed(1)


## Fixed parameter settings.
n <- 50000
A_lower <- 2
A_upper <- 9
f_MVC <- function(mu) {
    0.15 + 6 / 2^mu
}


## Random generators.

# For generating signal intensities
# based on normal distributions.
randomNormal <- function(mu, t, m) {
    num <- length(mu)
    sd <- sqrt(t)
    vapply(1:m, function(k) {
        rnorm(num, mean = mu, sd = sd)
    }, numeric(num))
}

# For generating a matrix of "normalized" signal intensities.
randomTable <- function(m = 36, gamma = 0.75, d0 = 28) {
    mu <- runif(n, A_lower, A_upper)
    t <- gamma * f_MVC(mu) * d0 / rchisq(n, d0)
    randomNormal(mu, t, m)
}


## Utility for applying HyperChIP.
library(MAnorm2)
hypervariableAnalysis <- function(d) {
    # Construct a bioCond.
    cond <- bioCond(d)

    # Fit a mean-variance curve.
    cond <- fitMeanVarCurve(list(cond))[[1]]

    # Apply the parameter estimation framework of HyperChIP.
    # All genomic regions rather than low-intensity ones are used.
    cond <- estParamHyperChIP(cond, prob = 1)

    # Perform statistical tests.
    res <- varTestBioCond(cond)

    # Get one-sided p-values for identifying HVRs.
    df <- attr(res, "df")
    pf(res$fold.change, df[1], df[2], lower.tail = FALSE)
}


## Perform statistical simulation.




# Scenario 1: default parameter settings.
pval_1 <- hypervariableAnalysis(randomTable())

# Scenario 2: decrease the total number of samples.
pval_2 <- hypervariableAnalysis(randomTable(m = 36 / 2))

# Scenario 3: increase the global signal variability.
pval_3 <- hypervariableAnalysis(randomTable(gamma = 0.75 * 2))

# Scenario 4: decrease the number of prior degrees of freedom.
pval_4 <- hypervariableAnalysis(randomTable(d0 = 28 / 2))


## Inspect the distribution of p-values.

# Utility for plotting histograms.
plotHistogram <- function(pval, main) {
    hist(pval, freq = FALSE, col = "#E41A1C", main = main,
         xlab = expression(paste(italic(P), "-value")))
    abline(h = 1, lty = 2)
}

pdf("HyperChIP_statistical_simulation.pdf", width = 3.3, height = 3.3, pointsize = 9)
par(mar = c(4, 4, 2, 2), cex.main = 1, font.main = 1)
plotHistogram(pval_1, main = "Default parameter settings")
plotHistogram(pval_2, main = expression(italic(m) == 36 / 2))
plotHistogram(pval_3, main = expression(gamma == 0.75 %*% 2))
plotHistogram(pval_4, main = expression(italic(d[0]) == 28 / 2))
dev.off()



## For performing statistical simulation.
set.seed(100)


## Fixed parameter settings.

n <- 30000
m1 <- 3
gamma_1 <- 1
d0 <- 10
f_MVC <- function(mu) {
    5e-3 + 20 / 2^mu
}

A_lower <- 3
A_upper <- 11
v_M <- 4


## Random generators.

# For generating signal intensities
# based on the normal distribution.
randomNormal <- function(mu, t, m) {
    num <- length(mu)
    sd <- sqrt(t)
    vapply(1:m, function(k) {
        rnorm(num, mean = mu, sd = sd)
    }, numeric(num))
}

# For generating a matrix of "normalized" signal intensities
# as well as the associated true labels.
randomTable <- function(diff_prop = 0.1, m2 = m1, gamma_2 = gamma_1) {
    # For differential genomic intervals.
    num <- max(2, round(n * diff_prop))
    A <- runif(num, A_lower, A_upper)
    M <- rnorm(num, mean = 0, sd = sqrt(v_M))
    label <- rep(-1, num)
    label[M > 0] <- 1
    mu_1 <- A - M / 2
    mu_2 <- A + M / 2
    t1 <- gamma_1 * f_MVC(mu_1) * d0 / rchisq(num, d0)
    t2 <- gamma_2 * f_MVC(mu_2) * d0 / rchisq(num, d0)
    ret <- cbind(randomNormal(mu_1, t1, m1),
                 randomNormal(mu_2, t2, m2))

    # For non-differential genomic intervals.
    num <- max(2, n - num)
    label <- c(label, rep(0, num))
    mu_1 <- runif(num, A_lower, A_upper)
    mu_2 <- mu_1
    temp <- f_MVC(mu_1) * d0 / rchisq(num, d0)
    t1 <- gamma_1 * temp
    t2 <- gamma_2 * temp
    ret <- rbind(ret, cbind(randomNormal(mu_1, t1, m1),
                            randomNormal(mu_2, t2, m2)))

    # Assemble the results.
    colnames(ret) <- c(paste("s_1", 1:m1, sep = "_"),
                       paste("s_2", 1:m2, sep = "_"))
    cbind(data.frame(label = label), as.data.frame(ret))
}


## Utility for applying MAnorm2.
library(MAnorm2)
diffAnalysis <- function(d) {
    # Construct a bioCond for each group of samples.
    conds <- list(c1 = bioCond(d[2:(1 + m1)], name = "c1"),
                  c2 = bioCond(d[(2 + m1):ncol(d)], name = "c2"))

    # Fit a mean-variance curve.
    conds <- fitMeanVarCurve(conds, method = "parametric", init.coef = c(0.1, 10))

    # Perform statistical tests for identifying differential signals.
    diffTest(conds[[1]], conds[[2]])
}


## Perform statistical simulation.

# Scenario 1: default parameter settings.
d1 <- randomTable()
res1 <- diffAnalysis(d1)

# Scenario 2: different global within-group variability.
d2 <- randomTable(gamma_2 = 3 * gamma_1)
res2 <- diffAnalysis(d2)

# Scenario 3: different group sizes.
d3 <- randomTable(m2 = 3 * m1)
res3 <- diffAnalysis(d3)

# Scenario 4: increase the proportion of differential genomic intervals.
d4 <- randomTable(diff_prop = 0.3)
res4 <- diffAnalysis(d4)


## Inspect the distribution of p-values for non-differential genomic intervals.

# Utility for plotting histograms.
plotHistogram <- function(x, main) {
    hist(x, freq = FALSE, col = "#E41A1C", main = main,
         xlab = "P-values for non-differential intervals")
    abline(h = 1, lty = 2)
}

pdf("MAnorm2_statistical_simulation.pdf", width = 3.3, height = 3.3, pointsize = 9)
par(mar = c(4, 4, 2, 2), cex.main = 1, font.main = 1)
plotHistogram(res1$pval[d1$label == 0], main = "Default parameter settings")
plotHistogram(res2$pval[d2$label == 0], main = expression(gamma[2] == 3 %.% gamma[1]))
plotHistogram(res3$pval[d3$label == 0], main = expression(m[2] == 3 %.% m[1]))
plotHistogram(res4$pval[d4$label == 0], main = "Proportion of differential intervals: 30%")
dev.off()




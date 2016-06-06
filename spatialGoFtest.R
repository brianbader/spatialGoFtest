## Code written as of 6-6-2016
## Function to generate a 'spatial' goodness-of-fit test under unspecified dependence
## For a given fitted multivariate (spatial) model, transform the marginal distributions 
## to uniform scale.
##
## Here, the test statistic is the maximum value of the one sample Kolmogorov-Smirnov 
## test statistics performed for each (uniform) marginal distribution.
##
## The function generates monte carlo versions of the testing statistic, preserving 
## the dependence structure (see reference). The null hypothesis is that the 
## marginal distributions at each site follows the Uniform(0, 1) distribution. If the 
## model is correctly specified, the transformation (via the cumulative distribution 
## function) should yield approximately uniformly distributed data.
##
## See example below.
##
## References: 
##
## Iman, R.L., and W.J. Conover. (1982). A Distribution-Free Approach to Inducing Rank 
## Correlation Among Input Variables. Communications in Statistics--Volume B, Simulation 
## and Computation, 11(3), 311--334.
##
## Also see function 'simulateMvMatrix' in package EnvStats.
##
##
##
##
## Function input:
## umat: Matrix of transformed data (to the uniform scale). Each column represents a site 
## or grouping, each row represents an observation.
##
## bootnum: The number of bootstrap replicates of the null distribution to take.
## 
## plot: Whether to plot the empirical CDF functions of each site/group.
## 
## allowParallel: Whether to allow parallel computations of the bootstrap replicates (using 
## package 'parallel').
## 
## numCores: If 'allowParallel = TRUE', the number of cores to use in parallel.
## 
##
##
##
## Example: Generate non-stationary bi-variate normal data, but 
## choose to fit with a misspecified stationary model...
## usePackage(c("MASS"))
## set.seed(7)
## Sigma <- matrix(c(10,3,3,2),2,2)
## num <- 100
## mat <- matrix(0, nrow = num, ncol = 2)
##
## for(i in 1:num) {
##   m1 <- 0 + i / 60
##   m2 <- 1 + i / 70
##   mat[i, ] <- mvrnorm(n = 1, c(m1, m2), Sigma)
## }
#
## mat[, 1] <- pnorm(mat[, 1], mean = 0, sd = sqrt(Sigma[1, 1]))
## mat[, 2] <- pnorm(mat[, 2], mean = 1, sd = sqrt(Sigma[2, 2]))
## 
## spatialGoFtest(mat, bootnum = 500, plot = TRUE)
##
##
##
##
## Function that checks if a package is installed and if not, installs it
## Then loads all requested packages using 'require' command
usePackage <- function(p) {
  newPackages <- p[!(p %in% installed.packages()[, "Package"])]
  if(length(newPackages))
    install.packages(newPackages, dependencies = TRUE)
  cat("Packages successfully loaded:\n")
  sapply(p, require, character.only = TRUE, quietly = TRUE)
}

packageList <- c("EnvStats", "parallel", "ggplot2")
usePackage(packageList)


spatialGoFtest <- function(umat, bootnum, plot = TRUE, allowParallel = FALSE, numCores = 1) {
  if(any(umat > 1) | any(umat < 0))
    stop("Data must be on uniform (0, 1) scale")
  cormat <- cor(umat, method = "spearman")
  stat <- max(apply(umat, 2, function(x) ks.test(x, "punif")$statistic))
  nr <- nrow(umat)
  nc <- ncol(umat)
  gen.list <- rep(list(list(min = 0, max = 1)), nc)
  bootsim <- function(umat, cormat, nr, nc, gen.list) {
    umat.new <- simulateMvMatrix(n = nr, cor.mat = cormat, 
                                 distributions = rep("unif", nc), param.list = gen.list)
    max(apply(umat.new, 2, function(x) ks.test(x, "punif")$statistic))
  }
  if(allowParallel == TRUE) {
    cl <- makeCluster(numCores)
    clusterExport(cl, c('simulateMvMatrix'))
    fun <- function(cl) {
      parSapply(cl, 1:bootnum, function(i,...) {bootsim(umat, cormat, nr, nc, gen.list)})
    }
    teststat <- fun(cl)
    stopCluster(cl)
  } else {
    teststat <- replicate(bootnum, bootsim(umat, cormat, nr, nc, gen.list))
  }
  if(plot) {
    umat <- as.data.frame(umat)
    umat <- reshape(umat, varying = list(names(umat)), direction = "long", timevar = "Site")
    umat$Site <- as.factor(umat$Site)
    ks_plot <- ggplot(umat, aes(V1, colour = Site)) + 
      stat_ecdf() +
      stat_function(fun = punif, args = list(0, 1), colour = "black") + 
      labs(x = "Empirical", y = "Theoretical", title = "Empirical vs. Theoretical CDF") + 
      theme(legend.position = "bottom")
    print(ks_plot)
  }
  p <- (sum(teststat > stat) + 1) / (bootnum + 2)
  out <- list(p, stat)
  names(out) <- c("p.value", "test.stat")
  out
}
# simulation of high dimensional case
# desire: factorial design, interaction term

n <- 40 # samples, in four groups
p <- 1000 # number of features, say genes
pi0 <- .95 # proportion of null genes

# suppose a 2x2 design
# condition: control or treated
# group: A or B
condition <- factor(rep(rep(c("ctr","trt"),each=n/4),2))
group <- factor(rep(c("A","B"),each=n/2))

# suppose there is a shift across the two groups which is not of interest
grp.eff <- 1

# we are interested in differences in the condition effect across group
cond.grp.interx.eff <- 2

# some continuous data:
# in classic genomics style, rows are features, columns are samples
x <- matrix(rnorm(n*p),ncol=n)
# add the group shift
x[,group=="B"] <- x[,group=="B"] + grp.eff
# add the interaction for (1-pi0) of the genes
idx <- condition=="trt" & group=="B"
x[1:(p-pi0*p),idx] <- x[1:(p-pi0*p),idx] + cond.grp.interx.eff

# plot the data
cols <- colorRampPalette(c("white","darkblue"))(99)
image2 <- function(x) image(t(x[nrow(x):1,]),col=cols)
image2(x)

# the sample groups
table(condition, group)

# define a function for later
shuf <- function(l) sample(which(l),sum(l))
# need this package for row-wise Wilcoxon,
# install with:
# install.packages("BiocManager")
# BiocManager::install("siggenes")
library(siggenes)

# function to compute row-wise wilcoxon on matrix 'z'
# 'cl' is a binary vector (0/1) for the classes
computeWilcoxon <- function(z, cl) {
  zzz <- capture.output({w <- rowWilcoxon(z, cl)})
  w
}

# a function to compute a proposed test statistic for the interaction
computeTestStat <- function(x, npseudo, null=FALSE) {
  xa <- x[,group=="A"]
  xb <- x[,group=="B"]
  # the condition vector for 'xa' or 'xb'
  cond <- condition[group=="A"]
  # the group vector once we compute pair differences
  grp <- group[condition=="ctr"]
  if (null) {
    grp <- grp[sample(n/2)] # for permuting across group variable
  }
  # siggenes wilcoxon requires binary vector
  cl <- as.integer(grp) - 1
  stats <- replicate(npseudo, {
    # pseudo-pairing and compute differences
    da <- xa[,shuf(cond=="trt")] - xa[,shuf(cond=="ctr")]
    db <- xb[,shuf(cond=="trt")] - xb[,shuf(cond=="ctr")]
    # test the differences across grouping variable
    deltas <- cbind(da, db)
    w <- computeWilcoxon(deltas, cl)
  })
  # average over the pseudo-pairing
  rowMeans(stats)
}

# this is our observed test statistics
stat <- computeTestStat(x, npseudo=20)
plot(stat)

### 2 ways to approximate the null distribution

# first approach, permute the pairs across group,
# but always ensure that the pairs are within an original group.
# perform one pseudo-pairing per permutation
# (instead of many pseudo-pairings and averaging)
nperms <- 30
nulls <- replicate(nperms, computeTestStat(x, npseudo=1, null=TRUE))
plot(c(stat, as.vector(nulls[,1:5])),
     col=rep(c("black","brown"),each=p))

# plot the nulls in the original data and the estimated null distribution
plot(density(stat[(p-pi0*p+1):p]))
lines(density(nulls), col="brown")

### alternatively, ignore the grouping for permutations

nperms <- 30
nulls <- replicate(nperms, {
  computeTestStat(x[,sample(n)], npseudo=20)
  })
plot(c(stat, as.vector(nulls[,1:5])),
     col=rep(c("black","brown"),each=p))

# plot the nulls in the original data and the estimated null distribution
plot(density(stat[(p-pi0*p+1):p]))
lines(density(nulls), col="brown")

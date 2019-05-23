# simulation of high dimensional case
# desire: factorial design, interaction term

n <- 40 # in four groups
p <- 1000
pi0 <- .95

condition <- factor(rep(rep(c("ctr","trt"),each=n/4),2))
group <- factor(rep(c("A","B"),each=n/2))

grp.eff <- 1
cond.grp.interx.eff <- 2

x <- matrix(rnorm(n*p),ncol=n)
x[,group=="B"] <- x[,group=="B"] + grp.eff
idx <- condition=="trt" & group=="B"
x[1:(p-pi0*p),idx] <- x[1:(p-pi0*p),idx] + cond.grp.interx.eff

cols <- colorRampPalette(c("white","darkblue"))(99)
image2 <- function(x) image(t(x[nrow(x):1,]),col=cols)

image2(x)

table(condition,group)

shuf <- function(l) sample(which(l),sum(l))
library(siggenes)

nreps <- 20

computeWilcoxon <- function(z, cl) {
  zzz <- capture.output({w <- rowWilcoxon(z, cl)})
  w
}

computeTestStat <- function(x, npseudo, null=FALSE) {
  xa <- x[,group=="A"]
  xb <- x[,group=="B"]
  cond <- condition[group=="A"]
  grp <- group[condition=="ctr"]
  if (null) {
    grp <- grp[sample(n/2)] # for permuting across group
  }
  cl <- as.integer(grp) - 1
  stats <- replicate(npseudo, {
    # pseudo-pairing
    da <- xa[,shuf(cond=="trt")] - xa[,shuf(cond=="ctr")]
    db <- xb[,shuf(cond=="trt")] - xb[,shuf(cond=="ctr")]
    deltas <- cbind(da, db)
    w <- computeWilcoxon(deltas, cl)
  })
  # average over the pseudo-pairing
  rowMeans(stats)
}

stat <- computeTestStat(x, npseudo=20)
plot(stat)

nperms <- 5 # this is 30 in our method, the higher the better
nulls <- replicate(nperms, computeTestStat(x, npseudo=1, null=TRUE))
plot(c(stat, as.vector(nulls)),
     col=rep(c("black","brown"),each=p))

plot(density(stat[(p-pi0*p+1):p]))
lines(density(nulls), col="brown")

### alternatively, ignore the grouping for permutations

nperms <- 20
nulls <- replicate(nperms, {
  computeTestStat(x[,sample(n)], npseudo=20)
  })
plot(c(stat, as.vector(nulls)),
     col=rep(c("black","brown"),each=p))

plot(density(stat[(p-pi0*p+1):p]))
lines(density(nulls), col="brown")

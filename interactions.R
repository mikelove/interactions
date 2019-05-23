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

computeTestStat <- function(x, nreps, null=FALSE) {
  xa <- x[,group=="A"]
  xb <- x[,group=="B"]
  cond <- condition[group=="A"]
  grp <- group[condition=="ctr"]
  if (null) {
    grp <- grp[sample(n/2)] # for permuting across group
  }
  cl <- as.integer(grp) - 1
  stats <- replicate(nreps, {
    # pseudo-pairing
    da <- xa[,shuf(cond=="trt")] - xa[,shuf(cond=="ctr")]
    db <- xb[,shuf(cond=="trt")] - xb[,shuf(cond=="ctr")]
    deltas <- cbind(da, db)
    zzz <- capture.output({w <- rowWilcoxon(deltas, cl)})
    w
  })
  # average over the pseudo-pairing
  rowMeans(stats)
}

stat <- computeTestStat(x, nreps=20)
plot(stat)

nperms <- 20
nulls <- replicate(nperms, computeTestStat(x, nreps=1, null=TRUE))
plot(c(stat, as.vector(nulls[,1:5])),
     col=rep(c("black","brown"),each=p))

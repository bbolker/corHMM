library(devtools)
load_all("~/R/pkgs/corHMM")

data(primates)
phy <- multi2di(primates[[1]])
data <- primates[[2]]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1, return.devfun=TRUE)

## return deviance function
MK_3state$devfun()
ee <- environment(MK_3state$devfun)
p <- ee$p

nllfun <- function(p) {
  MK_3state$devfun(p=p)
}
names(p) <- parnames(nllfun) <- paste("p",1:4)

n
library(bbmle)
m1 <- mle2(minuslogl=nllfun, start=p, vecpar=TRUE)

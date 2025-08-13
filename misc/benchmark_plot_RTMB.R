library(tidyverse)
theme_set(theme_bw())

## exploring differences ...
dd <- readRDS("benchmark1.rds")
if (!is.data.frame(dd)) dd <- do.call(rbind, dd)
print(nrow(dd))
tail(dd)

with(dd, summary(RTMB_loglik-orig_loglik))
dd2 <- dd |>
  mutate(ldiff = RTMB_loglik-orig_loglik) |>
  arrange(ldiff) |>
  mutate(n = seq(n()))

filter(dd2, ldiff<0) |> View()

ggplot(dd2, aes(n, ldiff)) + geom_point()

## will have to hack corHMM further: store convergence codes from optimizer?
## explore: which cases are bad?
## worst for smallest/overparameterized cases (3 traits/20 taxa, 2 traits/41 taxa)

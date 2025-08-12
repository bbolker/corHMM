library(tidyverse)
theme_set(theme_bw())

## exploring differences ...
dd <- readRDS("benchmark1.rds")
if (!is.data.frame(dd)) dd <- do.call(rbind, dd)
print(nrow(dd))
tail(dd)

with(dd, summary(RTMB_loglik-orig_loglik))
dd |>
  mutate(ldiff = RTMB_loglik-orig_loglik) |>
  arrange(ldiff) |>
  View()

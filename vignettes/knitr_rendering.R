library(knitr)
opts_knit$set(base.dir = 'vignettes', fig.path = 'vignettes')
knit("vignettes/realModelExample.Rmd.orig", output = "vignettes/realModelExample.Rmd")

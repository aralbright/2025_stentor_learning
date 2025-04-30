library("sleuth")
library("tidyr")
library("ggplot2")
library("plyr")
library("dplyr")
library("UpSetR")
library("sva")

# load data
sample_id <- dir(file.path("/Volumes/albright_postdoc/deepa/exp2/20240129_kallisto_output/"))
kal_dirs <- file.path("/Volumes/albright_postdoc/deepa/exp2/20240129_kallisto_output/",sample_id)

# samples to onditions
s2c <- read.table(file.path("/Volumes/albright_postdoc/deepa/exp2/bin", "conditions.csv"), header = TRUE, sep=",")
s2c <- s2c[, colSums(is.na(s2c) | s2c == "") < nrow(s2c)]
s2c <- s2c[rowSums(is.na(s2c) | s2c == "") < ncol(s2c), ]
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# filter data 
s2c <- s2c[s2c$condition %in% c('ut_19', 'ut_f_10', 'ut_f_90'), ]

# prep sleuth object 
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# fit models
so <- sleuth_fit(so, ~ condition, "full")
so <- sleuth_fit(so, ~ 1, "reduced")
so <- sleuth_lrt(so, "reduced", "full")

# shiny to explore data
sleuth_live(so)

# save results table
# all q-values are equal to 1, therefore these controls are redundant
res <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

write.csv(res, file = "exp2_RNA_controlremoval_res.csv")



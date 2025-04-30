##########################################################################
# RNA Substrates of Memory and Learning (and Forgetting!) in Stentor coeruleus 
##########################################################################

library("sleuth")
library("tidyr")
library("ggplot2")
library("plyr")
library("dplyr")

# load data
sample_id <- dir(file.path("/Volumes/albright_postdoc/deepa/exp2/20240129_kallisto_output/"))
kal_dirs <- file.path("/Volumes/albright_postdoc/deepa/exp2/20240129_kallisto_output/",sample_id)

# samples to conditions
s2c <- read.table(file.path("/Volumes/albright_postdoc/deepa/exp2/bin", "conditions.csv"), header = TRUE, sep=",")
s2c <- s2c[, colSums(is.na(s2c) | s2c == "") < nrow(s2c)]
s2c <- s2c[rowSums(is.na(s2c) | s2c == "") < ncol(s2c), ]
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# remove redudant controls
s2c <- s2c[!s2c$condition %in% c('ut_f_10', 'ut_f_90'), ]


#### with time 0 

# prep sleuth object 
so_0 <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# fit models
so_0 <- sleuth_fit(so_0, ~ condition, "full")
so_0 <- sleuth_fit(so_0, ~ 1, "reduced")
so_0 <- sleuth_lrt(so_0, "reduced", "full")

# view results
res_0 <- sleuth_results(so_0, 'reduced:full', test_type = 'lrt')
write.csv(res_0, file = "exp2_rna_sleuth_lrt.csv")

res_0_obsnorm <- sleuth_to_matrix(so_0, 'obs_norm', 'est_counts')
write.csv(res_0_obsnorm, file = "exp2_rna_sleuth_obsnorm.csv")

#### withOUT time 0 
s2c <- s2c[!s2c$condition %in% c('t_0'), ]

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~ condition, "full")
so <- sleuth_fit(so, ~ 1, "reduced")
so <- sleuth_lrt(so, "reduced", "full")

sleuth_live(so)

res <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
write.csv(res, file = "exp2_rna_no0_sleuth_lrt.csv")

res_obsnorm <- sleuth_to_matrix(so, 'obs_norm', 'est_counts')
write.csv(res_obsnorm, file = "exp2_rna_no0_sleuth_obsnorm.csv")

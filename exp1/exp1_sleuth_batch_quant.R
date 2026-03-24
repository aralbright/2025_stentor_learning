##########################################################################
# RNA Substrates of Memory and Learning in Stentor coeruleus 
##########################################################################

library("sleuth")
library("tidyr")
library("ggplot2")
library("plyr")
library("dplyr")

sample_id <- dir(file.path("/Volumes/albright_postdoc/deepa/kallisto_output"))

kal_dirs <- file.path("/Volumes/albright_postdoc/deepa/kallisto_output",sample_id)

s2c <- read.table(file.path("/Volumes/albright_postdoc/deepa", "conditions.csv"), header = TRUE, sep=",")
#s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Extract the last three letters of each sample name
sample_last_three <- substr(s2c$sample, nchar(s2c$sample) - 2, nchar(s2c$sample))

# Match the last three letters of each path to the sample names in s2c and populate a new column with the matching paths
s2c$path <- sapply(sample_last_three, function(x) kal_dirs[grep(paste0(x, "$"), kal_dirs)])

# View the result
print(s2c)


#####################################
# 19 hour training 
#####################################

# filter data
samp_19 <- c("trained_19", "untrained_19")

s2c_19 <- filter(s2c, s2c$condition %in% samp_19) 

# prep sleuth object 
so_19 <- sleuth_prep(s2c_19, extra_bootstrap_summary = TRUE)

sleuth_live(so_19) # for exploration + PCA 

# remove outlier 
s2c_19 <- s2c_19[!s2c_19$sample == "R1c",]

# prep so 
so_19 <- sleuth_prep(s2c_19, extra_bootstrap_summary = TRUE)

sleuth_live(so_19) # for exploration + PCA

# want model set up certain way
so_19$sample_to_covariates$condition <- factor(so_19$sample_to_covariates$condition)
so_19$sample_to_covariates$condition <- relevel(so_19$sample_to_covariates$condition, ref = "untrained_19")

# check for batch effects 
so_19 <- sleuth_fit(so_19, ~1, 'reduced')
so_19 <- sleuth_fit(so_19, ~batch, 'batch')

so_19 <- sleuth_lrt(so_19, 'reduced', 'batch')
sleuth_batch_19 <- sleuth_results(so_19, 'reduced:batch', 'lrt', show_all = FALSE)

# no detected batch effects, fit condition model and wald test 
so_19 <- sleuth_fit(so_19, ~condition, 'full')
so_19 <- sleuth_wt(so_19, 'conditiontrained_19', which_model = 'full')


sleuth_19_wald <- sleuth_results(so_19, 'conditiontrained_19', 'wt', show_all = FALSE)
sleuth_19_sig <- dplyr::filter(sleuth_19_wald, qval <= 0.2)

write.csv(sleuth_19_wald, "/Volumes/albright_postdoc/deepa/exp1_sleuth_19_wald_v2.csv")

write.csv(sleuth_19_sig, "/Volumes/albright_postdoc/deepa/exp1_sleuth_19_wald_sig_v2.csv")

#####################################
# 4 hour training 
#####################################

# filter data
samp_4 <- c("trained_4", "untrained_4")

s2c_4 <- filter(s2c, s2c$condition %in% samp_4) 

# prep sleuth object 
so_4 <- sleuth_prep(s2c_4, extra_bootstrap_summary = TRUE)

#sleuth_live(so_4) # for exploration 

# remove outlier 
s2c_4 <- s2c_4[!s2c_4$sample == "B1a",]

# prep so 
so_4 <- sleuth_prep(s2c_4, extra_bootstrap_summary = TRUE)

#sleuth_live(so_4) # for exploration 

# want model set up certain way
so_4$sample_to_covariates$condition <- factor(so_4$sample_to_covariates$condition)
so_4$sample_to_covariates$condition <- relevel(so_4$sample_to_covariates$condition, ref = "untrained_4")

# check for batch effects 
so_4 <- sleuth_fit(so_4, ~1, 'reduced')
so_4 <- sleuth_fit(so_4, ~batch, 'batch')

so_4 <- sleuth_lrt(so_4, 'reduced', 'batch')
sleuth_batch_4 <- sleuth_results(so_4, 'reduced:batch', 'lrt', show_all = FALSE)

# no detected batch effects, fit condition model and wald test 
so_4 <- sleuth_fit(so_4, ~condition, 'full')
so_4 <- sleuth_wt(so_4, 'conditiontrained_4', which_model = 'full')

sleuth_4_wald <- sleuth_results(so_4, 'conditiontrained_4', 'wt', show_all = FALSE)
sleuth_4_sig <- dplyr::filter(sleuth_4_wald, qval <= 0.2)

write.csv(sleuth_4_wald, "/Volumes/albright_postdoc/deepa/exp1_sleuth_4_wald_v2.csv")

write.csv(sleuth_4_sig, "/Volumes/albright_postdoc/deepa/exp1_sleuth_4_wald_sig_v2.csv")



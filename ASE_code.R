##ASE analysis continuation in R. The major analysis is conduct a binomial test and fdr correct the p values

library(data.table)
#always setwd to avoid wierd errors
setwd("/scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/")
##this file is the product of ls *_wasp_counter > counts_files.txt and basically has all the needed files
counts_files <- read.delim("/scratch/aa9gj/wrk_backup/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/counts_files.txt", header=FALSE)
##read in all the files from the directory
counts_data <- list()
for (i in seq_along(counts_files$V1)) {
  counts_data[[i]] <- read.delim(file = counts_files$V1[i], header=TRUE)
}
##split them to their respective files
split_df_counts<-function(list){
  for (i in 1:length(list)){
    assign(paste0(counts_files$V1[i]), list[[i]], envir = .GlobalEnv)
  }
}
split_df_counts(counts_data)
##function to analyze the data. the x is one of those counter wasp files. notice it conducts a binomial test and then fdr correct it.
ASE_function<-function(x) {
  x$percent_ref <- x$refCount/x$totalCount
  x <- filter(x, percent_ref > 0.10 & percent_ref < 0.90)
  #x<-filter(x, variantID != ".")
  x<-filter(x, totalCount >= 20)
  for (i in 1:length(x$contig)) {
    temp = binom.test(x$altCount[i], x$totalCount[i], p = 0.5)
    x$binom_p[i] = temp$p.value
    x$binom_success[i]= temp$estimate
    x$fdr[i] <- p.adjust(x$binom_p[i], method = "fdr", n=length(x$binom_p))
  }
  x_filterd<-filter(x, fdr < 0.05)
  return(x_filterd)
}
##perform on all samples
my_ASE_results <- list()
for (i in seq_along(counts_files$V1)) {
  my_ASE_results[[i]] <- ASE_function(counts_data[[i]])
}
##split into their respective results files
split_df_ASE<-function(list){
  for (i in 1:length(list)){
    assign(paste0(counts_files$V1[i], "ASE_res"), list[[i]], envir = .GlobalEnv)
  }
}
split_df_ASE(my_ASE_results)

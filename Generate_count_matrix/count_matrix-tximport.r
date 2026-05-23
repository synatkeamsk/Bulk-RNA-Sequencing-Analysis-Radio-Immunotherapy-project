
#' Quantification aligned sample to count matrix 

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("readr")
BiocManager::install("DESeq2")

# Create sample metadata 
samples <- data.frame(
  sample = c("sample1","sample2","sample3", "sample4","sample5","sample6", "sample7","sample8","sample9", "sample10","sample11","sample12"),
  condition = c("1d_post_rt","3d_post_rt", "6d_post_rt", "9d_post_rt"),
  files = file.path("kallisto_results", c("sample1","sample2","sample3", "sample4","sample5","sample6", "sample7","sample8","sample9", "sample10","sample11","sample12"), "abundance.tsv")
)

samples <- data.frame(
  sample = paste0("sample", 1:12),
  files = file.path("kallisto_results", paste0("sample", 1:12), "abundance.tsv")
)

#' Run tximport (kallisto to gene counts)
library(tximport)
library(readr)

files <- setNames(samples$files, samples$sample)

txi <- tximport(
  files,
  type = "kallisto",
  txOut = FALSE   # IMPORTANT: gives gene-level counts
)


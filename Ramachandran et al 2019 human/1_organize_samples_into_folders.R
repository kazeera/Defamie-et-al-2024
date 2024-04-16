# Note - this script puts each sample into its own 10x folder
# ** Works on unzipped raw GSE data, no need to rerun
# we analysed 67,494 human cells from healthy (n = 5) and cirrhotic (n = 5) liver

# May 4, 2021
# Reformat gz files into individual folder for each sample run that can be read by the Read10X() function in R

# The files in each folder currently look like this:
# # [1] "GSM2759555_5wk-2-barcodes.tsv.gz" "GSM2759555_5wk-2-genes.tsv.gz"    "GSM2759555_5wk-2-matrix.mtx.gz"  

# But they need to be this:
# [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx" 

# Note: The read10x only works with .gz files now
# "genes.tsv.gz" is supposed to be called "features.tsv.gz" now 
library(dplyr)
library(kazutils)

# Set to directory with all files
setwd("Data/GSE136103_RAW")

# Get sample IDs
x <- list.files(pattern="genes.tsv", include.dirs = F) %>%
  data.frame(ID = get_nth_part(., "_", 1),
             Sample = gsub("_genes.tsv.gz", "", .), stringsAsFactors = F)
x$Sample <- paste(get_nth_part( x$Sample, "_", 2), 
                  get_nth_part( x$Sample, "_", 3), sep="_")

# > x                                   .         ID            Sample
# 1   GSM4041150_healthy1_cd45+_genes.tsv.gz GSM4041150    healthy1_cd45+
# 2   GSM4041151_healthy1_cd45-A_genes.tsv.gz GSM4041151   healthy1_cd45-A
# 3   GSM4041152_healthy1_cd45-B_genes.tsv.gz GSM4041152   healthy1_cd45-B
# 4   GSM4041153_healthy2_cd45+_genes.tsv.gz GSM4041153    healthy2_cd45+
# 5   GSM4041154_healthy2_cd45-_genes.tsv.gz GSM4041154    healthy2_cd45-
# 6   GSM4041155_healthy3_cd45+_genes.tsv.gz GSM4041155    healthy3_cd45+
  
# In each of these directories, rename the files
lapply(1:nrow(x), FUN = function(i){
  # Create a new directory with sample name - eg. "M4_Vehicle" instead of "23790_M4_Vehicle" 
  new_dir <- create_folder(x$Sample[i])
  # List all the files in the folder
  files <- list.files(pattern = x$ID[i])
  # Specify the order of files
  files <- files[c(grep("barcodes", files), grep("genes", files), grep("matrix", files))]
  # Make the new file names
  new.files <- paste(new_dir, "/",  c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"), sep = "")
  # Rename the files
  file.rename(files, new.files)
})

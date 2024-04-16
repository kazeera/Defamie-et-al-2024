rm(list=ls()[ls()!= "sc"])

#' Export data needed for integration with other datasets

# Load libraries and functions
library(dplyr)
library(plyr)
library(Seurat)

out_dir = "Plots"; paper_author = "Yang"

# Load Seurat data
sc <- readRDS("GSE171904_2021-01-04_seurat.rds")

# Save to file for integration with other datasets
df <- cbind(PDGFRA=FetchData(sc, "Pdgfra"),
            sc@meta.data[,c("sample", "List1_scType", "List2_scType", "List3_scType")])
head(df)
saveRDS(df, "5_Yang_Pdgfr_scType.rds")


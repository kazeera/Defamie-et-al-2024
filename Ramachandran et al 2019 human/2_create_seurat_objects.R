library(kazutils)
library (Seurat)
library(dplyr)
# https://satijalab.org/seurat/v3.1/merge_vignette.html
options(max.print= 45, width = 10)

# Get sample directories
setwd("Data/GSE136103_RAW")
data_dirs <- list.files()
seurat_obs <- list()

# QC in paper's methods section: 
# cells < 300 genes or mitochondrial gene content >30% of the total UMI count. 
# Genes expressed in fewer than three cells in a sample were excluded
for(data_dir in data_dirs){
  # Read data
  sc.data <- Read10X(data.dir = data_dir)  
  # Create object 
  sc <- CreateSeuratObject(counts = sc.data, project = data_dir, min.cells = 3, min.features = 300)
  # Mitochondrial content
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-") ## "mt" instead of "MT" for mouse
  
  # Save RData to sample's directory
  saveRDS(sc, file = sprintf("%s.rds", data_dir))
  seurat_obs <- c(seurat_obs, sc)
}

names(seurat_obs) <- data_dirs
save(seurat_obs, data_dirs, file="../seurat_obs.RData")


# https://bioinformatics.stackexchange.com/questions/9225/how-to-merge-more-than-two-sample-in-seurat

# optional but probably a good idea
# rename cells using object names as prefix
for (i in names(seurat_obs)) {
  seurat_obs[[i]] <- RenameCells(seurat_obs[[i]],
                                         add.cell.id = i)
}
# merge function in Seurat v>3.0
merged_combined <- merge(x = seurat_obs[[1]], y =  seurat_obs[2:length(seurat_obs)])
# Add to metadata
merged_combined$Condition <- get_nth_part(merged_combined$orig.ident, "_", 1)
merged_combined$Condition  %>% unique
merged_combined$Condition <- ifelse(grepl("healthy", merged_combined$Condition), "healthy", "cirrhotic")

tail(merged_combined@meta.data)
# orig.ident nCount_RNA nFeature_RNA percent.mt Condition
# healthy5_cd45+_TTTGTCACATGCCCGA-1 healthy5_cd45+       2351          855   1.999149   healthy
# healthy5_cd45+_TTTGTCAGTCGAGATG-1 healthy5_cd45+       2306          825   5.333912   healthy
# healthy5_cd45+_TTTGTCAGTGATGTCT-1 healthy5_cd45+       2766          999   2.856110   healthy
dim(merged_combined) # rows (genes) x columns (cells)
# [1] 23039 61665

# QC in paper's methods section: 
# cells < 300 genes or mitochondrial gene content >30% of the total UMI count. 
merged_combined <- subset(merged_combined, percent.mt <= 30)
# > dim(sc2) [1] 23039 60925
saveRDS(merged_combined, file="../merged_samp_Seurat.rds")

rownames(merged_combined)[grepl("PDGF", rownames(merged_combined))]

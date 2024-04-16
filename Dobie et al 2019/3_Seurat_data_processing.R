library(Seurat)
library(cowplot)
library(ggplot2)
library(kazutils) # from kazeera's github
library(dplyr)

# Load normal data matrix
sc <- readRDS("../Data/merged_samp_Seurat.rds")
sc <- subset(sc, orig.ident != "CCl4 72h_HSC")
# # Normalize
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# Scale and PCA
sc <- ScaleData(sc, features = rownames(sc))
sc <- RunPCA(sc, verbose = F)

#You can also achieve this using the ElbowPlot function
win.graph()
ElbowPlot(sc) # levels out at dim=20

# Cluster the cells
sc <- FindNeighbors(sc, dims = 1:15, verbose = FALSE)
sc <- FindClusters(sc, resolution = 0.5, verbose = FALSE)
# Run non-linear dimensional reduction
sc <- RunTSNE(sc, dims = 1:15, verbose = FALSE, check_duplicates = FALSE)
sc <- RunUMAP(sc, dims = 1:15, verbose = FALSE)

# Add experiment to metadata and set new identities
sc$Expt <- sc$orig.ident 
sc$Condition <- as.character(sc$Expt) %>% get_nth_part("_", 1)
Idents(sc) <- sc$Condition 
sc$Treatment <- as.character(sc$Condition) %>% get_nth_part(" ", 1)


# Save processed Seurat object
saveRDS(sc, file = "processed_Seurat_ob_noCCl4.rds")

## 3. Process data  -----------------------------------------
library(Seurat)
library(cowplot)
library(ggplot2)
library(kazutils) # from kazeera's github
library(dplyr)

# Load data matrix
# Add experiment to metadata and set new identities
# ## Mouse
# sc <- readRDS("~/liver/sc_mouse.rds")
# sc$species <- "mouse"
# sc$timepoint  <- factor(sc$timepoint, levels=c("E11", "E11.5", "E13", "E14.5", "E16", "E17.5"))
# sigs <- read.xlsx("~/liver/Mouse_cell_signature_curated.xlsx") %>% as.list()

## Human
sc <- readRDS("~/liver/sc_human.rds")
sc$species <- "human"
sc$timepoint  <- factor(sc$timepoint, levels=c("5w","6w","7w","9w","12w","13w","14w","16w","19w"))
sigs <- read.xlsx("~/liver/Human_cell_signature_curated.xlsx") %>% as.list()

# Apply identities, join and filter data
sigs <- lapply(sigs, function(x) x[!is.na(x)])
sc <- JoinLayers(sc)
Idents(sc) <- sc$timepoint 

dim(sc)
# VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "species")
sc <- subset(sc, nFeature_RNA > 200 & percent.mt < 30 )
# VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "species")
dim(sc)

# Normalize, scale, run linear dim red
sc <- sc %>%
  Seurat::NormalizeData(verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = VariableFeatures(object = sc), verbose = FALSE)

# # Run HARMONY
# sc <- sc %>% RunHarmony("timepoint", plot_convergence = TRUE)
# # Access harmony embeddings
# harmony_embeddings <- Embeddings(sc, 'harmony')
# harmony_embeddings[1:5, 1:5]
# # Cluster the cells
# sc <- sc %>% 
#   RunUMAP(dims = 1:12, reduction = "harmony") %>%  # reduction = "harmony"
#   FindNeighbors(dims = 1:12, reduction = "harmony") %>% 
#   FindClusters(resolution = 0.5) %>% 
#   identity()

# Cluster the cells
sc <- sc %>% 
  FindNeighbors(dims = 1:12) %>% 
  FindClusters(resolution = 0.5) 

# Run non-linear dimensional reduction
sc <- RunTSNE(sc, dims = 1:12, verbose = FALSE, check_duplicates = FALSE)
sc <- RunUMAP(sc, dims = 1:12, verbose = FALSE)

# Add gene set scores
sc <- AddModuleScore(sc, sigs)
colnames(sc@meta.data)[grep("Cluster", colnames(sc@meta.data))] <- names(sigs)

# Save processed Seurat object
# saveRDS(sc, "~/liver/sc_mouse_Seurat.rds")
# sc_human=sc
saveRDS(sc, "~/liver/sc_human_Seurat.rds")

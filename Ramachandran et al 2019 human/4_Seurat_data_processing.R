library(Seurat)
library(Nebulosa)
library(cowplot)
library(dplyr)
# https://github.com/satijalab/seurat/issues/4436
# remove.packages("Matrix")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
library(Matrix)

# Read in 
sc <- readRDS("../Data/merged_samp_Seurat.rds")
dim(sc) # [1] 23039 60925
sc$Condition %>% unique # "cirrhotic" "healthy" 

# Normalize
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# Scale and PCA
sc <- ScaleData(sc, features = rownames(sc))
sc <- RunPCA(sc, verbose = F)

# #Next, we need to determine the dimensionality of the dataset, this is time consuming for big data
# sc <- JackStraw(sc, num.replicate = 100) # takes like 5 mins
# sc <- ScoreJackStraw(sc, dims = 1:20)
# 
#  #To visualize the p-value of each PC dimension use
# JackStrawPlot(sc, dims = 1:20)
# #^this shows that p-values drop off a lot in the PC10-12 range, thus its our choice how many to include in downstream anaylses

#You can also achieve this using the ElbowPlot function
win.graph()
ElbowPlot(sc) # levels out at dim=20

# Cluster the cells
sc <- FindNeighbors(sc, dims = 1:20, verbose = FALSE)
sc <- FindClusters(sc, resolution = 0.5, verbose = FALSE)
# Run non-linear dimensional reduction
sc <- RunTSNE(sc, dims = 1:20, verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:20, verbose = FALSE)

# Save processed Seurat object
saveRDS(sc, file = "processed_Seurat_ob.rds")

# Create UMAPs
out_dir <- kazutils::create_folder("Plots")

jpeg(sprintf("%s/0_total_umap.jpeg", out_dir), height = 480, width = 900)
plot1 <- DimPlot(sc, reduction = "tsne", pt.size = 0.5)+ ggtitle(sprintf("%s healthy total liver cells", ncol(sc)))
plot2 <- DimPlot(sc, group.by = "seurat_clusters", pt.size = 0.5)
plot_grid(plot1,plot2)
dev.off()




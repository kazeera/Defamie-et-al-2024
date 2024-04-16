library(Seurat)
library(RColorBrewer)
library(ggplot2)

# Load Seurat data
sc <- readRDS("GSE171904_2021-01-04_seurat.rds")

# Get celltype metadata
df <- Seurat::FetchData(sc, "celltype")
colnames(df) <- "Yang_celltype"

# Write to file
saveRDS(df, file="6_Yang_paper_celltype_metadata.rds")

# Plot feature plot
pdf("Plots/6_Yang_FeaturePlot_Sox9_Pdgfra.pdf", width = 10, height= 5)
FeaturePlot(sc, features = c("Pdgfra", "Sox9"))
dev.off()

# Plot feature plot
pdf("Plots/6_Yang_DimPlot_celltype.pdf", width = 7.5, height= 6.5)
DimPlot(sc, group.by = "celltype", label = T) + ggtitle("Paper's celltypes")
dev.off()

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
# Plot feature plot
pdf("Plots/6_Yang_VlnPlot_Pdgfra_celltype.pdf", width = 10, height= 5)
VlnPlot(sc,features = "Pdgfra", group.by = "celltype", pt.size = 0)
dev.off()

# Plot feature plot
pdf("Plots/6_Yang_VlnPlot_Pdgfra_celltype2.pdf", width = 10, height= 5)
VlnPlot(sc,features = "Pdgfra", group.by = "celltype", pt.size = 0) +
  stat_summary(fun.y = median.stat, geom='crossbar', size = 0.3, colour = "black") 
dev.off()


# Colors
cell_colors <- c(oil="#00474D", ccl4="#E51E1E", bdl="#faa7a7")

# Plot sample plots
pdf("Plots/6_Yang_DimPlot_condition_split.pdf", width = 15, height= 5)
DimPlot(sc, split.by = "sample", group.by = "sample", cols = cell_colors)
dev.off()


pdf("Plots/6_Yang_DimPlot_condition.pdf", width = 7.5, height= 6.5)
DimPlot(sc, group.by = "sample", cols = cell_colors)
dev.off()

# Other
DimPlot(sc, group.by = "seurat_cluster")
VlnPlot(sc, group.by = "List3_scType", split.by = "sample", "Pdgfra")
DimPlot(sc, group.by = "celltype", label = T)
FeaturePlot(sc, "Sox9")

rm(list=ls()[ls()!= "sc"])

#' Goal - Reassign celltypes using simpler/better R package called "scType" (compared to "cellassign")
# Note: gs = gene set

# Load libraries and functions
library(kazutils) # devtools::install_github("kazeera/kazutils")
library(dplyr)
library(plyr)
library(Seurat)
library(ggplot2)

# Load Seurat data
sc <- readRDS("processed_Seurat_ob.rds")

# New annotation combining treated groups
group_colors <- c(healthy="#00474D", cirrhotic="#E51E1E")
# scales::show_col(group_colors)


# Read in genes of interest
markers <- readRDS("../../GeneLists/markers_of_interest_Virginie_updated20221107.rds") %>% 
  toupper

## Visualize
font_size <- 10
line_size <- 1
group_metacolumn <- "Condition"
gene_boxplot <- "DLK1"
paper_author <- "Ramachandran-human"
gs_name <- "List2_scType"

# Create output folder
out_dir <- kazutils::create_folder("Plots")

# Condition UMAP
pdf(sprintf("Plots/4_DimPlot_%s_Condition.pdf", paper_author), width = 7.5, height= 6.5)
DimPlot(sc, label = FALSE, group.by = "Condition", cols = group_colors)
dev.off()

# Condition UMAP split
pdf(sprintf("Plots/4_DimPlot_%s_Condition_split.pdf", paper_author), width = 10, height= 5)
DimPlot(sc, group.by = "Condition", split.by = "Condition", cols = group_colors)
dev.off()

# Order of celltypes - remove Bipotent.progenitor
celltype_order <- c("Hepatocyte", "Cholangiocyte", "Endothelia", "Mesenchyme", "Fibroblast", "HSC", "VSMC", "Macrophage", "B.cell", "T.cell", "ILC", "Denditric.cell", "Plasma.cell")

# Get name in Seurat meta data
# Subset out Unknown
sc2 <- subset(sc, List2_scType != "Unknown" & List2_scType != "Bipotent.progenitor")

# Factor based on list
celltype_order2 <- celltype_order[celltype_order %in% sc2$List2_scType]
sc2$List2_scType <- factor(x = sc2$List2_scType, levels = celltype_order2)

# Stacked violin plot
pdf(sprintf("%s/4_%s_%s_vln.pdf", out_dir, paper_author, "List2_scType"), height = 9, width = 7)
VlnPlot(sc2, features = rev(markers), split.by = group_metacolumn, cols = group_colors, group.by = "List2_scType", pt.size = 0.1, stack = T,flip = T)+ ggtitle("List2_scType")
dev.off()

# Plot feature plot
pdf("Plots/4_Ramachandran-human_DimPlot_List2_scType.pdf", width = 7.5, height= 6.5)
DimPlot(sc, group.by = "List2_scType", label = T) 
dev.off()

# # Subset to fibroblasts only
# sc_fibro <- subset(sc2, List2_scType == "Fibroblast")
# 
# # Make data values table
# df3 <- data.frame(group = sc_fibro[[group_metacolumn]],
#                   celltype = sc_fibro$List2_scType,
#                   gene = Seurat::FetchData(object = sc_fibro, slot = "data", vars = c("ADAM10", "ADAM17", "DLK1")))
# 
# # Make correlation scatter in high res
# pdf("Plots/4_Ramachandran-human_Scatter_DLK1_ADAM17.pdf", width = 7.5, height= 6.5)
# FeatureScatter(sc_fibro, "ADAM17", "DLK1", pt.size = 2, group.by = gs_name, cols = "#CB1B4FFF")+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s ADAM17+DLK1+, %s ADAM17+ only, %s DLK1+ only", sum(df3$gene.ADAM17 != 0 & df3$gene.DLK1 != 0), sum(df3$gene.ADAM17 != 0 & df3$gene.DLK1 == 0), sum(df3$gene.ADAM17 == 0 & df3$gene.DLK1 != 0))))
# dev.off()
# 
# pdf("Plots/4_Ramachandran-human_Scatter_DLK1_ADAM10.pdf", width = 7.5, height= 6.5)
# FeatureScatter(sc_fibro, "ADAM10", "DLK1", pt.size = 2, group.by = gs_name, cols = "#CB1B4FFF")+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s ADAM10+DLK1+, %s ADAM10+ only, %s DLK1+ only", sum(df3$gene.ADAM10 != 0 & df3$gene.DLK1 != 0), sum(df3$gene.ADAM10 != 0 & df3$gene.DLK1 == 0), sum(df3$gene.ADAM10 == 0 & df3$gene.DLK1 != 0))))
# dev.off()



# Visualize/analyze data 
library(Seurat)
library(cowplot)
library(ggplot2)
library(kazutils) # from kazeera's github
library(dplyr)
library(RColorBrewer)
library(viridis)
library(svglite)

# Create color palette
color_pal <- c("#AA336A", brewer.pal(12, "Paired"), "#00008B", "#FF5733", "#F08000", "#06768d", "#023020", "#80461B", brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))
color_pal <- color_pal[!color_pal %in% c("#1F78B4", "#E78AC3", "#FC8D62", "#FF7F00", "#E41A1C", "#A65628", "#B2DF8A", "")]
scales::show_col(color_pal)
color_pal_current <- color_pal[c(1:4, 9, 6:8, 5, 10:12)] # color pal for Virginie - switch color between CD68 and Timp1

## Mouse ----------------------
sc_mouse <- readRDS("~/liver/sc_mouse_Seurat.rds")
sc <- sc_mouse

# Filter out doublets
sc <- subset(sc, cluster != "doublet")

# Combine celltype groups
celltype_map <- c("Hepatoblast/Hepatocyte"="Hepatoblast/Hepatocyte",
                  "Cholangiocyte"="Cholangiocyte", 
                  "Erythroid progenitor"="ErythroidLineage", "Erythroblast"="ErythroidLineage",
                  "Primitive erythrocyte"="ErythroidLineage", "Early erythrocyte"="ErythroidLineage",
                  "Myeloid/lymphoid/megakaryoid cell"="HematopoieticLineage", 
                  "Hepatic stellate cell"="HepaticStellate",
                  "Hematopoietic stem and progenitor cell"="HematopoieticLineage",
                  "Kupffer cell"="Kupffer",
                  "Liver endothelial cell"="Endothelial",
                  "Mesothelial cell"="Mesothelial",
                  "Septum transversumal cell"="HepaticStellate")

# Map and reorder
sc$cluster_combined <- plyr::mapvalues(sc$cluster, names(celltype_map), unname(celltype_map))
celltype_order <- c("Hepatoblast/Hepatocyte", "Cholangiocyte", "ErythroidLineage", 
                    "HematopoieticLineage", "Kupffer", "HepaticStellate", "Mesothelial", "Endothelial")
sc$cluster_combined <- factor(sc$cluster_combined, levels=celltype_order)

# Add to original metadata
df <- FetchData(sc, "cluster_combined")
sc_mouse <- AddMetaData(sc_mouse, df)
# saveRDS(sc_mouse, "~/liver/sc_mouse_Seurat.rds")

# Final cell type marker gene list
genes_final <- c("Hnf4a", "Sox9", "Gypa", "Fcnb", "Cd68", "Pdgfrb", "Pdpn", "Flt4", "Timp1", "Timp2", "Timp3", "Timp4")

# Final violin plots
p1 <- VlnPlot(sc, genes_final, group.by = "cluster_combined", cols = color_pal_current, 
              stack=T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90))
p2 <- VlnPlot(sc, genes_final, group.by = "cluster_combined", cols = color_pal_current, 
              stack=T, flip = T, split.by = "timepoint") + NoLegend() + theme(axis.text.x = element_text(angle = 90), legend.position = "top")
# Final Tsne plots
p3 <- DimPlot(sc, cols=rev(color_pal)[-6], group.by = "cluster_combined", reduction = "tsne", raster = F)
p4 <- DimPlot(sc, cols=color_pal, group.by = "timepoint", reduction = "tsne", raster = F)

pdf("~/liver/results/241120_mouse_violin.pdf", height = 4, width = 6)
print(p1)
print(p2)
dev.off()

svglite("~/liver/results/241120_mouse_tsne1.svg", width = 6, height = 4)
print(p3)
dev.off()

svglite("~/liver/results/241120_mouse_tsne2.svg", width = 4.5, height = 4)
print(p4)
dev.off()

## Human ----------------------
# Read data
sc_human <- readRDS("~/liver/sc_human_Seurat.rds")

sc <- sc_human
# Filter out doublets
sc <- subset(sc, cluster != "doublet")

# Combine celltype groups
celltype_map <- c("Hepatoblast/Hepatocyte"="Hepatoblast/Hepatocyte",
                  "Cholangiocyte"="Cholangiocyte", 
                  "Erythroid progenitor"="ErythroidLineage", "Erythroblast"="ErythroidLineage",
                  "Primitive erythrocyte"="ErythroidLineage", "Early erythrocyte"="ErythroidLineage",
                  "Myeloid/lymphoid/megakaryoid cell"="HematopoieticLineage", 
                  "Hepatic stellate cell"="HepaticStellate",
                  "Hematopoietic stem and progenitor cell"="HematopoieticLineage",
                  "Kupffer cell"="Kupffer",
                  "Liver endothelial cell"="Endothelial",
                  "Mesothelial cell"="Mesothelial",
                  "Septum transversumal cell"="HepaticStellate")

# Map and reorder
sc$cluster_combined <- plyr::mapvalues(sc$cluster, names(celltype_map), unname(celltype_map))
celltype_order <- c("Hepatoblast/Hepatocyte", "Cholangiocyte", "ErythroidLineage", 
                    "HematopoieticLineage", "Kupffer", "HepaticStellate", "Mesothelial", "Endothelial")
sc$cluster_combined <- factor(sc$cluster_combined, levels=celltype_order)

# Add to original metadata
df <- FetchData(sc, "cluster_combined")
sc_human <- AddMetaData(sc_human, df)
# saveRDS(sc_human, "~/liver/sc_human_Seurat.rds")

# Final cell type marker gene list
genes_final <- c("Hnf4a", "Sox9", "Gype", "Plek", "Cd68", "Pdgfrb", "Pdpn", "Flt1", "Timp1", "Timp2", "Timp3", "Timp4") %>% toupper

# Final Violin plot showing celltype markers
p1 <- VlnPlot(sc, genes_final, group.by = "cluster_combined", cols = color_pal_current, 
              stack=T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90))

p2 <- VlnPlot(sc, genes_final, group.by = "cluster_combined", cols = color_pal_current, 
              stack=T, flip = T, split.by = "timepoint") + theme(axis.text.x = element_text(angle = 45), legend.position = "top")

# Final Tsne plots
p3 <- DimPlot(sc, cols=rev(color_pal)[-6], group.by = "cluster_combined", reduction = "tsne", raster = F)
p4 <- DimPlot(sc, cols=color_pal, group.by = "timepoint", reduction = "tsne", raster = F)

pdf("~/liver/results/241120_human_vlnplots.pdf", height = 4, width = 6)
# pdf("~/241119_human_plots.pdf", height = 5, width = 6)
print(p1)
print(p2)
dev.off()

svglite("~/liver/results/241120_human_dimplot1.svg", width = 6, height = 4)
print(p3)
dev.off()

svglite("~/liver/results/241120_human_dimplot2.svg", width = 4.5, height = 4)
print(p4)
dev.off()

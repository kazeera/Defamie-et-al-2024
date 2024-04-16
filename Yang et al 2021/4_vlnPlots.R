# Load libraries and functions
library(kazutils) # devtools::install_github("kazeera/kazutils")
library(dplyr)
library(plyr)
library(Seurat)
library(ggplot2)
library(cowplot)

# Load Seurat data
sc <- readRDS("GSE171904_2021-01-04_seurat.rds")

# Create output folder
out_dir <- kazutils::create_folder("Plots")

# Define group colors
group_colors <- c(oil="#00474D",ccl4="#E51E1E",bdl="#f78686")
sc$sample2 <- factor(sc$sample, levels = c("oil", "ccl4", "bdl"))

# Stacked violin plot
jpeg(sprintf("%s/4_%s_umap_treatment.jpeg", out_dir, paper_author), height =  900, width = 900)
# pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
p1 <- DimPlot(sc,group.by = "sample2", cols = group_colors, pt.size = 0.02)
print(cowplot::plot_grid(p1)); 
dev.off()

# Stacked violin plot
jpeg(sprintf("%s/4_%s_umap_Dlk.jpeg", out_dir, paper_author), height =  900, width = 900)
# pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
p1 <- FeaturePlot(sc, "Dlk1", pt.size = 0.5) # see how to show expressing cells on top
print(cowplot::plot_grid(p1)); 
dev.off()

# Specify parameters
font_size <- 12
line_size <- 1
group_metacolumn <- "sample2"
gene_boxplot <- "Dlk1"
paper_author <- "Yang"

# Read in genes of interest
markers <- readRDS("../../GeneLists/markers_of_interest_Virginie_updated20221107.rds")
markers <- markers[markers!="Dll3"] # Warning: All cells have the same value of Dll3
# markers <- c("Cnn1", "Cd34", "Lrat", markers)

# Order of celltypes - remove Bipotent.progenitor
celltype_order <- c("Hepatocyte", "Cholangiocyte", "Bipotent.progenitor", "Endothelia", "Mesenchyme", "Fibroblast", "HSC", "VSMC", "Macrophage", "B.cell", "T.cell", "ILC", "Denditric.cell", "Plasma.cell")

# Iteratively, create UMAPs and violin plots
for(gs_name in c("List2longnoBip", "List2long")){
  # Get name in Seurat meta data
  gs_name <- paste0(gs_name, "_scType")
  
  # Subset out Unknown
  sc$temp <- sc[[gs_name]]
  sc2 <- subset(sc, temp != "Unknown")
  
  # Factor based on list
  celltype_order2 <- celltype_order[celltype_order %in% sc2$temp]
  sc2$temp <- factor(x = sc2$temp, levels = celltype_order2)
  
  # Stacked violin plot
  jpeg(sprintf("%s/4_%s_%s_vln.jpeg", out_dir, paper_author, gs_name), height = 1000, width = 700)
  # pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
  plot2 <- VlnPlot(sc2, features = rev(markers), split.by = group_metacolumn, cols = group_colors, group.by = "temp", pt.size = 0.1, stack = T,flip = T)+ ggtitle(gs_name)
  print(plot2); dev.off()
  
  # # UMAP
  jpeg(sprintf("%s/4_%s_%s_umap.jpeg", out_dir, paper_author, gs_name), height = 500, width = 500)
  plot3 <- DimPlot(sc2, label = TRUE, group.by = gs_name)
  print(plot3); dev.off()
  
  # Make data values table
  df3 <- data.frame(group = sc2[[group_metacolumn]],
                    celltype = sc2$temp, 
                    gene = Seurat::FetchData(object = sc2, slot = "data", vars = c(gene_boxplot)))
  colnames(df3) <- c("group","celltype", "gene")
  
  # Remove 0 values
  df3 <- df3[df3$gene != 0, ]
  
  # Count Dlk cells
  counts_Dlk <- df3 %>%
    group_by(group, celltype)  %>%
    tally() %>% 
    data.frame()
  
  # Count total cells
  counts_total <- sc2@meta.data[,c(gs_name, group_metacolumn)] %>%
    table %>%
    reshape2::melt()
  
  # Add MergeID columns
  counts_Dlk$MergeID <- paste(counts_Dlk$group, counts_Dlk$celltype, sep="_")
  counts_total$MergeID <- paste(counts_total[,group_metacolumn], counts_total[,gs_name], sep="_")
  
  # Merge 
  counts_merged <- merge(counts_Dlk, counts_total, by="MergeID", all=T)
  counts_merged <- counts_merged[!is.na(counts_merged$n),]
  
  # Find proportion
  counts_label <- paste0(counts_merged$MergeID, ": ", counts_merged$n, "/", counts_merged$value) %>%
    paste(collapse = "; ")
  
  # Make boxplot
  plot4 <- ggplot(df3, aes(x = celltype, y = gene, color=group)) +
    geom_point(size=2,alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height =  0)) +
    scale_color_manual(values = group_colors) +
    labs(
      title = gene_boxplot,
      subtitle = sprintf("%s et al dataset; Crossbar = median", paper_author),
      caption = stringr::str_wrap(paste("Dlk/Total cells:", counts_label)),
      y = "Expression",
      x = gs_name
    ) +
    # scale_y_continuous(trans = "log10") +
    theme(
      panel.background = element_blank(), # remove background color and lines
      axis.line = element_line(colour = "black", size = line_size), # increase the axis-line thickness and change the color to blac
      # Ticks
      axis.ticks = element_line(colour = "black", size = line_size), # increase the tick thickness)
      axis.ticks.length = unit(.25, "cm"),
      # Axes labels
      axis.text = element_text(colour = "black", size = font_size),
      axis.text.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0), hjust = 1, vjust = 1, angle = 45), # increase space between x axis title and labels
      axis.text.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
      # axes tick labels
      axis.title = element_text(colour = "black", size = font_size, face = "bold"), # axes title labels
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), # increase space between x axis title and labels
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      # legend
      legend.position = "top",
      legend.text = element_text(colour = "black", size = font_size),
      legend.title = element_blank()
    )+ stat_summary(fun = median, fun.min = median, fun.max = median,position = position_dodge(width = 0.75), geom = "crossbar", width = 0.5)
  ggsave(plot4, filename = sprintf("%s/4_%s_%s_%s.png", out_dir, paper_author, gs_name, gene_boxplot))
  
  # Subset to fibroblasts only
  sc_fibro <- subset(sc2, temp == "Fibroblast")
  # Make data values table
  df3 <- data.frame(group = sc_fibro[[group_metacolumn]],
                    celltype = sc_fibro$temp,
                    gene = Seurat::FetchData(object = sc_fibro, slot = "data", vars = c("Adam10", "Adam17", "Dlk1")))
  p1 <- FeatureScatter(sc_fibro, "Adam17", "Dlk1", pt.size = 2, group.by = gs_name)+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s Adam17+Dlk1+, %s Adam17+ only, %s Dlk1+ only", sum(df3$gene.Adam17 != 0 & df3$gene.Dlk1 != 0), sum(df3$gene.Adam17 != 0 & df3$gene.Dlk1 == 0), sum(df3$gene.Adam17 == 0 & df3$gene.Dlk1 != 0))))
  p2 <- FeatureScatter(sc_fibro, "Adam10", "Dlk1", pt.size = 2, group.by = gs_name)+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s Adam10+Dlk1+, %s Adam10+ only, %s Dlk1+ only", sum(df3$gene.Adam10 != 0 & df3$gene.Dlk1 != 0), sum(df3$gene.Adam10 != 0 & df3$gene.Dlk1 == 0), sum(df3$gene.Adam10 == 0 & df3$gene.Dlk1 != 0))))
  
  # correlation scatter
  jpeg(sprintf("%s/4_%s_%s_Adam.vs.Dlk.jpeg", out_dir, paper_author, gs_name), height =  480, width = 900)
  # pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
  print(cowplot::plot_grid(p1, p2)); 
  dev.off()
  
}


# Mar 27, 2023
# Print pdf ====================================================
# Subset out Unknown
sc2 <- subset(sc, List2longnoBip_scType != "Unknown")

# Order of celltypes - remove Bipotent.progenitor
celltype_order <- c("Hepatocyte", "Cholangiocyte", "Endothelia", "Mesenchyme", "Fibroblast", "HSC", "VSMC", "Macrophage", "B.cell", "T.cell", "ILC", "Denditric.cell", "Plasma.cell")

# Factor based on list
celltype_order2 <- celltype_order[celltype_order %in% sc2$List2longnoBip_scType]
sc2$List2longnoBip_scType <- factor(x = sc2$List2longnoBip_scType, levels = celltype_order2)

# Stacked violin plot
pdf(sprintf("%s/4_%s_%s_vln.pdf", out_dir, paper_author, "List2longnoBip_scType"), height = 9, width = 7)
# pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
plot2 <- VlnPlot(sc2, features = rev(markers), split.by = group_metacolumn, cols = group_colors, group.by = "List2longnoBip_scType", pt.size = 0.1, stack = T,flip = T)+ ggtitle("List2longnoBip_scType")
print(plot2); dev.off()

# Subset to fibroblasts only
sc_fibro <- subset(sc2, List2longnoBip_scType == "Fibroblast")
gs_name <- "List2longnoBip_scType"

# Make data values table
df3 <- data.frame(group = sc_fibro[[group_metacolumn]],
                  celltype = sc_fibro$List2longnoBip_scType,
                  gene = Seurat::FetchData(object = sc_fibro, slot = "data", vars = c("Adam10", "Adam17", "Dlk1")))

# Make correlation scatter in high res
pdf("Plots/4_Yang_Scatter_Dlk1_Adam17.pdf", width = 7.5, height= 6.5)
FeatureScatter(sc_fibro, "Adam17", "Dlk1", pt.size = 2, group.by = gs_name, cols = "#CB1B4FFF")+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s Adam17+Dlk1+, %s Adam17+ only, %s Dlk1+ only", sum(df3$gene.Adam17 != 0 & df3$gene.Dlk1 != 0), sum(df3$gene.Adam17 != 0 & df3$gene.Dlk1 == 0), sum(df3$gene.Adam17 == 0 & df3$gene.Dlk1 != 0))))
dev.off()

pdf("Plots/4_Yang_Scatter_Dlk1_Adam10.pdf", width = 7.5, height= 6.5)
FeatureScatter(sc_fibro, "Adam10", "Dlk1", pt.size = 2, group.by = gs_name, cols = "#CB1B4FFF")+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s Adam10+Dlk1+, %s Adam10+ only, %s Dlk1+ only", sum(df3$gene.Adam10 != 0 & df3$gene.Dlk1 != 0), sum(df3$gene.Adam10 != 0 & df3$gene.Dlk1 == 0), sum(df3$gene.Adam10 == 0 & df3$gene.Dlk1 != 0))))
dev.off()

# Plot feature plot
pdf("Plots/4_Yang_DimPlot_List2longnoBip_scType.pdf", width = 7.5, height= 6.5)
DimPlot(sc, group.by = "List2longnoBip_scType", label = T) 
dev.off()

pdf("Plots/4_Yang_VlnPlot_Dlk1_List2longnoBip_scType.pdf", width = 10, height= 5)
VlnPlot(sc,features = "Dlk1", group.by = "List2longnoBip_scType", pt.size = 1)
dev.off()

pdf("Plots/4_Yang_VlnPlot_Adam17_List2longnoBip_scType.pdf", width = 10, height= 5)
VlnPlot(sc,features = "Adam17", group.by = "List2longnoBip_scType", pt.size = 1)
dev.off()

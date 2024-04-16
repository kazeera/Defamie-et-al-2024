rm(list=ls()[ls()!= "sc"])
#' Goal - Reassign celltypes using simpler/better R package called "scType" (compared to "cellassign")
# Note: gs = gene set

# Load libraries and functions
library(kazutils) # devtools::install_github("kazeera/kazutils")
library(dplyr)
library(plyr)
library(Seurat)

# Load Seurat data
sc <- readRDS("GSE171904_2021-01-04_seurat.rds")

# Read in genes of interest
markers <- readRDS("../../GeneLists/markers_of_interest_Virginie_updated20221107.rds")

# New annotation combining treated groups
sc$sample2 <- plyr::mapvalues(sc$sample, c("oil", "ccl4", "bdl"), c("untreated", "treated", "treated"))
sc$sample2 <- factor(sc$sample2, levels = c( "treated", "untreated"))
group_colors <- c( treated="#E51E1E",untreated="#00474D")
# scales::show_col(group_colors)
# group_colors <- c(untreated="#078992", treated="#f8766d")

## Visualize
library(ggplot2)

font_size <- 12
line_size <- 1
group_metacolumn <- "sample2"
gene_boxplot <- "Dlk1"
paper_author <- "Yang"

# Create output folder
out_dir <- kazutils::create_folder("Plots")

# Order of celltypes - remove Bipotent.progenitor
celltype_order <- c("Hepatocyte", "Cholangiocyte", "Bipotent.progenitor", "Endothelia", "Mesenchyme", "Fibroblast", "HSC", "VSMC", "Macrophage", "B.cell", "T.cell", "ILC", "Denditric.cell", "Plasma.cell")

# Change to Denditric.cell
sc$List2_scType[ sc$List2_scType == "Denditric.cells"] <- "Denditric.cell"
sc$List3_scType[ sc$List3_scType == "Denditric.cells"] <- "Denditric.cell"

# Iteratively, create UMAPs and violin plots
for(gs_name in c("List2longnoBip", "List2longnoBip")){
  # Get name in Seurat meta data
  gs_name <- paste0(gs_name, "_scType")
  
  # Subset out Unknown
  sc$temp <- sc[[gs_name]]
  sc2 <- subset(sc, temp != "Unknown")
  
  # Factor based on list
  celltype_order2 <- celltype_order[celltype_order %in% sc2$temp]
  sc2$temp <- factor(x = sc2$temp, levels = celltype_order2)
  
  # Stacked violin plot
  jpeg(sprintf("%s/2_%s_%s_vln.jpeg", out_dir, paper_author, gs_name), height = 1000, width = 700)
  # pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
  plot2 <- VlnPlot(sc2, features = rev(markers), split.by = group_metacolumn, cols = group_colors, group.by = "temp", pt.size = 0.1, stack = T,flip = T)  # + ggforce::geom_sina(jitt)
  print(plot2); dev.off()
  
  # UMAP
  jpeg(sprintf("%s/2_%s_%s_umap.jpeg", out_dir, paper_author, gs_name), height = 500, width = 500)
  plot3 <- DimPlot(sc2, label = TRUE, group.by = gs_name)
  print(plot3); dev.off()
  
  # Make data values table
  df3 <- data.frame(group = sc2[[group_metacolumn]],
                    celltype = sc2$temp, 
                    gene = Seurat::FetchData(object = sc2, slot = "data", vars = c(gene_boxplot)))
  colnames(df3) <- c("group","celltype", "gene")

  # Make boxplot
  plot4 <- ggplot(df3, aes(x = celltype, y = gene, fill = group)) +
    geom_boxplot(aes(fill = group)) +
    geom_point(alpha = 0.5, position = position_jitter(w = 0.1, h = 0)) + #position = position_jitterdodge()
    scale_fill_manual(values = group_colors) +
    labs(
      title = gene_boxplot,
      # subtitle = gs_name,
      y = "log10 Expression",
      x = gs_name
    ) +
    scale_y_continuous(trans = "log10") +
    theme(
      panel.background = element_blank(), # remove background color and lines
      plot.title = element_text(colour = "black", size = font_size),
      # plot.subtitle = element_text(colour = "black", size = font_size / 2),
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
      legend.title = element_text(colour = "black", size = font_size, face = "bold")
    )
  ggsave(plot4, filename = sprintf("%s/2_%s_%s_%s.png", out_dir, paper_author, gs_name, gene_boxplot), width = 4, height = 5)
}

saveRDS(sc, file = "GSE171904_2021-01-04_seurat.rds")


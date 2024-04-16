rm(list=ls()[ls()!= "sc"])

library(viridis)

#' Goal - Reassign celltypes using simpler/better R package called "scType" (compared to "cellassign")
# Note: gs = gene set

# Load libraries and functions
library(kazutils) # devtools::install_github("kazeera/kazutils")
library(dplyr)
library(plyr)
library(Seurat)
library(ggplot2)
library(reshape2)

# Load Seurat data
sc <- readRDS("processed_Seurat_ob_noCCl4.rds")

# Read in genes of interest
markers <- readRDS("../../GeneLists/markers_of_interest_Virginie_updated20221107.rds")
markers <- markers[markers!="Dll3"] # Warning: All cells have the same value of Dll3
markers <- c("Cnn1", "Cd34", "Lrat", markers)

# Colours
group_colors <- c(Healthy="#00474D", Treated="#E51E1E")
sc$Treatment <- plyr::mapvalues(sc$Treatment, c("CCl4", "Healthy"), c("Treated", "Healthy"))

# Cell colors
scales::show_col(c(rocket(5), treated="#E51E1E",untreated="#00474D"))
scales::show_col(rocket(5)) # picked 3 colors in the middle
cell_colors <- c(Fibroblast="#CB1B4FFF", HSC="#F4875EFF", VSMC="#611F53FF")

# Plot
jpeg(sprintf("%s/3_Treatment_umap.jpeg", out_dir), height = 480, width = 480)
plot1 <- DimPlot(sc, group.by = "Treatment", pt.size = 0.5, cols = group_colors)+ ggtitle(sprintf("%s total; %s healthy; %s treated", ncol(sc), sum(sc$Treatment == "Healthy"), sum(sc$Treatment == "Treated")))
print(cowplot::plot_grid(plot1))
dev.off()

## Visualize

font_size <- 12
line_size <- 1
group_metacolumn <- "Treatment"
gene_boxplot <- "Dlk1"
paper_author <- "Dobie"

# Create output folder
out_dir <- kazutils::create_folder("Plots")

# Order of celltypes - remove Bipotent.progenitor
celltype_order <- c("Hepatocyte", "Cholangiocyte", "Endothelia", "Mesenchyme", "VSMC", "Fibroblast", "HSC", "activatedHSC", "Macrophage", "B.cell", "T.cell", "ILC", "Denditric.cell", "Plasma.cell")

# Iteratively, create UMAPs and violin plots
for(gs_name in c("List2mesenchyme")){
  # Get name in Seurat meta data
  gs_name <- paste0(gs_name, "_scType")
  
  # Subset out Unknown
  sc$temp <- sc[[gs_name]]
  sc2 <- subset(sc, temp != "Unknown" & temp != "Bipotent.progenitor")
  
  # Factor based on list
  celltype_order2 <- celltype_order[celltype_order %in% sc2$temp]
  sc2$temp <- factor(x = sc2$temp, levels = celltype_order2)
  
  # Stacked violin plot
  jpeg(sprintf("%s/3_%s_%s_vln.jpeg", out_dir, paper_author, gs_name), height = 1000, width = 700)
  # pdf(sprintf("%s/1_%s2.pdf", out_dir, gs_name), height=6, width = 4)
  plot2 <-  VlnPlot(sc2, features = rev(markers), split.by = group_metacolumn, cols = group_colors, group.by = "temp", pt.size = 0.1, stack = T,flip = T)  + ggtitle(gs_name)# + ggforce::geom_sina(jitt)
  print(plot2); dev.off()
  
  # UMAP
  jpeg(sprintf("%s/3_%s_%s_umap.jpeg", out_dir, paper_author, gs_name), height = 600, width = 600)
  plot3 <- DimPlot(sc2, label = FALSE, group.by = gs_name, cols = cell_colors)
  print(plot3); dev.off()
  
  # Dlk1 not found
  # Make data values table
  df3 <- data.frame(group = sc2[[group_metacolumn]],
                    celltype = sc2$temp,
                    gene = Seurat::FetchData(object = sc2, slot = "data", vars = c(gene_boxplot)))
  # Rename columns
  colnames(df3) <- c("group","celltype", "gene")
  
  # Remove 0 values
  df3 <- df3[df3$gene != 0, ]
  
  # Count Dlk cells
  counts_Dlk <- df3 %>%
    group_by(group, celltype)  %>%
    tally() %>% 
    data.frame()
  
  # Count total cells
  counts_total <- sc2@meta.data[,c(gs_name, "Treatment")] %>%
    table %>%
    reshape2::melt()
  
  # Add MergeID columns
  counts_Dlk$MergeID <- paste(counts_Dlk$group, counts_Dlk$celltype, sep="_")
  counts_total$MergeID <- paste(counts_total$Treatment, counts_total$List2mesenchyme_scType, sep="_")
  
  # Merge 
  counts_merged <- merge(counts_Dlk, counts_total, by="MergeID", all=T)
  
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
  ggsave(plot4, filename = sprintf("%s/3_%s_%s_%s.png", out_dir, paper_author, gs_name, gene_boxplot))
}




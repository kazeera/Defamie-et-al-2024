library("Nebulosa")
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# Load Seurat data
sc <- readRDS("processed_Seurat_ob_noCCl4.rds")


# Cell colors
scales::show_col(c(rocket(5), treated="#E51E1E",untreated="#00474D"))
sc$Treatment <- plyr::mapvalues(sc$Treatment, c("CCl4", "Healthy"), c("Treated", "Healthy"))

# Colours
group_colors <- c(Healthy="#00474D", Treated="#E51E1E")
cell_colors <- c(Fibroblast="#CB1B4FFF", HSC="#F4875EFF", VSMC="#611F53FF")

# Plot
pdf("Plots/5_Dobie_DimPlot_Treatment_umap.pdf", width = 7.5, height= 6.5)
DimPlot(sc, group.by = "Treatment", cols=group_colors)+ggtitle(sprintf("%s total; %s healthy; %s treated", ncol(sc), sum(sc$Treatment == "Healthy"), sum(sc$Treatment == "Treated")))
dev.off()


# Plot
pdf("Plots/5_Dobie_DimPlot_scType_umap.pdf", width = 7.5, height= 6.5)
DimPlot(sc, group.by = "List2mesenchyme_scType", cols=cell_colors)
dev.off()
# BiocManager::install("Nebulosa")

# Plot
pdf("Plots/5_Dobie_FeaturePlot_Dlk1_Adams.pdf", width = 15, height= 5)
FeaturePlot(sc, c("Dlk1", "Adam17", "Adam10"), ncol = 3)
dev.off()

# Density plot
pdf("Plots/5_Dobie_DensityPlot_Dlk1_Adams.pdf", width = 15, height= 5)
plot_density(sc, c("Dlk1", "Adam17","Adam10"))
dev.off()

# Make scatter plots
gs_name <- "List2mesenchyme_scType"
group_metacolumn <- "Treatment"

# Subset to fibroblast
sc_fibro <- subset(sc, List2mesenchyme_scType == "Fibroblast")

# Make data values table
df3 <- data.frame(group = sc_fibro[[group_metacolumn]],
                  celltype = sc_fibro[[gs_name]],
                  gene = Seurat::FetchData(object = sc_fibro, slot = "data", vars = c("Adam10", "Adam17", "Dlk1")))

# Density plot
pdf("Plots/5_Dobie_Scatter_Dlk1_Adam17.pdf", width = 7.5, height= 6.5)
FeatureScatter(sc_fibro, "Adam17", "Dlk1", pt.size = 2, group.by = gs_name, cols = "#CB1B4FFF")+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s Adam17+Dlk1+, %s Adam17+ only, %s Dlk1+ only", sum(df3$gene.Adam17 != 0 & df3$gene.Dlk1 != 0), sum(df3$gene.Adam17 != 0 & df3$gene.Dlk1 == 0), sum(df3$gene.Adam17 == 0 & df3$gene.Dlk1 != 0))))
dev.off()

pdf("Plots/5_Dobie_Scatter_Dlk1_Adam10.pdf", width = 7.5, height= 6.5)
FeatureScatter(sc_fibro, "Adam10", "Dlk1", pt.size = 2, group.by = gs_name, cols = "#CB1B4FFF")+ labs(title = gs_name, caption = stringr::str_wrap(sprintf("n fibroblasts: %s Adam10+Dlk1+, %s Adam10+ only, %s Dlk1+ only", sum(df3$gene.Adam10 != 0 & df3$gene.Dlk1 != 0), sum(df3$gene.Adam10 != 0 & df3$gene.Dlk1 == 0), sum(df3$gene.Adam10 == 0 & df3$gene.Dlk1 != 0))))
dev.off()

# Make ggplot2
font_size <- 12
line_size <- 1
paper_author <- "Dobie"
out_dir <- "Plots"


for(gene_boxplot in c("Dlk1", "Adam10", "Adam17")){
  # Make data values table{
  df3 <- data.frame(group = sc[["Treatment"]],
                    celltype = sc$List2mesenchyme_scType,
                    gene = Seurat::FetchData(object = sc, slot = "data", vars = c(gene_boxplot))) 
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
  counts_total <- sc@meta.data[,c('List2mesenchyme_scType', "Treatment")] %>%
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
  ggplot(df3, aes(x = celltype, y = gene, color=group)) +
    geom_point(size=1.5,alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, jitter.height =  0)) +
    scale_color_manual(values = group_colors) +
    labs(
      title = gene_boxplot,
      subtitle = sprintf("%s et al dataset; Crossbar = median", paper_author),
      caption = stringr::str_wrap(paste(gene_boxplot, "/ Total cells:", counts_label)),
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
  ggsave(filename = sprintf("%s/5_%s_%s_%s.pdf", out_dir, paper_author, gs_name, gene_boxplot), width = 7.5, height = 5.5, units = "in")
}

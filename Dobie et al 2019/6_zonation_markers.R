# Look at zonation markers

library(kazutils) # devtools::install_github("kazeera/kazutils")
library(dplyr)
library(Seurat)
library(HGNChelper)
library(ggplot2)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Subset to only treated HSCs
sc3$Treatment <- plyr::mapvalues(sc3$Treatment, c("CCl4", "Healthy"), c("Treated", "Healthy"))
sc3 <- subset(sc, List2mesenchyme_scType == "HSC" & Treatment == "Treated")

# Define input gene set folder (we have to do this iteratively)
gs_folder <- "../../GeneLists/scType DB mouse"

# gs_file <- list.files(gs_folder)[1]
for(gs_file in c("HSCzonation.xlsx")){
  # Prepare gene set
  xl_file <- paste(gs_folder, gs_file,sep = "/")
  gs_list = gene_sets_prepare(xl_file, "Liver")
  
  gs_name <- gsub(" ", "", gs_file) %>% gsub(".xlsx", "", .) %>% paste0(., "_scType")
  
  # Run scType to assign cell types
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  es.max = sctype_score(scRNAseqData = sc3@assays$RNA@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  # es.max[1:6,1:5]
  # dim(es.max) # rows(cell type) = 15, columns(cell_ID)=33168
  
  # merge by cluster
  cL_results = do.call("rbind", lapply(unique(sc3$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sc3@meta.data[sc3$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc3$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
   
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  # 
  sc3@meta.data[,gs_name] = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    sc3@meta.data[sc3@meta.data$seurat_clusters == j, gs_name] = as.character(cl_type$type[1])
  }
  
  rm(es.max)
}
table(sc3$HSCzonation_scType)

# Central  Portal 
# 9    5423 
sc4 <- subset(sc3, Dlk1 > 0)

table(sc4$HSCzonation_scType)
#  this many treated HSCs are DLk1+
# Portal 
# 31 
jpeg(sprintf("%s/4_Dlk_HSCzonation_umap.jpeg", out_dir), height = 480, width = 480)
print(DimPlot(sc3, group.by = "HSCzonation_scType")+labs(caption="DLk1+ on total treated HSCs = Central 0/9, Portal 31/5423"))
dev.off()


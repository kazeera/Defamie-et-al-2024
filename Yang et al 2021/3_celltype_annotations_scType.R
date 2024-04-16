rm(list=ls()[ls()!= "sc"])

#' Goal - Reassign celltypes using simpler/better R package called "scType" (compared to "cellassign")
# Note: gs = gene set

# Load libraries and functions
library(kazutils) # devtools::install_github("kazeera/kazutils")
library(dplyr)
library(Seurat)
library(HGNChelper)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Load Seurat data
sc <- readRDS("GSE171904_2021-01-04_seurat.rds")

# Define input gene set folder (we have to do this iteratively)
gs_folder <- "../../GeneLists/scType DB mouse"

# gs_file <- list.files(gs_folder)[1]
for(gs_file in c("List 2 long no Bip.xlsx", "List 2 long.xlsx")){
  # Prepare gene set
  xl_file <- paste(gs_folder, gs_file,sep = "/")
  gs_list = gene_sets_prepare(xl_file, "Liver")
  
  gs_name <- gsub(" ", "", gs_file) %>% gsub(".xlsx", "", .) %>% paste0(., "_scType")
  
  # Run scType to assign cell types
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  es.max = sctype_score(scRNAseqData = sc@assays$RNA@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  es.max[1:6,1:5]
  dim(es.max) # rows(cell type) = 15, columns(cell_ID)=33168
  
  # merge by cluster
  cL_results = do.call("rbind", lapply(unique(sc$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sc@meta.data[sc$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  # 
  sc@meta.data[,gs_name] = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    sc@meta.data[sc@meta.data$seurat_clusters == j, gs_name] = as.character(cl_type$type[1])
  }
  
  rm(es.max)
} 

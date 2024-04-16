# No data processing needed - normalized Seurat object is published
# ie: GSE171904_2021-01-04_seurat.rds

# Run once to reorganize file structures and make Seurat objects
library(kazutils)
load_packages(c("dplyr", "openxlsx", "Seurat"))

# 1. Reorganize all files from GSE into seperate folders for each expt
# Define data directories
setwd("Data")

# Get unique expts
dirs <- list.dirs(recursive = FALSE)

# In each of these directories, rename the files
lapply(dirs, FUN = function(dir){
  # List all the files in the folder
  files <- list.files(dir)
  # Specify the order of files
  files <- paste(dir, files[c(grep("barcodes", files), grep("features", files), grep("matrix", files))], sep="/")
  # Make the new file names
  new.files <- paste(dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"), sep = "/")
  # Rename the files
  file.rename(files, new.files)
})

## 2. Make Seurat objects -----------------------------------------

# https://satijalab.org/seurat/v3.1/merge_vignette.html
options(max.print= 45, width = 10)

# Initialize list for Seurat objects
seurat_obs <- list()

# For each expt, create a Seurat object based on the cutoffs indicated in the paper
for(dir in dirs){
  # Label for expt
  label <- get_nth_part(dir, "/", 2)
  
  # Read in data
  sc.data <- Read10X(data.dir = dir) 
  
  # QC from paper: Genes expressed in fewer than three cells in a sample were excluded, 
  # as were cells that expressed fewer than 300 genes or mitochondrial gene content > 30% of the total UMI count. 
  sc <- CreateSeuratObject(counts = sc.data, project = label)
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  
  # Save RData to sample's directory
  seurat_obs <- c(seurat_obs, sc)
}
names(seurat_obs) <- get_nth_part(dirs, "/", 2)

# Rename cells using object names (ie. expt) as prefix
for (i in names(seurat_obs)) {
  seurat_obs[[i]] <- RenameCells(seurat_obs[[i]], add.cell.id = i)
}

# Merge function in Seurat v>3.0
sc <- merge(x = seurat_obs[[1]], y = seurat_obs[2:3])

# Add celltype metadata from paper / published Seurat object
setwd("..")
df <- readRDS("6_Yang_paper_celltype_metadata.rds")
sc <- Seurat::AddMetaData(sc, df)

# Add metadata
sc$condition <- get_nth_part(colnames(sc), "_", 1)

# Visualize QC metrics as a violin plot
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "condition")

# Subset to cells 
range(sc$nCount_RNA)
sc <- subset(sc, cells = rownames(df)) #percent.mt < 50 & nCount_RNA > 200 & nFeature_RNA > 3)
sum(!rownames(df) %in% colnames(sc)) # 32 not found

saveRDS(sc, file="1_merged_Seurat.rds")
 
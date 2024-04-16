# Run once to reorganize file structures and make Seurat objects
library(kazutils)
load_packages(c("dplyr", "openxlsx", "Seurat"))

# 1. Reorganize all files from GSE into seperate folders for each expt
# Define data directories
gse <- "GSE137720"
setwd(sprintf("Data/%s_RAW", gse))
files <- list.files()

# Get unique expts
expts <- paste(get_nth_part(files, "_", 1), get_nth_part(files, "_", 2), get_nth_part(files, "_", 3), sep = "_") %>%
  unique

# Make each expt into its own file
lapply(expts, function(expt){
  dir.create(expt)
  file.copy(files[grep(expt, files)], to = expt, copy.mode = T)
})
# Note: Delete original files

# Get 10x dir
dirs_10x <- list.dirs()[grep("10x", list.dirs())]

# In each of these directories, rename the files
lapply(dirs_10x[1], FUN = function(dir){
  # List all the files in the folder
  files <- list.files(path = dir)
  # Specify the order of files
  files <- paste(dir, files[c(grep("barcodes", files), grep("genes", files), grep("matrix", files))], sep="/")
  # Make the new file names
  new.files <- paste(dir, c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"), sep = "/")
  # Rename the files
  file.rename(files, new.files)
})

# Delete copied files from current directory

## 2. Make Seurat objects -----------------------------------------
setwd("..")

# https://satijalab.org/seurat/v3.1/merge_vignette.html
options(max.print= 45, width = 10)

# Directories
dirs_10x <- list.dirs()[grep("10x", list.dirs())]
out_dir <- create_folder("GSE137720_Seurat_objects")
save.image()

# Read annotations for files
df <- read.xlsx("GSE137720_info.xlsx", sheet=1)
rownames(df) <- paste(df$GSM, df$Experiment, sep="_")

# Create labels for each expt based on annotations
labels <- df[get_nth_part(dirs_10x, "/", 3), "Characteristics"] %>% 
  gsub("_liver tissue", "", .)
names(labels) <- get_nth_part(dirs_10x, "/", 3)

# Initialize list for Seurat objects
seurat_obs <- list()

# For each expt, create a Seurat object based on the cutoffs indicated in the paper
for(data_dir in dirs_10x){
  # Label for expt
  label <- labels[get_nth_part(data_dir, "/", 3)]
  
  # Read in data
  sc.data <- Read10X(data.dir = data_dir) 
  
  # QC from paper: Genes expressed in fewer than three cells in a sample were excluded, 
  # as were cells that expressed fewer than 300 genes or mitochondrial gene content > 30% of the total UMI count. 
  sc <- CreateSeuratObject(counts = sc.data, project = label, min.cells = 3, min.features = 300)
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  # VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  sc <- subset(sc, percent.mt < 30)
  
  # Save RData to sample's directory
  saveRDS(sc, file = sprintf("%s/%s.rds", out_dir, label))
  seurat_obs <- c(seurat_obs, sc)
}
names(seurat_obs) <- labels

# Save individual Seurat objects
expts <- names(seurat_obs)
save(seurat_obs,  expts, file="../seurat_obs.RData")

# Rename cells using object names (ie. expt) as prefix
for (i in names(seurat_obs)) {
  seurat_obs[[i]] <- RenameCells(seurat_obs[[i]], add.cell.id = i)
}

# Merge function in Seurat v>3.0
merged_combined <- merge(x = seurat_obs[[1]], y = seurat_obs[2:5])
# > tail(merged_combined@meta.data)
# orig.ident nCount_RNA nFeature_RNA percent.mt
# CCl4 6wk_PF_TTTGTCATCCCAGGTG-1 CCl4 6wk_PF      10903         3177   1.724296
# CCl4 6wk_PF_TTTGTCATCGCGATCG-1 CCl4 6wk_PF       4754         1960   4.795961
# CCl4 6wk_PF_TTTGTCATCGGATGGA-1 CCl4 6wk_PF       2232         1219   2.688172
# CCl4 6wk_PF_TTTGTCATCGGCTACG-1 CCl4 6wk_PF       3267         1384   7.499235
saveRDS(merged_combined, file="merged_samp_Seurat.rds")

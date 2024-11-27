## 1. Reorganize all files from GSE into separate folders for each expt -----------------------------------------
# Run once to reorganize file structures and make Seurat objects
library(kazutils)
load_packages(c("dplyr", "openxlsx", "Seurat"))

# # Counts data
# counts <- read.csv("~/liver/read_count_mSTRT-seq_mouse_fetal_erythrocytes.csv")
# # Metadata - more than one sheet
# metadata <- read.xlsx("~/liver/cell_metadata.xlsx") 

# Define data directories
main <- "~/liver/"
main_dirs <- list.dirs(main, recursive = F)

# Get unique expts
dirs <- list.dirs(main_dirs[grep("mouse|w", main_dirs)], recursive = FALSE)

# In each of these directories, rename the files
lapply(dirs, FUN = function(dir){
  # # Get data directory
  # dir <- paste0(dir, "/rawdata")
  # List all the files in the folder
  files <- list.files(dir)
  # Specify the order of files
  files <- paste(dir, files[c(grep("barcodes", files), grep("genes", files), grep("matrix", files))], sep="/")
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
  label <- get_nth_part(dir, "/", 10)
  label <- gsub("mouse ", "", label)
  # Read in data
  sc.data <- Read10X(data.dir = dir) 
  # Create object
  sc <- CreateSeuratObject(counts = sc.data, project = label)
  
  # Label mitochondrial gene content percentage of the total UMI count. 
  if(grepl("mouse", dir)){
    sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  } else {
    sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  }
  
  # Save object
  seurat_obs <- c(seurat_obs, sc)
}
names(seurat_obs) <- get_nth_part(dirs, "/", 10) %>% gsub("mouse ", "", .)

# Rename cells using object names (ie. expt) as prefix
for (i in names(seurat_obs)) {
  seurat_obs[[i]] <- RenameCells(seurat_obs[[i]], add.cell.id = i)
}
saveRDS(seurat_obs, file="~/liver/seurat_obs.RData")

# Merge function in Seurat v>3.0
sc_human <- merge(x = seurat_obs[[1]], y = seurat_obs[2:9])
sc_mouse <- merge(x = seurat_obs[[10]], y = seurat_obs[11:12])
# Add metadata
sc_human$timepoint <- get_nth_part(colnames(sc_human), "_", 1)
sc_mouse$timepoint <- get_nth_part(colnames(sc_mouse), "_", 1)

# Add celltype metadata from paper / published Seurat object
# Human metadata
human_df <- read.xlsx("~/liver/cell_metadata.xlsx", sheet = "Human, 10x genomics", startRow = 2) 
rownames(human_df) <- gsub("W", "", human_df$SampleName)
rownames(human_df) <- gsub("_", "w_", rownames(human_df))
rownames(human_df) <- paste0(rownames(human_df), "-1")
all(colnames(sc_human) %in% rownames(human_df))
sc_human <- Seurat::AddMetaData(sc_human, human_df)

# Mouse metadata
mouse_df <- read.xlsx("~/liver/cell_metadata.xlsx", sheet = "Mouse, 10x genomics", startRow = 2) 
rownames(mouse_df) <- gsub(".0", "", mouse_df$SampleName)
rownames(mouse_df) <- paste0(rownames(mouse_df), "-1")
all(colnames(sc_mouse) %in% rownames(mouse_df))
sc_mouse <- Seurat::AddMetaData(sc_mouse, mouse_df)

# Save data
saveRDS(sc_mouse, file="~/liver/sc_mouse.rds")
saveRDS(sc_human, file="~/liver/sc_human.rds")

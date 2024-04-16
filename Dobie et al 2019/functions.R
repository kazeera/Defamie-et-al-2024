# Check whether there are rows that are all NAs (is a problem!!)
check_rows_NA <- function(df){
  all_NA_indices <- apply(df, 1, function(x) all(is.na(x)))
  # x <- sprintf("any NA: %s", any(is.na(all_NA_indices)))
  x <- sprintf("any all rows NA: %s", any(all_NA_indices, na.rm = T))
  print(x)
}

install_Bioconductor_pkg <- function(pkg){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(pkg)
}
#' Provide a list of file names and regular expression to get the file name
get_file <- function(file_list, regex){
  file_list[grep(regex, file_list)]
}

# Make feature plot
make_feature_plot <- function(sc, genes, file_label = "", label=F, rdn = "umap", lbl = F, pt.size = 0.6, ncol = 3, h = 6, w = 9, out_dir = ".", analysis_name ="", PDF = T){
  if(PDF){
    pdf(sprintf("%s/feature_%s.pdf", out_dir, file_label), height = h, width = w)
  } else {
    png(sprintf("%s/feature_%s.png", out_dir, file_label), height = h, width = w, units = "in", res = 500)
  }
  print(FeaturePlot(sc, reduction = rdn, label = lbl, pt.size = pt.size, ncol = ncol, features = genes))
  dev.off()
}

# Make heatmap
make_heatmap <- function(sc, genes, file_label="", out_dir=".", PDF = T){
  if(PDF){
    pdf(sprintf("%s/heatmap_%s.pdf", out_dir, file_label))
  } else {
    png(sprintf("%s/heatmap_%s.png", out_dir, file_label))
  }
  print(DoHeatmap(sc, features = genes))
  dev.off()
}

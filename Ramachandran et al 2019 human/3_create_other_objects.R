# -------------------------
# Mesenchyme signature:
markers <- list(mesenchyme_Ramachandran = c("PDGFRA",	"ACTA2", "COL1A1", "COL1A2", "COL3A1", "DES", "DCN"),
                mesenchyme_alt = c("PDGFRA", "VIM", "GATA4", "DES"))
saveRDS(markers, file="markers_mesenchymal.rds")

# Mesothelial signature:
saveRDS(c("GPM6A", "KRT18", "PDGFRA", "DES", "WT1", "CD200", "PDPN"), 
        file="markers_mesothelial.rds")
# -------------------------

# Mesenchymal signatures subtypes in paper -------------------------
# TODO - long to wide format so each celltype in its own column
df <- read.csv("../Ramachandran et al. 2019/Paper/Ramachandran_Marker genes for unsupervised clustering of mesenchymal cells.csv")

# Cell types
cts <- df$cluster %>% unique

# Gene markers
markers_h <- lapply(cts, function(x){
  df$gene[df$cluster == x]
})
names(markers_h) <- gsub("Mesenchyme \\(|\\)", "", cts)
saveRDS(markers_h, file="markers_Ramachandran_mesenchymal_subclasses_human.rds")

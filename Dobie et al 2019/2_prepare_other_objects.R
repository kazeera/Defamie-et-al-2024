# informatics.jax.org
# Hi, These are few hepatoblast markers : 
# Dlk1, Epcam+ (Epcam), CD133+ (Prom1, Cd133), E cadherin (Cdh1, Ecad), Nope+ (Nope, Igdcc4),
# CD13+ (Anpep, Cd13), Liv2, Hnf4?? (Hnf4a), T??RII (Tgfbr1), 
# CVK19 (Krt19), Sox9, osteopontin (Spp1), MIC1-1C3 (oval cell marker??), 
# Trop2 (Tacstd2), Foxl1, Lgr5
# CD133
# Save vectors of genes of interest to an RData file ------------
## Define genes of interest
markers <- c("Dlk1", "Epcam", "Prom1", "Cdh1", "Igdcc4", "Anpep", "Liv2", "Hnf4a", "Tgfbr1", "Krt19", "Sox9", "Spp1", "MIC1-1C3", "Tacstd2", "Foxl1", "Lgr5")
markers[!markers %in% rownames(data)] # "Epcam"    "Prom1"    "Cdh1"     "Liv2"     "MIC1-1C3" "Tacstd2" 


saveRDS(markers, file="hepatoblast_markers.rds")


# Zonation markers
library(openxlsx)
library(dplyr)

markers_zonation <- read.xlsx("Paper/1-s2.0-S2211124719313245-mmc2.xlsx", "markers_zonation") %>%
  lapply(function(x) x[!is.na(x)])
save(markers_zonation, file="Dobie2019_markers_zonation.RData")

# FB, HSC, VSMC signatures - only take top 10?
# markers_celltypes <- read.xlsx("Paper/1-s2.0-S2211124719313245-mmc2.xlsx", "markers_uninjured") 
# markers_celltypes %>%
#   group_by(cluster) %>%
#   .[,"gene"]
# Zonation markers
library(openxlsx)
library(dplyr)

df <- read.xlsx(xlsxFile = "D:/Kazeera/0000 Current Projects/Liver/20210129 sc Dobie et al 2019/Paper/1-s2.0-S2211124719313245-mmc2.xlsx", "markers_uninjured+fibrotic") 

# Cell types
cts <- df$cluster %>% unique
# table(df$cluster)
# FB  HSC VSMC 
# 149  294  176 

# Gene markers
markers_mesenchymal <- lapply(cts, function(x){
  df$gene[df$cluster == x]
})
names(markers_mesenchymal) <- cts
save(markers_mesenchymal, file="Dobie2019_markers_mesenchymal.RData")



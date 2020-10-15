rm(list = ls())

library(ape)
library(picante)
library(phytools)
library(apTreeshape)
library(raster)

source('./scripts/all_functions.R')

load('./data/continent/continent_new.Rdata')

# reading in the tree
shark_tree <- read.nexus('./data/611_taxa_tree topology for Emmaline.txt')

# creating species list from range data
species_names_poly <- dir('./data/polygon')
species_names_poly <- sub(pattern = '.json', "", species_names_poly)

# removing unnecessary tips
shark_tips <- shark_tree$tip.label
shark_bin <- as.character(sapply(sapply(shark_tips, function(x) strsplit(x, '_', 
                                                                         fixed = T)), function(y) 
                                                                           paste(y[[1]], y[[2]])))
shark_bin_clean <- gsub("'", '', shark_bin, fixed = T)
shark_bin_clean <- gsub('.', '', shark_bin_clean, fixed = T)
shark_tree$tip.label <- shark_bin_clean
shark_tree
shark_bin_discard <- grep('[2-5]', shark_bin_clean, value = T)
cf_names <- grep('cf', shark_tree$tip.label, value = T)
shark_tree_clean <- drop.tip(shark_tree, shark_bin_discard)
shark_tree_clean$tip.label = sub('cf', '', shark_tree_clean$tip.label)
sum(species_names_poly %in% shark_tree_clean$tip.label)
drop_names <- shark_tree_clean$tip.label[!(shark_tree_clean$tip.label 
                                           %in% species_names_poly)]
shark_tree_clean <- drop.tip(shark_tree_clean, drop_names)

# averaging duplicates and dropping one set
duplicate_tips <- which(duplicated(shark_tree_clean$tip.label))
species <- shark_tree_clean$tip.label[duplicate_tips]
edge_lengths <- shark_tree_clean$edge.length[duplicate_tips]
edge_df <- data.frame(species, edge_lengths)
edge_lengths2 <- shark_tree_clean$edge.length[duplicate_tips-1]
edge_df[,3] <- edge_lengths2
avg_edge <- rowMeans(edge_df[,2:3])
shark_tree_clean$edge.length[duplicate_tips-1] <- avg_edge
shark_tree_clean <- drop.tip(shark_tree_clean, duplicate_tips)
save(shark_tree_clean, file = './data/shark_tree_clean.Rdata')

# Applying mini_tree function
# test mini_tree using iucn
iucn_tree <- mini_tree('./data/IUCN')

# Carcharhiniforme tree
car_tree <- mini_tree('./data/Carcharhiniformes')

# Lamniforme tree
lam_tree <- mini_tree('./data/Lamniformes')

# Applying com_mat function in order to get mrd and psv
# total tree community matrix
load('./data/raster/sp_res_stack.Rdata')
sp_mat_list <- com_mat('./data/polygon', shark_tree_clean, sp_res_stack)

# Carcharhiniforme community matrix
load('./data/raster/car_res_stack.Rdata')
car_mat_list <- com_mat('./data/Carcharhiniformes', car_tree, car_res_stack)

# Lamniforme community matrix
load('./data/raster/lam_res_stack.Rdata')
lam_mat_list <- com_mat('./data/Lamniformes', lam_tree, lam_res_stack)

# Applying mrd_runplot for initial use
# Note: mrd rasters must be divided by their specific max value in order to standardize
# mrd_runplot test on total species
mrd_raster_list <- mrd_runplot(shark_tree_clean, sp_mat_list, continents,
                               './figures/mrd_rasters.pdf', species_richness)
save(mrd_raster_list, file = './data/raster/mrd_raster_list.Rdata')
load('./data/raster/mrd_raster_list.Rdata')

# mrd_runplot Carcharhiniformes
car_mrd_list <- mrd_runplot(car_tree, car_mat_list, continents,
                               './figures/Carcharhiniforme_mrd.pdf', species_richness)
save(car_mrd_list, file = './data/raster/car_mrd_list.Rdata')
load('./data/raster/car_mrd_list.Rdata')

# mrd_runplot Lamniformes
lam_mrd_list <- mrd_runplot(lam_tree, lam_mat_list, continents,
                               './figures/Lamniforme_mrd.pdf', species_richness)
save(lam_mrd_list, file = './data/raster/lam_mrd_list.Rdata')
load('./data/raster/lam_mrd_list.Rdata')

# dividing each mrd raster by its max to standardize the value
for (i in 1:6) {
max_wanted <- max(mrd_raster_list[[i]]@data@values, na.rm = T)
new_ras <- calc(mrd_raster_list[[i]], fun = function(x) x/max_wanted)
mrd_raster_list[[i]] <- new_ras
}
save(mrd_raster_list, file = './data/raster/mrd_raster_list.Rdata')

for (i in 1:6) {
  max_wanted <- max(car_mrd_list[[i]]@data@values, na.rm = T)
  new_ras <- calc(car_mrd_list[[i]], fun = function(x) x/max_wanted)
  car_mrd_list[[i]] <- new_ras
}
save(car_mrd_list, file = './data/raster/car_mrd_list.Rdata')

for (i in 1:6) {
  max_wanted <- max(lam_mrd_list[[i]]@data@values, na.rm = T)
  new_ras <- calc(lam_mrd_list[[i]], fun = function(x) x/max_wanted)
  lam_mrd_list[[i]] <- new_ras
}
save(lam_mrd_list, file = './data/raster/lam_mrd_list.Rdata')

pdf('./figures/mrd_rasters.pdf')
for (i in 1:6) {
  plot(mrd_raster_list[[i]])
  plot(continents, add = T, col = "black")
}
dev.off()

pdf('./figures/Carcharhiniforme_mrd.pdf')
for (i in 1:6) {
  plot(car_mrd_list[[i]])
  plot(continents, add = T, col = "black")
}
dev.off()

pdf('./figures/Lamniforme_mrd.pdf')
for (i in 1:6) {
  plot(lam_mrd_list[[i]])
  plot(continents, add = T, col = "black")
}
dev.off()

# Applying psv_runplot
# psv_runplot total species test
psv_raster_list <- psv_runplot(shark_tree_clean, sp_mat_list, continents,
                               './figures/psv_rasters.pdf', species_richness)
save(psv_raster_list, file = './data/raster/psv_raster_list.Rdata')
load('./data/raster/psv_raster_list.Rdata')

# psv_runplot Carcharhiniformes
car_psv_list <- psv_runplot(car_tree, car_mat_list, continents,
                            './figures/Carcharhiniforme_psv.pdf', species_richness)
save(car_psv_list, file = './data/raster/car_psv_list.Rdata')
load('./data/raster/car_psv_list.Rdata')

# psv_runplot Lamniformes 
lam_psv_list <- psv_runplot(lam_tree, lam_mat_list, continents,
                            './figures/Lamniforme_psv.pdf', species_richness)
save(lam_psv_list, file = './data/raster/lam_psv_list.Rdata')
load('./data/raster/lam_psv_list.Rdata')

# rooting the shark tree
plot(shark_tree_clean, show.tip.label = F)
#identify(shark_tree_clean)
#is.rooted(shark_tree_clean)
shark_tree_clean <- root(shark_tree_clean, outgroup = 257, resolve.root = T)
is.rooted(shark_tree_clean)
is.binary(shark_tree_clean)
shark_tree_clean <- multi2di(shark_tree_clean)

# beta statistic for entire tree
sharkb <- maxlik.betasplit(shark_tree_clean, confidence.interval = "profile")
sharkb

# beta across shark_tree_clean = -0.8868239 conf: -1.1103004:-0.621985

# Car_tree beta statistic
plot(car_tree, show.tip.label = F)
is.rooted(car_tree)
car_tree <- multi2di(car_tree)
car_b <- maxlik.betasplit(car_tree, confidence.interval = "profile")
car_b

# beta across Carcharhiniforme subtree = -0.9441983 conf: -1.2345835:-0.5738432

# Lam_tree beta statistic
plot(lam_tree, show.tip.label = F)
is.rooted(lam_tree)
lam_tree <- multi2di(lam_tree)
lam_b <- maxlik.betasplit(lam_tree, confidence.interval = "profile")
lam_b

# beta across Lamniformes subtree = -0.6364018 conf: -2:6.558514

# Applying beta_runplot
# total species beta_runplot
beta_raster_list <- beta_runplot(shark_tree_clean, sp_mat_list, continents,
                                 './figures/beta_rasters.pdf', species_richness)
save(beta_raster_list, file = './data/raster/beta_raster_list.Rdata')
load('./data/raster/beta_raster_list.Rdata')

library(ape)
library(phytools)
library(picante)
library(apTreeshape)

# reading in the tree
shark_tree <- read.nexus('./data/611_taxa_tree topology for Emmaline.txt')
shark_tree
plot(shark_tree)
plot(shark_tree, type ='f', cex=.25, show.tip.label = F)

# creating species list from range data
species_names_poly <- dir('./data/polygon')
species_names_poly <- sub(pattern = '.json', "", species_names_poly)

# removing unnecessary tips
shark_tips <- shark_tree$tip.label
shark_bin <- as.character(sapply(sapply(shark_tips, function(x) strsplit(x, '_', 
                                                                         fixed = T)), function(y) paste(y[[1]], y[[2]])))
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


# mini_tree function to make subtrees based on clades
# x = file directory to species polygons
# y = name of new tree
mini_tree <- function(x) {
   names_poly <- dir(x)
   names_poly <- sub(pattern = '.json', "", names_poly)
   keep <- shark_tree_clean$tip.label %in% names_poly
   drop <- shark_tree_clean$tip.label[!keep]
   new_tree <- drop.tip(shark_tree_clean, drop)
   return(new_tree)
}

# test mini_tree using iucn
iucn_tree <- mini_tree('./data/IUCN')

# Carcharhiniforme tree
car_tree <- mini_tree('./data/Carcharhiniformes')

# Lamniforme tree
lam_tree <- mini_tree('./data/Lamniformes')

# com_mat function makes community matrix for phylogenetic diversity
# x = file directory of species
# y = clade tree
# z = res_stack of clade
# output value is community matrix for specified clade, name accordingly
com_mat <- function(x, y, z) {
names_poly <- dir(x)
names_poly <- sub(pattern = '.json', "", names_poly)
test_species <- names_poly %in% y$tip.label
i <- 1
mat_list <- vector("list", length = length(z))
for (j in 1:6) {
  mat_list[[j]] = matrix(NA, ncol = sum(test_species), 
                         nrow=length(z[[j]][[i]]))
  colnames(mat_list[[j]]) = names_poly[test_species]
  icol = 1
  for (i in which(test_species)) {
    mat_list[[j]][ , icol] = extract(z[[j]][[i]], 
                                     1:nrow(mat_list[[j]]))
    icol = icol + 1
  }
  mat_list[[j]] = ifelse(is.na(mat_list[[j]]), 0, mat_list[[j]])
}
return(mat_list)
}

# total tree community matrix
sp_mat_list <- com_mat('./data/polygon', shark_tree_clean, sp_res_stack)

# Carcharhiniforme community matrix
car_mat_list <- com_mat('./data/Carcharhiniformes', car_tree, car_res_stack)

# Lamniforme community matrix
lam_mat_list <- com_mat('./data/Lamniformes', lam_tree, lam_res_stack)


# mrd_runplot function for mean root distance
# x = clade tree
# y = community matrix for clade
# z = file name for pdf
# output value is the mrd raster list, name and save accordingly
mrd_runplot <- function(x, y, z) {
phylo_bl1 <- compute.brlen(x, 1)
all_dist <- dist.nodes(phylo_bl1)
root_dist <- all_dist[length(x$tip.label) + 1, 
                      1:length(x$tip.label)]
tips_to_root <- data.frame(spp.name=x$tip.label, root_dist)

mrd_test_list <- vector("list", length = length(y)) 
for (i in 1:length(y)) {
  allspecies = colnames(y[[i]])
  mrd_test_list[[i]] = rep(0, nrow(y[[i]]))
  for(j in 1:nrow(y[[i]])) {
    sp_list = data.frame(spp.name = allspecies[y[[i]][j, ] == 1])
    if (nrow(sp_list) > 0) {
      root_dist_tot <- merge(sp_list, tips_to_root, sort = F)
      mrd_test_list[[i]][j] <- mean(root_dist_tot$root_dist)
    }
  }
}

load('./data/raster/species_richness.Rdata')
raster_list <- vector("list", length = 6)
pdf(z)
for (i in 1:6) {
  mrd_raster <- species_richness[[i]]
  mrd_raster@data@values <- mrd_test_list[[i]]
  plot(mrd_raster)
  plot(continents, add = T, col = "black")
  raster_list[[i]] <- mrd_raster
}
dev.off()
return(raster_list)
}

# mrd_runplot test on total species
mrd_raster_list <- mrd_runplot(shark_tree_clean, sp_mat_list, 
                               './figures/mrd_rasters.pdf')
save(mrd_raster_list, file = './data/raster/mrd_raster_list.Rdata')
load('./data/raster/mrd_raster_list.Rdata')

# mrd_runplot Carcharhiniformes
car_mrd_list <- mrd_runplot(car_tree, car_mat_list, 
                               './figures/Carcharhiniforme_mrd.pdf')
save(car_mrd_list, file = './data/raster/car_mrd_list.Rdata')
load('./data/raster/car_mrd_list.Rdata')

# mrd_runplot Lamniformes
lam_mrd_list <- mrd_runplot(lam_tree, lam_mat_list, 
                               './figures/Lamniforme_mrd.pdf')
save(lam_mrd_list, file = './data/raster/lam_mrd_list.Rdata')
load('./data/raster/lam_mrd_list.Rdata')

# psv_runplot function for phylogenetic species diversity 
# x = clade tree
# y = community matrix for clade
# z = pdf ile
# output is raster list, name and save accordingly
psv_runplot <- function(x, y, z) {
psv_test_list <- vector("list", length = 6)
for (i in 1:6) {
  psv_test <- psv(y[[i]], x)
  psv_test_list[[i]] <- psv_test
}

load('./data/raster/species_richness.Rdata')
pdf(z)
raster_list <- vector("list", length = 6)
for (i in 1:6) {
  psv_raster <- species_richness[[i]]
  psv_raster@data@values <- psv_test_list[[i]]$PSVs
  plot(psv_raster)
  plot(continents, add = T, col = "black")
  raster_list[[i]] <- psv_raster
}
dev.off()
return(raster_list)
}

# psv_runplot total species test
psv_raster_list <- psv_runplot(shark_tree_clean, sp_mat_list, 
                               './figures/psv_rasters.pdf')
save(psv_raster_list, file = './data/raster/psv_raster_list.Rdata')
load('./data/raster/psv_raster_list.Rdata')

# psv_runplot Carcharhiniformes
car_psv_list <- psv_runplot(car_tree, car_mat_list, 
                            './figures/Carcharhiniforme_psv.pdf')
save(car_psv_list, file = './data/raster/car_psv_list.Rdata')
load('./data/raster/car_psv_list.Rdata')

# psv_runplot Lamniformes 
lam_psv_list <- psv_runplot(lam_tree, lam_mat_list, 
                            './figures/Lamniforme_psv.pdf')
save(lam_psv_list, file = './data/raster/lam_psv_list.Rdata')
load('./data/raster/lam_psv_list.Rdata')

# rooting the shark tree
plot(shark_tree_clean, show.tip.label = F)
identify.phylo(shark_tree_clean)
is.rooted(shark_tree_clean)
shark_tree_clean <- root(shark_tree_clean, outgroup = 257, resolve.root = T)
is.rooted(shark_tree_clean)
is.binary(shark_tree_clean)
shark_tree_clean <- multi2di(shark_tree_clean)

# beta statistic for entire tree
sharkb <- maxlik.betasplit(shark_tree_clean, confidence.interval = "profile")
sharkb

# beta across shark_tree_clean = -0.8729644

# Car_tree beta statistic
plot(car_tree, show.tip.label = F)
is.rooted(car_tree)
car_tree <- multi2di(car_tree)
car_b <- maxlik.betasplit(car_tree, confidence.interval = "profile")
car_b

# beta across Carcharhiniforme subtree = -1.014664


# beta_runplot
# I'll write a function but I don't think I'll be able to apply this one to at 
# least the Lamniformes due to the small tree size
# x = clade tree
# y = community matrix
# z = pdf file
# output is beta raster list, name and save accordingly

beta_runplot <- function(x, y, z) {
beta_list <- vector("list", length = length(y))
for (i in 1:length(y)) {
  allspecies = colnames(y[[i]])
  beta_list[[i]] = rep(NA, nrow(y[[i]]))
  for(j in 1:nrow(y[[i]])) {
    sp_list = data.frame(spp.name = allspecies[y[[i]][j, ] == 1])
    if (nrow(sp_list) > 1) {
      drop_these <- x$tip.label[!(x$tip.label %in% sp_list$spp.name)]
      temp_tree <- drop.tip(x, drop_these)
      temp_tree <- multi2di(temp_tree)
      b <- maxlik.betasplit(temp_tree)
      beta_list[[i]][j] <- b$max_lik
    }
  }
}

load('./data/raster/species_richness.Rdata')
raster_list <- vector("list", length = 6)
pdf(z)
for (i in 1:6) {
  beta_raster <- species_richness[[i]]
  beta_raster@data@values <- beta_list[[i]]
  plot(beta_raster)
  plot(continents, add = T, col = "black")
  raster_list[[i]] <- beta_raster
}
dev.off()
return(raster_list)
}

# total species beta_runplot
beta_raster_list <- beta_runplot(shark_tree_clean, sp_mat_list, 
                                 './figures/beta_rasters.pdf')
save(beta_raster_list, file = './data/raster/beta_raster_list.Rdata')
load('./data/raster/beta_raster_list.Rdata')

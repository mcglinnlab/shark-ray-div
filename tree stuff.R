library(ape)
library(phytools)
library(picante)

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
drop_names <- shark_tree_clean$tip.label[!(shark_tree_clean$tip.label %in% species_names_poly)]
shark_tree_clean <- drop.tip(shark_tree_clean, drop_names)
shark_tree_clean$tip.label
genus_names <- sapply(strsplit(species_names_poly, ' '), function(x) x[[1]])
tree_genera <- sapply(strsplit(shark_tree_clean$tip.label, ' '), function(x) x[[1]])
sum(unique(genus_names) %in% unique(tree_genera))
length(unique(genus_names))
genus_names[!(genus_names %in% tree_genera)]

# rooting the shark tree
plot(shark_tree_clean, show.tip.label = F)
identify.phylo(shark_tree_clean)
is.rooted(shark_tree_clean)
shark_tree_clean <- root(shark_tree_clean, node = 575, resolve.root = T)
is.rooted(shark_tree_clean)
is.binary(shark_tree_clean)

# potential method to make tree ultrmetric
is.binary.tree(shark_tree_clean)
attempt_tree <- multi2di(shark_tree)
is.binary.tree(attempt_tree)
is.rooted(attempt_tree)
attempt_tree <- chronopl(attempt_tree, lambda = 0)
plot(attempt_tree, show.tip.label = F)

# Community matrix for phylogenetic diversity
test_species <- species_names_poly %in% shark_tree_clean$tip.label

i <- 1
mat_list <- vector("list", length = 7)
for (j in 1:7) {
    mat_list[[j]] = matrix(NA, ncol = sum(test_species), 
                            nrow=length(raster_res_list[[j]][[i]]@data@values))
    colnames(mat_list[[j]]) = species_names_poly[test_species]
    icol = 1
    for (i in which(test_species)) {
         mat_list[[j]][ , icol] = raster_res_list[[j]][[i]]@data@values
         icol = icol + 1
    }
    mat_list[[j]] = ifelse(is.na(mat_list[[j]]), 0, mat_list[[j]])
}

# Faith's phylogenetic diversity test
pd_test_list <- vector("list", length = 7) 
for (i in 1:7) {
  pd_test <- pd(mat_list[[i]], shark_tree_clean)
  pd_test_list[[i]] <- pd_test
}

# phylogenetic species diversity metrics
psv_test_list <- vector("list", length = 7)
for (i in 1:7) {
     psv_test <- psv(mat_list[[i]], shark_tree_clean)
     psv_test_list[[i]] <- psv_test
}

psr_test_list <- vector("list", length = 7)
for (i in 1:7) {
  psr_test <- psr(mat_list[[i]], shark_tree_clean)
  psr_test_list[[i]] <- psr_test
}

# phylogenetic diversity rasters
pdf('./figures/pd_rasters.pdf')
pd_raster_list <- vector("list", length = 7)
for (i in 1:7) {
     pd_raster <- raster_res_list[[i]][[1]]
     pd_raster@data@values <- pd_test_list[[i]]$PD
     pd_raster <- mask(pd_raster, mask_ras_list[[i]])
     plot(pd_raster)
     pd_raster_list[[i]] <- pd_raster
}
dev.off()

save(pd_raster_list, file = './data/raster/pd_raster_list.Rdata')
load('./data/raster/pd_raster_list.Rdata')

pdf('./figures/psv_rasters.pdf')
psv_raster_list <- vector("list", length = 7)
for (i in 1:7) {
     psv_raster <- raster_res_list[[i]][[1]]
     psv_raster@data@values <- psv_test_list[[i]]$PSVs
     psv_raster <- mask(psv_raster, mask_ras_list[[i]])
     plot(psv_raster)
     psv_raster_list[[i]] <- psv_raster
}
dev.off()

save(psv_raster_list, file = './data/raster/psv_raster_list.Rdata')
load('./data/raster/psv_raster_list.Rdata')

pdf('./figures/psr_rasters.pdf')
psr_raster_list <- vector("list", length = 7)
for (i in 1:7) {
  psr_raster <- raster_res_list[[i]][[1]]
  psr_raster@data@values <- psr_test_list[[i]]$PSR
  psr_raster <- mask(psr_raster, mask_ras_list[[i]])
  plot(psr_raster)
  psr_raster_list[[i]] <- psr_raster
}
dev.off()

save(psr_raster_list, file = './data/raster/psr_raster_list.Rdata')
load('./data/raster/psr_raster_list.Rdata')


# Tree randomization
randotree_list <- vector("list", length = 50)
for (i in 1:50) {
     randotree <- shark_tree_clean
     randomtips <- sample(randotree$tip.label, length(randotree$tip.label))
     randotree$tip.label <- randomtips
     randotree_list[[i]] <- randotree
}

# Faith's phylogenetic diversity test
random_pd_list <- vector("list", length = 7) 
for (i in 1:7) {
  pd_test <- pd(mat_list[[i]], randotree_list[[50]])
  random_pd_list[[i]] <- pd_test
}

# phylogenetic species diversity metrics
random_psv_list <- vector("list", length = 7)
for (i in 1:7) {
  psv_test <- psv(mat_list[[i]], randotree_list[[50]])
  random_psv_list[[i]] <- psv_test
}

random_psr_list <- vector("list", length = 7)
for (i in 1:7) {
  psr_test <- psr(mat_list[[i]], randotree_list[[50]])
  random_psr_list[[i]] <- psr_test
}

# Random phylogenetic diversity rasters
pdf('./figures/random_pd_rasters.pdf')
randompd_raster_list <- vector("list", length = 7)
for (i in 1:7) {
  pd_raster <- raster_res_list[[i]][[1]]
  pd_raster@data@values <- random_pd_list[[i]]$PD
  pd_raster <- mask(pd_raster, mask_ras_list[[i]])
  plot(pd_raster)
  randompd_raster_list[[i]] <- pd_raster
}
dev.off()

save(randompd_raster_list, file = './data/raster/randompd_raster_list.Rdata')
load('./data/raster/randompd_raster_list.Rdata')

pdf('./figures/random_psv_rasters.pdf')
randompsv_raster_list <- vector("list", length = 7)
for (i in 1:7) {
  psv_raster <- raster_res_list[[i]][[1]]
  psv_raster@data@values <- random_psv_list[[i]]$PSVs
  psv_raster <- mask(psv_raster, mask_ras_list[[i]])
  plot(psv_raster)
  randompsv_raster_list[[i]] <- psv_raster
}
dev.off()

save(randompsv_raster_list, file = './data/raster/randompsv_raster_list.Rdata')
load('./data/raster/randompsv_raster_list.Rdata')

pdf('./figures/random_psr_rasters.pdf')
randompsr_raster_list <- vector("list", length = 7)
for (i in 1:7) {
  psr_raster <- raster_res_list[[i]][[1]]
  psr_raster@data@values <- random_psr_list[[i]]$PSR
  psr_raster <- mask(psr_raster, mask_ras_list[[i]])
  plot(psr_raster)
  randompsr_raster_list[[i]] <- psr_raster
}
dev.off()

save(randompsr_raster_list, file = './data/raster/randompsr_raster_list.Rdata')
load('./data/raster/randompsr_raster_list.Rdata')
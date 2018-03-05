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

shark_tree_clean$tip.label
genus_names <- sapply(strsplit(species_names_poly, ' '), function(x) x[[1]])
tree_genera <- sapply(strsplit(shark_tree_clean$tip.label, ' '), function(x) x[[1]])
sum(unique(genus_names) %in% unique(tree_genera))
length(unique(genus_names))
genus_names[!(genus_names %in% tree_genera)]

# Community matrix for phylogenetic diversity
test_species <- species_names_poly %in% shark_tree_clean$tip.label

i <- 1
mat_list <- vector("list", length = length(sp_res_stack))
for (j in 1:6) {
    mat_list[[j]] = matrix(NA, ncol = sum(test_species), 
                            nrow=length(sp_res_stack[[j]][[i]]@data@values))
    colnames(mat_list[[j]]) = species_names_poly[test_species]
    icol = 1
    for (i in which(test_species)) {
         mat_list[[j]][ , icol] = sp_res_stack[[j]][[i]]@data@values
         icol = icol + 1
    }
    mat_list[[j]] = ifelse(is.na(mat_list[[j]]), 0, mat_list[[j]])
}

# Mean root distance
phylo_bl1 <- compute.brlen(shark_tree_clean, 1)
all_dist <- dist.nodes(phylo_bl1)
root_dist <- all_dist[length(shark_tree_clean$tip.label) + 1, 
                      1:length(shark_tree_clean$tip.label)]
tips_to_root <- data.frame(spp.name=shark_tree_clean$tip.label, root_dist)

mrd_test_list <- vector("list", length = length(mat_list)) 
for (i in 1:length(mat_list)) {
    allspecies = colnames(mat_list[[i]])
    mrd_test_list[[i]] = rep(0, nrow(mat_list[[i]]))
    for(j in 1:nrow(mat_list[[i]])) {
        sp_list = data.frame(spp.name = allspecies[mat_list[[i]][j, ] == 1])
        if (nrow(sp_list) > 0) {
            root_dist_tot <- merge(sp_list, tips_to_root, sort = F)
            mrd_test_list[[i]][j] <- mean(root_dist_tot$root_dist)
        }
    }
}

#MRD rasters
mrd_raster_list <- vector("list", length = 6)
pdf('./figures/mrd_rasters.pdf')
for (i in 1:6) {
     mrd_raster <- species_richness[[i]]
     mrd_raster@data@values <- mrd_test_list[[i]]
     plot(mrd_raster)
     plot(continents, add = T, col = "black")
     mrd_raster_list[[i]] <- mrd_raster
}
dev.off()

save(mrd_raster_list, file = './data/raster/mrd_raster_list.Rdata')
load('./data/raster/mrd_raster_list.Rdata')

# Faith's phylogenetic diversity test
pd_test_list <- vector("list", length = 6) 
for (i in 1:6) {
     pd_test <- pd(mat_list[[i]], shark_tree_clean)
     pd_test_list[[i]] <- pd_test
}

# phylogenetic species diversity metrics
psv_test_list <- vector("list", length = 6)
for (i in 1:6) {
     psv_test <- psv(mat_list[[i]], shark_tree_clean)
     psv_test_list[[i]] <- psv_test
}

psr_test_list <- vector("list", length = 6)
for (i in 1:6) {
  psr_test <- psr(mat_list[[i]], shark_tree_clean)
  psr_test_list[[i]] <- psr_test
}

# phylogenetic diversity rasters
pdf('./figures/pd_rasters.pdf')
pd_raster_list <- vector("list", length = 6)
for (i in 1:6) {
     pd_raster <- sp_res_stack[[i]][[1]]
     pd_raster@data@values <- pd_test_list[[i]]$PD
     plot(pd_raster)
     plot(continents, add = T, col = "black")
     pd_raster_list[[i]] <- pd_raster
}
dev.off()

save(pd_raster_list, file = './data/raster/pd_raster_list.Rdata')
load('./data/raster/pd_raster_list.Rdata')

pdf('./figures/psv_rasters.pdf')
psv_raster_list <- vector("list", length = 6)
for (i in 1:6) {
     psv_raster <- sp_res_stack[[i]][[1]]
     psv_raster@data@values <- psv_test_list[[i]]$PSVs
     plot(psv_raster)
     plot(continents, add = T, col = "black")
     psv_raster_list[[i]] <- psv_raster
}
dev.off()

save(psv_raster_list, file = './data/raster/psv_raster_list.Rdata')
load('./data/raster/psv_raster_list.Rdata')

pdf('./figures/psr_rasters.pdf')
psr_raster_list <- vector("list", length = 6)
for (i in 1:6) {
  psr_raster <- sp_res_stack[[i]][[1]]
  psr_raster@data@values <- psr_test_list[[i]]$PSR
  plot(psr_raster)
  plot(continents, add = T, col = "black")
  psr_raster_list[[i]] <- psr_raster
}
dev.off()

save(psr_raster_list, file = './data/raster/psr_raster_list.Rdata')
load('./data/raster/psr_raster_list.Rdata')

# rooting the shark tree
plot(shark_tree_clean, show.tip.label = F)
identify.phylo(shark_tree_clean)
is.rooted(shark_tree_clean)
shark_tree_clean <- root(shark_tree_clean, outgroup = 257, resolve.root = T)
is.rooted(shark_tree_clean)
is.binary(shark_tree_clean)
shark_tree_clean <- multi2di(shark_tree_clean)

# beta statistic
sharkb <- maxlik.betasplit(shark_tree_clean, confidence.interval = "profile")
sharkb

i <- 1
j <- 9783


beta_list <- vector("list", length = length(mat_list))
for (i in 1:length(mat_list)) {
     allspecies = colnames(mat_list[[i]])
     beta_list[[i]] = rep(NA, nrow(mat_list[[i]]))
     for(j in 1:nrow(mat_list[[i]])) {
         sp_list = data.frame(spp.name = allspecies[mat_list[[i]][j, ] == 1])
         if (nrow(sp_list) > 1) {
            drop_these <- shark_tree_clean$tip.label[!(shark_tree_clean$tip.label 
                                                    %in% sp_list$spp.name)]
            temp_tree <- drop.tip(shark_tree_clean, drop_these)
            temp_tree <- multi2di(temp_tree)
            b <- maxlik.betasplit(temp_tree)
            beta_list[[i]][j] <- b$max_lik
         }
    }
}

save(beta_list, file = './data/raster/beta_list.Rdata')

#beta rasters
beta_raster_list <- vector("list", length = length(raster_res_list))
pdf('./figures/beta_rasters.pdf')
for (i in 1:9) {
  beta_raster <- raster_res_list[[i]][[1]]
  beta_raster@data@values <- beta_list[[i]]
  plot(beta_raster)
  plot(continents, add = T, col = "black")
  beta_raster_list[[i]] <- beta_raster
}
dev.off()

save(beta_raster_list, file = './data/raster/beta_raster_list.Rdata')
load('./data/raster/beta_raster_list.Rdata')





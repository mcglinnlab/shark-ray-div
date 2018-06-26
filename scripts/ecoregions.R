library(rgdal)
library(sp)
library(maps)
library(maptools)
library(raster)

# creating and transforming realms
regions <- readOGR(dsn = "./data/continent", layer = "meow_ecos")
class(regions)
plot(regions)
regions$REALM
pull1 <- grep("Tropical Atlantic", regions$REALM)
tropical_atlantic <- regions[pull1,]
plot(tropical_atlantic)
pull2 <- grep("Central Indo-Pacific", regions$REALM)
central_indopacific <- regions[pull2,]
plot(central_indopacific)
tropical_atlantic <- spTransform(tropical_atlantic, CRS("+proj=cea +units=km"))
central_indopacific <- spTransform(central_indopacific, CRS("+proj=cea +units=km"))

# cropping species richness based on realms

# realm_crop function
# x = output list
# y = raster list
# z = realm spatial object
realm_crop <- function(x, y) {
  temporary <- vector("list", length = 6)
  for (i in 1:6) {
    ras <- intersect(x[[i]], y)
    temporary[[i]] <- ras
  }
  return(temporary)
}

trop_atlantic_list <- realm_crop(species_richness, tropical_atlantic)
save(trop_atlantic_list, file = './data/raster/trop_atlantic_list.Rdata')
load('./data/raster/trop_atlantic_list.Rdata')

trop_atlantic_stack <- realm_crop(sp_res_stack, tropical_atlantic)
save(trop_atlantic_stack, file = './data/raster/trop_atlantic_stack.Rdata')
load('./data/raster/trop_atlantic_stack.Rdata')

indopacific_list <- realm_crop(species_richness, central_indopacific)
save(indopacific_list, file = './data/raster/indopacific_list.Rdata')
load('./data/raster/indopacific_list.Rdata')

indopacific_stack <- realm_crop(sp_res_stack, central_indopacific)
save(indopacific_stack, file = './data/raster/indopacific_stack.Rdata')
load('./data/raster/indopacific_stack.Rdata')

# creating dir
# Pull_names function to pull out names of species in the realm
# x = res_stack
pull_names <- function(x) {
names_list <- vector("list", length = 6)
s <- vector("list", length = 6)
for (i in 1:6) {
  for (j in 1:534) {
    #print(any(trop_atlantic_stack[[i]][[j]]@data@values==1))
    if (sum(x[[i]][[j]]@data@values==1,na.rm=T)>0)
      {
      names_list[[i]][[j]] <- x[[i]][[j]]@data@names
      } else {names_list[[i]][[j]] <- NA}
  }
  s[[i]] <- names_list[[i]][which(!(is.na(names_list[[i]])))]
}
return(s)
}

trop_atlantic_names <- pull_names(trop_atlantic_stack)
indopacific_names <- pull_names(indopacific_stack)
    
dir.create('./data/tropical_atlantic')
dir.create('./data/indopacific')

ta_name_sub <- lapply(trop_atlantic_names[[6]], function(x)
  paste0(x, '.json'))
ta_name_sub <- gsub('_', ' ', ta_name_sub)
for (i in ta_name_sub) {
  file.copy(from = paste0('./data/polygon/', i), to = 'data/tropical_atlantic')
}

ip_name_sub <- lapply(indopacific_names[[1]], function(x)
  paste0(x, '.json'))
ip_name_sub <- gsub('_', ' ', ip_name_sub)
for (i in ip_name_sub) {
  file.copy(from = paste0('./data/polygon/', i), to = 'data/indopacific')
}

# phylogenetic metrics for realms

source('./scripts/Phylogenetic_metrics.R')

trop_atlantic_tree <- mini_tree('./data/tropical_atlantic')
indopacific_tree <- mini_tree('./data/indopacific')

trop_atlantic_mat <- com_mat('./data/tropical_atlantic', trop_atlantic_tree, 
                             trop_atlantic_stack)
indopacific_mat <- com_mat('./data/indopacific', indopacific_tree, 
                           indopacific_stack)

# mean root distance for ecoregions
  phylo_bl1 <- compute.brlen(trop_atlantic_tree, 1)
  all_dist <- dist.nodes(phylo_bl1)
  root_dist <- all_dist[length(trop_atlantic_tree$tip.label) + 1, 
                        1:length(trop_atlantic_tree$tip.label)]
  tips_to_root <- data.frame(spp.name=trop_atlantic_tree$tip.label, root_dist)
  
  mrd_test_list <- vector("list", length = length(trop_atlantic_mat)) 
  for (i in 1:length(trop_atlantic_mat)) {
    allspecies = colnames(trop_atlantic_mat[[i]])
    mrd_test_list[[i]] = rep(0, nrow(trop_atlantic_mat[[i]]))
    for(j in 1:nrow(trop_atlantic_mat[[i]])) {
      sp_list = data.frame(spp.name = allspecies[trop_atlantic_mat[[i]][j, ] == 1])
      if (nrow(sp_list) > 0) {
        root_dist_tot <- merge(sp_list, tips_to_root, sort = F)
        mrd_test_list[[i]][j] <- mean(root_dist_tot$root_dist)
      }
    }
  }
  
  load('./data/raster/trop_atlantic_list.Rdata')
  trop_atlantic_mrd <- vector("list", length = 6)
  pdf('./figures/trop_atlantic_mrd.pdf')
  for (i in 1:6) {
    mrd_raster <- trop_atlantic_list[[i]]
    mrd_raster@data@values <- mrd_test_list[[i]]
    plot(mrd_raster)
    plot(continents, add = T, col = "black")
    trop_atlantic_mrd[[i]] <- mrd_raster
  }
  dev.off()

save(trop_atlantic_mrd, file = './data/raster/trop_atlantic_mrd.Rdata')
load('./data/raster/trop_atlantic_mrd.Rdata')

phylo_bl1 <- compute.brlen(indopacific_tree, 1)
all_dist <- dist.nodes(phylo_bl1)
root_dist <- all_dist[length(indopacific_tree$tip.label) + 1, 
                      1:length(indopacific_tree$tip.label)]
tips_to_root <- data.frame(spp.name=indopacific_tree$tip.label, root_dist)

mrd_test_list <- vector("list", length = length(indopacific_mat)) 
for (i in 1:length(indopacific_mat)) {
  allspecies = colnames(indopacific_mat[[i]])
  mrd_test_list[[i]] = rep(0, nrow(indopacific_mat[[i]]))
  for(j in 1:nrow(indopacific_mat[[i]])) {
    sp_list = data.frame(spp.name = allspecies[indopacific_mat[[i]][j, ] == 1])
    if (nrow(sp_list) > 0) {
      root_dist_tot <- merge(sp_list, tips_to_root, sort = F)
      mrd_test_list[[i]][j] <- mean(root_dist_tot$root_dist)
    }
  }
}

load('./data/raster/indopacific_list.Rdata')
indopacific_mrd <- vector("list", length = 6)
pdf('./figures/indopacific_mrd.pdf')
for (i in 1:6) {
  mrd_raster <- indopacific_list[[i]]
  mrd_raster@data@values <- mrd_test_list[[i]]
  plot(mrd_raster)
  plot(continents, add = T, col = "black")
  indopacific_mrd[[i]] <- mrd_raster
}
dev.off()

indopacific_mrd <- mrd_runplot(indopacific_tree, indopacific_mat, 
                               './figures/indopacific_mrd.pdf')

save(indopacific_mrd, file = './data/raster/indopacific_mrd.Rdata')
load('./data/raster/indopacific_mrd.Rdata')
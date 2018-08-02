library(rgdal)
library(sp)
library(maps)
library(maptools)
library(raster)
library(rgeos)
library(doParallel)
library(foreach)
library(ape)
library(phytools)
library(picante)
library(apTreeshape)

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
central_indopacific <- spTransform(central_indopacific, CRS("+proj=cea +lat_0=0 +lon_0=125 +units=km"))

# may have to drop these ecoregions
#drop1 <- grep("Tonga Islands", regions$ECOREGION)
#drop1 <- regions[drop1,]
#drop2 <- grep("Fiji Islands", regions$ECOREGION)
#drop2 <- regions[drop2,]
#drop_indo <- rbind(drop1, drop2)
#test <- gDifference(central_indopacific, drop_indo)
#test <- spTransform(test, CRS("+proj=cea +units=km"))

# cropping species richness based on realms

# realm_crop function
# x = output list
# y = realm spatial object
realm_crop <- function(x, y) {
  temporary <- vector("list", length = 6)
  for (i in 1:6) {
    ras <- intersect(x[[i]], y)
    temporary[[i]] <- ras
  }
  return(temporary)
}

load('./data/raster/species_richness.Rdata')
load('./data/raster/sp_res_stack.Rdata')
continents <- shapefile('./data/continent/continent')
continents <- spTransform(continents, CRS("+proj=cea +units=km"))
continents_indo <- spTransform(continents, 
                               CRS("+proj=cea +lat_0=0 +lon_0=125 +units=km"))
oceans <- readOGR(dsn = "./data/Environment", layer = "ne_10m_ocean")
oceans <- spTransform(oceans, CRS("+proj=cea +lat_0=0 +lon_0=125 +units=km"))
oceans_indo <- raster(oceans)
res(oceans_indo) <- 110
save(oceans_indo, file = './data/raster/oceans_indo.Rdata')


trop_atlantic_list <- realm_crop(species_richness, tropical_atlantic)
save(trop_atlantic_list, file = './data/raster/trop_atlantic_list.Rdata')
load('./data/raster/trop_atlantic_list.Rdata')

trop_atlantic_stack <- realm_crop(sp_res_stack, tropical_atlantic)
save(trop_atlantic_stack, file = './data/raster/trop_atlantic_stack.Rdata')
load('./data/raster/trop_atlantic_stack.Rdata')

pdf('./figures/tropical_atlantic_richness.pdf')
for (i in 1:6) {
  plot(trop_atlantic_list[[i]])
  plot(continents, add = T, col = "black")
}
dev.off()

# recreating species richness in correct projection for indopacific
# richness rasterize function
# w = directory where all the shape files are ('./data/polygon')
# y = directory to be created ('./data/raster/sp')
# output value is a new res_stack, name it and save it accordingly
richness_rasterize <- function(w, y) {
  load('./data/raster/oceans_indo.Rdata')
  sp_poly_files <- dir(w)
  sp_raster_files <- sub("json", "grd", sp_poly_files)
  sp_raster_files <- sub(" ", "_", sp_raster_files)
  
  dir.create(y)
  
  cl <- makeCluster(24) 
  registerDoParallel(cl) 
  #clusterExport(cl, list('sp_poly_files', 'sp_raster_files'))
  
  foreach(i = seq_along(sp_poly_files),
          .packages = c("raster", "rgdal")) %dopar% {
            # read in 
            temp_poly = readOGR(dsn = paste0(w, "/", sp_poly_files[i]),
                                layer = "OGRGeoJSON")
            # convert projection to cea
            temp_poly = spTransform(temp_poly, CRS("+proj=cea +lat_0=0 +lon_0=125 +units=km"))
            # add field that will provide values when rasterized
            temp_poly@data$occur = 1
            # rasterize
            sp_raster = rasterize(temp_poly, oceans_indo, field = 'occur')
            writeRaster(sp_raster, 
                        filename = paste0(y, "/", sp_raster_files[i]),
                        datatype = "LOG1S", overwrite = TRUE)
          }
  stopCluster(cl)
  
  sp_stack = stack(sapply(sp_raster_files, function(x) 
    paste0(y, "/", x)))
  names(sp_stack) = sub('.grd', ' ', names(sp_stack))
  
  # create a list of stack at each resultion
  factor_val <- c(2, 4, 8, 16, 32)
  res_stack <- lapply(factor_val, function(x) 
    aggregate(sp_stack, fac = x, fun = sum) > 0)
  res_stack <- c(sp_stack, res_stack)
  return(res_stack)
}

for_indo_stack <- richness_rasterize('./data/polygon', './data/raster/sp')

# richness_plot function  
# creating a species richness layer for each resolution
# y = res_stack created in above function (sp_res_stack)
# w = file path for pdf ('./figures/species_richness_maps.pdf')
# output is the richness list, name and save accordingly
richness_plot <- function(y, w) {
  richness_list = lapply(y, function(x)
    calc(x, fun = sum, na.rm = T))
  
  pdf(w)
  for (i in 1:6) {
    test <- rasterize(continents_indo, richness_list[[i]], getCover = T)
    is.na(values(richness_list[[i]])) <- values(test) > 90
    plot(richness_list[[i]], 
         main=paste('resolution =', res(richness_list[[i]])))
    plot(continents_indo, add = T, col = "black")
  }
  dev.off()
  return(richness_list)
}

for_indo <- richness_plot(for_indo_stack, './figures/for_indo.pdf')

indopacific_list <- realm_crop(for_indo, central_indopacific)
save(indopacific_list, file = './data/raster/indopacific_list.Rdata')
load('./data/raster/indopacific_list.Rdata')

pdf('./figures/indopacific_richness.pdf')
for (i in 1:6) {
  plot(indopacific_list[[i]])
  #plot(continents_indo, add = T, col = "black")
}
dev.off()

indopacific_stack <- realm_crop(for_indo_stack, central_indopacific)
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

trop_atlantic_tree <- mini_tree('./data/tropical_atlantic')
indopacific_tree <- mini_tree('./data/indopacific')

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
  indopacific_mrd[[i]] <- mrd_raster
}
dev.off()


save(indopacific_mrd, file = './data/raster/indopacific_mrd.Rdata')
load('./data/raster/indopacific_mrd.Rdata')

# psv for regions
psv_test_list <- vector("list", length = 6)
  for (i in 1:6) {
    psv_test <- psv(trop_atlantic_mat[[i]], trop_atlantic_tree)
    psv_test_list[[i]] <- psv_test
}
  
load('./data/raster/trop_atlantic_list.Rdata')
pdf('./figures/trop_atlantic_psv.pdf')
trop_atlantic_psv <- vector("list", length = 6)
for (i in 1:6) {
    psv_raster <- trop_atlantic_list[[i]]
    psv_raster@data@values <- psv_test_list[[i]]$PSVs
    plot(psv_raster)
    plot(continents, add = T, col = "black")
    trop_atlantic_psv[[i]] <- psv_raster
}
  dev.off()
  
save(trop_atlantic_psv, file = './data/raster/trop_atlantic_psv.Rdata')  
load('./data/raster/trop_atlantic_psv.Rdata')


psv_test_list <- vector("list", length = 6)
for (i in 1:6) {
  psv_test <- psv(indopacific_mat[[i]], indopacific_tree)
  psv_test_list[[i]] <- psv_test
}

load('./data/raster/indopacific_list.Rdata')
pdf('./figures/indopacific_psv.pdf')
indopacific_psv <- vector("list", length = 6)
for (i in 1:6) {
  psv_raster <- indopacific_list[[i]]
  psv_raster@data@values <- psv_test_list[[i]]$PSVs
  plot(psv_raster)
  indopacific_psv[[i]] <- psv_raster
}
dev.off()

save(indopacific_psv, file = './data/raster/indopacific_psv.Rdata')  
load('./data/raster/indopacific_psv.Rdata')

# beta for regions
is.rooted(trop_atlantic_tree)
multi2di(trop_atlantic_tree)
maxlik.betasplit(trop_atlantic_tree, confidence.interval = "profile")
plot(trop_atlantic_tree)
identify.phylo(trop_atlantic_tree)
trop_atlantic_tree <- root(trop_atlantic_tree, outgroup = 128, resolve.root = T)

is.rooted(indopacific_tree)
multi2di(indopacific_tree)
maxlik.betasplit(indopacific_tree, confidence.interval = "profile")


# energy gradient for regions
load('./data/raster/temp_list.Rdata')

trop_atlantic_temp <- realm_crop(temp_list, tropical_atlantic)
indopacific_temp <- realm_crop(temp_list, central_indopacific)

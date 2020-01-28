rm(list = ls())

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

source('./scripts/all_functions.R')

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
# species_raster
# wanted_realm = realm spatial object
realm_crop <- function(species_raster, wanted_realm) {
  temporary <- vector("list", length = 6)
  for (i in 1:6) {
    ras <- intersect(species_raster[[i]], wanted_realm)
    temporary[[i]] <- ras
  }
  return(temporary)
}

load('./data/raster/species_richness.Rdata')
load('./data/raster/sp_res_stack.Rdata')
load('./data/continent/continent_new.Rdata')
continents_indo <- spTransform(continents, 
                               CRS("+proj=cea +lat_0=0 +lon_0=125 +units=km"))
oceans <- readOGR(dsn = "./data/Environment", layer = "ne_10m_ocean")
oceans_indo <- spTransform(oceans, CRS("+proj=cea +lat_0=0 +lon_0=125 +units=km"))
oceans_indo <- raster(oceans_indo)
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
  test <- rasterize(continents, trop_atlantic_list[[i]], getCover = T)
  trop_atlantic_list[[i]][values(test) > 0.9] <- NA
  plot(trop_atlantic_list[[i]])
  plot(continents, add = T, col = "black")
}
dev.off()

# recreating species richness in correct projection for indopacific
for_indo_stack <- richness_rasterize('./data/polygon', './data/raster/sp', oceans_indo, 
                                     "+proj=cea +lat_0=0 +lon_0=125 +units=km")

for_indo <- richness_plot(for_indo_stack, './figures/for_indo.pdf', continents_indo)

indopacific_list <- realm_crop(for_indo, central_indopacific)
save(indopacific_list, file = './data/raster/indopacific_list.Rdata')
load('./data/raster/indopacific_list.Rdata')

pdf('./figures/indopacific_richness.pdf')
for (i in 1:6) {
  test <- rasterize(continents_indo, indopacific_list[[i]], getCover = T)
  indopacific_list[[i]][values(test) > 0.9] <- NA
  plot(indopacific_list[[i]])
  plot(continents_indo, add = T, col = "black")
}
dev.off()

indopacific_stack <- realm_crop(for_indo_stack, central_indopacific)
save(indopacific_stack, file = './data/raster/indopacific_stack.Rdata')
load('./data/raster/indopacific_stack.Rdata')

# creating dir for phylogenetic functions
# Pull_names function to pull out names of species in the realm
# res_stack = res_stack of raster list wanted
pull_names <- function(res_stack) {
names_list <- vector("list", length = 6)
s <- vector("list", length = 6)
for (i in 1:6) {
  for (j in 1:534) {
    #print(any(trop_atlantic_stack[[i]][[j]]@data@values==1))
    if (sum(res_stack[[i]][[j]]@data@values==1,na.rm=T)>0)
      {
      names_list[[i]][[j]] <- res_stack[[i]][[j]]@data@names
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
load('./data/shark_tree_clean.Rdata')

trop_atlantic_tree <- mini_tree('./data/tropical_atlantic')
indopacific_tree <- mini_tree('./data/indopacific')

trop_atlantic_mat <- com_mat('./data/tropical_atlantic', trop_atlantic_tree, 
                             trop_atlantic_stack)
indopacific_mat <- com_mat('./data/indopacific', indopacific_tree, 
                           indopacific_stack)

trop_atlantic_mrd <- mrd_runplot(trop_atlantic_tree, trop_atlantic_mat, continents, 
                                 './figures/trop_atlantic_mrd.pdf', trop_atlantic_list)
indopacific_mrd <- mrd_runplot(indopacific_tree, indopacific_mat, continents_indo, 
                               './figures/indopacific_mrd.pdf', indopacific_list)

# dividing mrd by its max to standardize units
for (i in 1:6) {
  max_wanted <- max(trop_atlantic_mrd[[i]]@data@values, na.rm = T)
  new_ras <- calc(trop_atlantic_mrd[[i]], fun = function(x) x/max_wanted)
  trop_atlantic_mrd[[i]] <- new_ras
}
save(trop_atlantic_mrd, file = './data/raster/trop_atlantic_mrd.Rdata')

for (i in 1:6) {
  max_wanted <- max(indopacific_mrd[[i]]@data@values, na.rm = T)
  new_ras <- calc(indopacific_mrd[[i]], fun = function(x) x/max_wanted)
  indopacific_mrd[[i]] <- new_ras
}
save(indopacific_mrd, file = './data/raster/indopacific_mrd.Rdata')

# all indo plots still broken
pdf('./figures/trop_atlantic_mrd.pdf')
for (i in 1:6) {
  plot(trop_atlantic_mrd[[i]])
  plot(continents, add = T, col = "black")
}
dev.off()

pdf('./figures/indopacific_mrd.pdf')
for (i in 1:6) {
  plot(indopacific_mrd[[i]])
  plot(continents_indo, add = T, col = "black")
}
dev.off()

trop_atlantic_psv <- psv_runplot(trop_atlantic_tree, trop_atlantic_mat, continents, 
                                 './figures/trop_atlantic_psv.pdf', trop_atlantic_list)
indopacific_psv <- psv_runplot(indopacific_tree, indopacific_mat, continents_indo, 
                               './figures/indopacific_psv.pdf', indopacific_list)
  
save(trop_atlantic_psv, file = './data/raster/trop_atlantic_psv.Rdata')  
load('./data/raster/trop_atlantic_psv.Rdata')
save(indopacific_psv, file = './data/raster/indopacific_psv.Rdata')  
load('./data/raster/indopacific_psv.Rdata')

# beta for regions
is.rooted(trop_atlantic_tree)
trop_atlantic_tree <- multi2di(trop_atlantic_tree)
maxlik.betasplit(trop_atlantic_tree, confidence.interval = "profile")

# beta for trop atlantic = -0.97, conf interval -0.61 to -1.25
plot(trop_atlantic_tree)
#identify.phylo(trop_atlantic_tree)
trop_atlantic_tree <- root(trop_atlantic_tree, outgroup = 128, resolve.root = T)

is.rooted(indopacific_tree)
indopacific_tree <- multi2di(indopacific_tree)
maxlik.betasplit(indopacific_tree, confidence.interval = "profile")

# beta for indopacific = -0.79, conf interval -0.42 to -1.08

# energy gradient for regions
load('./data/raster/temp_list.Rdata')

trop_atlantic_temp <- realm_crop(temp_list, tropical_atlantic)
save(trop_atlantic_temp, file = './data/raster/trop_atlantic_temp.Rdata')
load('./data/raster/trop_atlantic_temp.Rdata')

indopacific_temp <- realm_crop(temp_list, central_indopacific)
save(indopacific_temp, file = './data/raster/indopacific_temp.Rdata')
load('./data/raster/indopacific_temp.Rdata')


# data analysis
data_list_trop <- vector("list", length = 6)
for (i in 1:6) {
  coo_trop <- coordinates(trop_atlantic_list[[i]])
  longitude_trop <- coo_trop[,1]
  latitude_trop <- coo_trop[,2]
  dat <- data.frame(longitude_trop, latitude_trop, values(trop_atlantic_list[[i]]), 
                    values(trop_atlantic_mrd[[i]]), values(trop_atlantic_psv[[i]]), 
                    values(trop_atlantic_temp[[i]]))
  colnames(dat) <- c("tropical_atlantic_longitude", "tropical_atlantic_latitude", 
                     "tropical_atlantic_richness", "tropical_atlantic_mrd", "tropical_atlantic_psv",
                     "tropical_atlantic_temperature")
  data_list_trop[[i]] <- dat
}
# disaggregating list into dataframe
total_df_trop <- dplyr::bind_rows(data_list_trop, .id = 'scale')
# remove rows when richness is zero
total_df_trop <- subset(total_df_trop, tropical_atlantic_richness != 0)
save(total_df_trop, file = "data/raster/total_df_trop.Rdata")


data_list_indo <- vector("list", length = 6)
for (i in 1:6) {
  coo_indo <- coordinates(indopacific_list[[i]])
  longitude_indo <- coo_indo[,1]
  latitude_indo <- coo_indo[,2]
  dat <- data.frame(longitude_indo, latitude_indo,
                    values(indopacific_list[[i]]), values(indopacific_mrd[[i]]),
                    values(indopacific_psv[[i]]), values(indopacific_temp[[i]]))
  colnames(dat) <- c("indopacific_longitude", "indopacific_latitude",
                     "indopacific_richness", "indopacific_mrd", "indopacifc_psv", "indopacific_temperature")
  data_list_indo[[i]] <- dat
}
# disaggregating list into dataframe
total_df_indo <- dplyr::bind_rows(data_list_indo, .id = 'scale')
# remove rows when richness is zero
total_df_indo <- subset(total_df_indo, indopacific_richness != 0)
save(total_df_indo, file = "data/raster/total_df_indo.Rdata")

# regressions for tropical atlantic
trop_atlantic_richnessVtemp <- make_plot(total_df_trop, 'tropical_atlantic_temperature', 'tropical_atlantic_richness', 
                                              "Temperature (°C)", "Tropical Atlantic Richness", 
                                              './figures/trop_atlantic_vs_temp.pdf')
save(trop_atlantic_richnessVtemp, file = './data/stats/trop_atlantic_richnessVtemp.Rdata')

trop_atlantic_richnessVmrd <- make_plot(total_df_trop, 'tropical_atlantic_richness', 'tropical_atlantic_mrd',
                                             "Tropical Atlantic Richness", 
                                             "Mean Root Distance (MRD) Tropical Atlantic",
                                             './figures/trop_atlantic_mrd_vs_richness.pdf')
save(trop_atlantic_richnessVmrd, file = './data/stats/trop_atlantic_richnessVmrd.Rdata')

trop_atlantic_tempVmrd <- make_plot(total_df_trop, 'tropical_atlantic_temperature', 'tropical_atlantic_mrd', 
                                         "Temperature (°C)", "Mean Root Distance (MRD) Tropical Atlantic",
                                         './figures/trop_atlantic_temp_vs_mrd.pdf')  
save(trop_atlantic_tempVmrd, file = './data/stats/trop_atlantic_tempVmrd.Rdata')

# regressions for indopacific
indopacific_richnessVtemp <- make_plot(total_df_indo, 'indopacific_temperature', 'indopacific_richness', 
                                                "Temperature (°C)", "Central Indopacific Richness",
                                                './figures/indopacific_richness_vs_temp.pdf')
save(indopacific_richnessVtemp, file = './data/stats/indopacific_richnessVtemp.Rdata')

indopacific_richnessVmrd <- make_plot(total_df_indo, 'indopacific_richness', 'indopacific_mrd',
                                           "Central Indopacific Richness", 
                                           "Mean Root Distance (MRD) Central Indopacific", 
                                           './figures/indopacific_mrd_vs_richness.pdf')
save(indopacific_richnessVmrd, file = './data/stats/indopacific_richnessVmrd.Rdata')

indopacific_tempVmrd <- make_plot(total_df_indo, 'indopacific_temperature', 'indopacific_mrd', "Temperature (°C)",
                                       "Mean Root Distance (MRD) Central Indopacific",
                                       './figures/indopacific_temp_vs_mrd.pdf')
save(indopacific_tempVmrd, file = './data/stats/indopacific_tempVmrd.Rdata')



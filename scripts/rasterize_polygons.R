library(rgdal)
library(sp)
library(maps)
library(maptools)
library(raster)

# read in the ocean
oceans <- readOGR(dsn = "./data/Environment", layer = "ne_10m_ocean")

# create a global raster layer
oceans <- spTransform(oceans, CRS("+proj=cea +units=km"))
plot(oceans)
oceans_raster <- raster(oceans)
res(oceans_raster) <- 110

# saving the world raster grid
save(oceans_raster, file = './data/raster/oceans_raster.Rdata')
load('./data/raster/oceans_raster.Rdata')

# code to change resolution of rasters
factor_val <- c(1, 2, 4, 8, 16, 32)

# making continents polygon
world <- map(database = "world", fill = T, plot = F)
continents <- map2SpatialPolygons(world, IDs = world$names, proj4string = CRS("+proj=longlat"))
continents <- spTransform(continents, CRS("+proj=cea +units=km"))

# raster species polygons
sp_files <- dir('./data/polygon')
sp_raster <- vector("list", length = length(sp_files))
raster_res_list <- vector("list", length = length(factor_val))
for (j in seq_along(factor_val)) {
     for (i in seq_along(sp_files)) {
          # read in 
          temp_poly = readOGR(dsn = paste("./data/polygon/", sp_files[i], 
                                          sep=''), layer = "OGRGeoJSON")
          # convert projection to cea
          temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
          # add field that will provide values when rasterized
          temp_poly@data$occur = 1
          # rasterize
          sp_raster[[i]] = rasterize(temp_poly, oceans_raster, field = 'occur')
     }
  raster_res_list[[j]] = sp_raster
  raster_res_list[[j]] = aggregate(sp_raster[[i]], fac = j, fun = sum) > 0
}

# now aggregate this finest spatial grid to the coarser resolutions

raster_res_list[[j]] = aggregate(sp_raster[[i]], fac=16, fun = sum) > 0)


# Note use rasterize( ... , getCover = TRUE) to estimate how much
# of a pixel the land is covering or vice-versa how much 
# ocean is covering. 

  
save(raster_res_list, file = './data/raster/raster_res_list.Rdata')
load('./data/raster/raster_res_list.Rdata')

# creating a species richness layer for each resolution
pdf('./figures/species_richness_maps_unmasked.pdf')
species_richness_list <- vector("list", length = length(raster_res_list))
for (i in seq_along(raster_res_list)) {
     sp_raster_stack <- stack(raster_res_list[[i]])
     species_richness <- calc(sp_raster_stack, fun = sum, na.rm = T)
     plot(species_richness, main=paste('resolution =', res(res_list[[i]])))
     plot(continents, add = T, col = "black")
     species_richness_list[[i]] <- species_richness
}
dev.off()

save(species_richness_list, file = './data/raster/species_richness.Rdata')
load('./data/raster/species_richness.Rdata')

# code to find the position of the max value
indx <- which.max(species_richness)
pos <- xyFromCell(species_richness, indx)
pos

# create a temperature raster
temp <- read.csv('./data/Environment/temp.csv')
head(temp)
temp$Meandepth <- rowMeans(temp[,3:87], na.rm = TRUE)
names(temp)
class(temp)
coordinates(temp) <- ~LONGITUDE + LATITUDE
class(temp)
proj4string(temp) <- "+proj=longlat +datum=WGS84"
temp <- spTransform(temp, CRS("+proj=cea +units=km"))
pdf('./figures/temperature_unmasked.pdf')
temp_list <- vector("list", length = length(factor_val))
for (i in seq_along(factor_val)) {
     temp_raster <- rasterize(temp, oceans_raster, 'Meandepth')
     temp_raster = aggregate(temp_raster, fac = i, fun = sum) > 0
     plot(temp_raster, main = paste('resolution =', res(res_list[[i]])))
     plot(continents, add = T, col = "black")
     temp_list[[i]] <- temp_raster
}
dev.off()

save(temp_list, file = './data/raster/temp_list.Rdata')
load('./data/raster/temp_list.Rdata')

# chlorophyll
chloro <- raster('./data/Environment/MY1DMM_CHLORA_2017-06-01_rgb_360x180.TIFF')
chloro <- rasterToPolygons(chloro)
chloro <- spTransform(chloro, CRS("+proj=cea +units=km"))
pdf('./figures/chlorophyll.pdf')
chloro_list <- vector("list", length = length(factor_val))
for (i in seq_along(factor_val)) {
     chloro_ras <- rasterize(chloro, oceans_raster, 
                             'MY1DMM_CHLORA_2017.06.01_rgb_360x180')
     chloro_ras <- aggregate(chloro_ras, fac = i, fun = sum) > 0
     values(chloro_ras)[values(chloro_ras) == 255] <- NA
     plot(chloro_ras, main = paste('resolution =', res(res_list[[i]])))
     plot(continents, add = T, col = "black")
     chloro_list[[i]] <- chloro_ras
}
dev.off()

save(chloro_list, file = './data/raster/chloro_list.Rdata')
load('./data/raster/chloro_list.Rdata')

# IUCN shark species richness
sp_files_IUCN <- dir('./data/IUCN')
IUCN_raster <- vector("list", length = length(sp_files_IUCN))
IUCN_res_list <- vector("list", length = length(factor_val))
for (j in seq_along(factor_val)) {
  for (i in seq_along(sp_files_IUCN)) {
    # read in 
    temp_poly = readOGR(dsn = paste("./data/polygon/", sp_files_IUCN[i], 
                                    sep=''), layer = "OGRGeoJSON")
    # convert projection to cea
    temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
    # add field that will provide values when rasterized
    temp_poly@data$occur = 1
    # rasterize
    IUCN_raster[[i]] = rasterize(temp_poly, oceans_raster, field = 'occur')
  }
  IUCN_res_list[[j]] = IUCN_raster
  IUCN_res_list[[j]] = aggregate(IUCN_raster[[i]], fac= j, fun = sum) > 0
}

save(IUCN_res_list, file = './data/raster/IUCN_res_list.Rdata')
load('./data/raster/IUCN_res_list.Rdata')

# creating an IUCN richness layer for each resolution
pdf('./figures/IUCN_richness_maps_unmasked.pdf')
IUCN_richness_list <- vector("list", length = length(IUCN_res_list))
for (i in seq_along(IUCN_res_list)) {
     sp_raster_stack <- stack(IUCN_res_list[[i]])
     IUCN_richness <- calc(sp_raster_stack, fun = sum, na.rm = T)
     plot(IUCN_richness, main = paste('resolution =', res(res_list[[i]])))
     plot(continents, add = T, col = "black")
     IUCN_richness_list[[i]] <- IUCN_richness
}
dev.off()

save(IUCN_richness_list, file = './data/raster/IUCN_richness_list.Rdata')
load('./data/raster/IUCN_richness_list.Rdata')

# making a value for latitude
latitude_list <- vector("list", length = length(res_list))
for (i in seq_along(res_list)) {
     oceans_p <- rasterToPoints(res_list[[i]])
     oceans_p_df <- data.frame(oceans_p)
     latitude <- oceans_p_df$y
     latitude_list[[i]] <- latitude
}

# salinity
salinity <- read.csv('./data/Environment/salinity.csv')
salinity$Meandepth <- rowMeans(salinity[,3:86], na.rm = TRUE)
coordinates(salinity) <- ~ Longitude + Latitude
proj4string(salinity) <- "+proj=longlat +datum=WGS84"
salinity <- spTransform(salinity, CRS("+proj=cea +units=km"))
pdf('./figures/salinity_unmasked.pdf')
salinity_list <- vector("list", length = length(factor_val))
for (i in seq_along(factor_val)) {
     salinity_raster <- rasterize(salinity, oceans_raster, 'Meandepth')
     salinity_raster <- aggregate(salinity_raster, fac = i, fun = sum) > 0
     plot(salinity_raster, main = paste('resolution =', res(res_list[[i]])))
     plot(continents, add = T, col = "black")
     salinity_list[[i]] <- salinity_raster
}
dev.off()

save(salinity_list, file = './data/raster/salinity_list')
load('./data/raster/salinity_list')

# distance from coast
coast_distance_list <- vector("list", length = length(res_list))
pdf('./figures/distance_from_coast_unmasked.pdf')
for (i in seq_along(res_list)) {
  distance_raster <- setValues(res_list[[i]], 0)
  distance_raster <- mask(distance_raster, mask_ras_list, inverse = T)
  rd <- distance(distance_raster)
  plot(rd, main = paste('resolution =', res(res_list[[i]])))
  plot(continents, add = T, col = "black")
  coast_distance_list[[i]] <- rd
}
dev.off()



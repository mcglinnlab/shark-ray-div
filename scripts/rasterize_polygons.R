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

# for loop to change resolution of oceans_raster
factor_val <- c(110, 220, 330, 440, 550, 770, 1100)
res_list <- vector("list", length = length(factor_val))
res_list <- lapply(res_list, function(x) oceans_raster)
for (i in seq_along(factor_val)) {
     # making new resolutions
     res(res_list[[i]]) <- factor_val[i]
}

# for loop to make a mask raster list
mask_ras_list <- vector("list", length = length(res_list))
for (i in seq_along(res_list)) {
     mask_ras <- rasterize(oceans, res_list[[i]], field = 1)
     mask_ras_list[[i]] <- mask_ras
}

# raster species polygons
sp_files <- dir('./data/polygon')
sp_raster <- vector("list", length = length(sp_files))
raster_res_list <- vector("list", length = length(res_list))
for (j in seq_along(res_list)) {
     for (i in seq_along(sp_files)) {
          # read in 
          temp_poly = readOGR(dsn = paste("./data/polygon/", sp_files[i], 
                                          sep=''), layer = "OGRGeoJSON")
          # convert projection to cea
          temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
          # add field that will provide values when rasterized
          temp_poly@data$occur = 1
          # rasterize
          sp_raster[[i]] = rasterize(temp_poly, res_list[[j]], field = 'occur')
     }
     raster_res_list[[j]] = sp_raster
}
  
save(raster_res_list, file = './data/raster/raster_res_list.Rdata')
load('./data/raster/raster_res_list.Rdata')

# creating a species richness layer for each resolution
species_richness_list <- vector("list", length = length(raster_res_list))
for (i in seq_along(raster_res_list)) {
     sp_raster_stack <- stack(raster_res_list[[i]])
     species_richness <- calc(sp_raster_stack, fun = sum, na.rm = T)
     species_richness <- mask(species_richness, mask_ras_list[[i]])
     plot(species_richness)
     plot(oceans, add = T)
     species_richness_list[[i]] <- species_richness
}

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
temp_list <- vector("list", length = length(res_list))
for (i in seq_along(res_list)) {
     temp_raster <- rasterize(temp, res_list[[i]], 'Meandepth')
     temp_raster <- mask(temp_raster, mask_ras_list[[i]])
     plot(temp_raster)
     plot(oceans, add = T)
     temp_list[[i]] <- temp_raster
}

save(temp_list, file = './data/raster/temp_list.Rdata')
load('./data/raster/temp_list.Rdata')

# chlorophyll
chloro <- raster('./data/Environment/MY1DMM_CHLORA_2017-06-01_rgb_360x180.TIFF')
chloro <- rasterToPolygons(chloro)
chloro <- spTransform(chloro, CRS("+proj=cea +units=km"))
chloro_list <- vector("list", length = length(res_list))
for (i in seq_along(res_list)) {
     chloro_ras <- rasterize(chloro, res_list[[i]], 
                             'MY1DMM_CHLORA_2017.06.01_rgb_360x180')
     values(chloro_ras)[values(chloro_ras) == 255] <- NA
     plot(chloro_ras)
     plot(oceans, add = T)
     chloro_list[[i]] <- chloro_ras
}

save(chloro_list, file = './data/raster/chloro_list.Rdata')
load('./data/raster/chloro_list.Rdata')

# IUCN shark species richness
sp_files_IUCN <- dir('./data/IUCN')
IUCN_raster <- vector("list", length = length(sp_files_IUCN))
IUCN_res_list <- vector("list", length = length(res_list))
for (j in seq_along(res_list)) {
  for (i in seq_along(sp_files_IUCN)) {
    # read in 
    temp_poly = readOGR(dsn = paste("./data/polygon/", sp_files_IUCN[i], 
                                    sep=''), layer = "OGRGeoJSON")
    # convert projection to cea
    temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
    # add field that will provide values when rasterized
    temp_poly@data$occur = 1
    # rasterize
    IUCN_raster[[i]] = rasterize(temp_poly, res_list[[j]], field = 'occur')
  }
  IUCN_res_list[[j]] = IUCN_raster
}

save(IUCN_res_list, file = './data/raster/IUCN_res_list.Rdata')
load('./data/raster/IUCN_res_list.Rdata')

# creating an IUCN richness layer for each resolution
IUCN_richness_list <- vector("list", length = length(IUCN_res_list))
for (i in seq_along(IUCN_res_list)) {
  sp_raster_stack <- stack(IUCN_res_list[[i]])
  IUCN_richness <- calc(sp_raster_stack, fun = sum, na.rm = T)
  IUCN_richness <- mask(IUCN_richness, mask_ras_list[[i]])
  plot(IUCN_richness)
  plot(oceans, add = T)
  IUCN_richness_list[[i]] <- IUCN_richness
}

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
salinity_list <- vector("list", length = length(res_list))
for (i in seq_along(res_list)) {
     salinity_raster <- rasterize(salinity, res_list[[i]], 'Meandepth')
     salinity_raster <- mask(salinity_raster, mask_ras_list[[i]])
     plot(salinity_raster)
     plot(oceans, add = T)
     salinity_list[[i]] <- salinity_raster
}

save(salinity_list, file = './data/raster/salinity_list')

library(rgdal)
library(sp)
library(raster)
library(doParallel) 
library(foreach)

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

# making continents polygon
continents <- shapefile('./data/polygon/continent/continent.shp')
continents <- spTransform(continents, CRS("+proj=cea +units=km"))

# rasterize species polygons to the 110 scale
# then write them to file
sp_poly_files <- dir("./data/polygon")
sp_raster_files <- sub("json", "grd", sp_poly_files)
sp_raster_files <- sub(" ", "_", sp_raster_files)

dir.create("./data/raster/sp")

cl <- makeCluster(24) 
registerDoParallel(cl) 
#clusterExport(cl, list('sp_poly_files', 'sp_raster_files'))

foreach(i = seq_along(sp_poly_files),
        .packages = c("raster", "rgdal")) %dopar% {
     # read in 
     temp_poly = readOGR(dsn = paste0("./data/polygon/", sp_poly_files[i]),
                                layer = "OGRGeoJSON")
     # convert projection to cea
     temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
     # add field that will provide values when rasterized
     temp_poly@data$occur = 1
     # rasterize
     sp_raster = rasterize(temp_poly, oceans_raster, field = 'occur')
     writeRaster(sp_raster, 
                 filename = paste0("./data/raster/sp/", sp_raster_files[i]),
                 datatype = "LOG1S", overwrite = TRUE)
}
stopCluster(cl)

# stack rasterized files and aggregate to other scales

sp_stack = stack(sapply(sp_raster_files, function(x) 
                 paste0("./data/raster/sp/", x)))
names(sp_stack) = sub('.grd', ' ', names(sp_stack))

# create a list of stack at each resultion
factor_val <- c(1, 2, 4, 8, 16, 32)
sp_res_stack = lapply(factor_val, function(x) 
                      aggregate(sp_stack, fac = x, fun = sum) > 0)

save(sp_res_stack, file = './data/raster/sp_res_stack.Rdata')
load('./data/raster/sp_res_stack.Rdata')

# creating a species richness layer for each resolution
species_richness = lapply(sp_res_stack, function(x)
                  calc(x, fun = sum, na.rm = T))

pdf('./figures/species_richness_maps.pdf')
for (i in seq_along(species_richness)) {
     plot(species_richness[[i]], 
          main=paste('resolution =', res(species_richness[[i]])))
     plot(continents, add = T, col = "black")
}
dev.off()


save(species_richness, file = './data/raster/species_richness.Rdata')
load('./data/raster/species_richness.Rdata')

# code to find the position of the max value
indx <- which.max(species_richness[[1]])
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


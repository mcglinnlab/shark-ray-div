library(rgdal)
library(sp)
library(maps)
library(maptools)
library(raster)

# create a global raster layer
world <- map(database = "world", fill =T)
world <- map2SpatialPolygons(world, IDs = world$names, CRS("+proj=longlat +datum=WGS84"))
world <- spTransform(world, CRS("+proj=cea +units=km"))
plot(world)
world_raster <- raster(world)
res(world_raster) <- 110

# saving the world raster grid
save(world_raster, file = './data/raster/world_raster.Rdata')

# for loop to change resolution
factor_val <- c(110, 220, 330, 440, 550, 770, 1100)
res_list <- vector("list", length = length(factor_val))
res_list <- lapply(res_list, function(x) world_raster)
for (i in seq_along(factor_val)) {
  # making new resolutions
  res(res_list[[i]]) <- factor_val[i]
}

# raster species polygons
sp_files = dir('./data/polygon')
sp_raster <- vector("list", length = length(sp_files))
raster_res_list <- vector("list", length = length(res_list))
for (j in seq_along(res_list)) {
     for (i in seq_along(sp_files)) {
    # read in 
    temp_poly = readOGR(dsn = paste("./data/polygon/", sp_files[i], sep=''),
                        layer = "OGRGeoJSON")
    # convert projection to cea
    temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
    # add field that will provide values when rasterized
    temp_poly@data$occur = 1
    # rasterize
    sp_raster[[i]] = rasterize(temp_poly, res_list[[j]], field = 'occur')
     }
  raster_res_list[j] = lapply(raster_res_list, function(x) sp_raster)
}
  
save(raster_res_list, file = './data/raster/raster_res_list.Rdata')
load('./data/raster/raster_res_list.Rdata')

# creating a species richness layer for each resolution
species_richness_list <- vector("list", length = length(raster_res_list))
for (i in seq_along(raster_res_list)) {
sp_raster_stack <- stack(raster_res_list[[i]])
species_richness <- calc(sp_raster_stack, fun = sum, na.rm = T)
plot(species_richness)
plot(world, add = T)
species_richness_list[i] <- lapply(species_richness_list, function(x) species_richness)
}

save(species_richness_list, file = './data/raster/species_richness.Rdata')

# code to find the position of the max value
indx <- which.max(species_richness)
pos <- xyFromCell(species_richness, indx)
pos




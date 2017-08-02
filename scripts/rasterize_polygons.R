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

# raster species polygons
sp_files = dir('./data/polygon')
sp_raster <- vector("list", length = length(sp_files))
for (i in seq_along(sp_files)) {
    # read in 
    temp_poly = readOGR(dsn = paste("./data/polygon/", sp_files[i], sep=''),
                        layer = "OGRGeoJSON")
    # convert projection to cea
    temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
    # add field that will provide values when rasterized
    temp_poly@data$occur = 1
    # rasterize
    sp_raster[[i]] = rasterize(temp_poly, world_raster, field = 'occur')
}
  
load('./data/raster/sp_raster.Rdata')
sp_raster_stack <- stack(sp_raster)
species_richness <- calc(sp_raster_stack, fun = sum, na.rm = T)
plot(species_richness, add = T)

save(species_richness, file = './data/raster/species_richness.Rdata')

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
trop_atlantic_list <- vector("list", length = 6)
for (i in 1:6) {
     ras <- intersect(species_richness[[i]], tropical_atlantic)
     trop_atlantic_list[[i]] <- ras
}
save(trop_atlantic_list, file = './data/raster/trop_atlantic_list.Rdata')
load('./data/raster/trop_atlantic_list.Rdata')

indopacific_list <- vector("list", length = 6)
for (i in 1:6) {
  ras <- intersect(species_richness[[i]], central_indopacific)
  indopacific_list[[i]] <- ras
}
save(indopacific_list, file = './data/raster/indopacific_list.Rdata')
load('./data/raster/indopacific_list.Rdata')
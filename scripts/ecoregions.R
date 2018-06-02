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
trop_atlantic_names <- vector("list", length = 6)
for (i in 1:6) {
  for (j in 1:534) {
    #print(any(trop_atlantic_stack[[i]][[j]]@data@values==1))
    if (sum(trop_atlantic_stack[[i]][[j]]@data@values==1,na.rm=T)>0)
      {
      trop_atlantic_names[[i]][[j]] <- trop_atlantic_stack[[i]][[j]]@data@names
      } else {trop_atlantic_names[[i]][[j]] <- NA}
  }
}



# phylogenetic metrics for realms
source('./scripts/Phylogenetic_metrics.R')



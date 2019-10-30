# clean and set up data code
library(rgdal)
library(sp)
library(raster)
library(ncdf4)
library(R.utils)

# read in the ocean
oceans <- readOGR(dsn = "./data/Environment", layer = "ne_10m_ocean")

# create a global raster layer
oceans <- spTransform(oceans, CRS("+proj=cea +units=km"))
oceans_raster <- raster(oceans)
res(oceans_raster) <- 110

# saving the world raster grid
save(oceans_raster, file = './data/raster/oceans_raster.Rdata')
load('./data/raster/oceans_raster.Rdata')

# making continents polygon
continents <- shapefile('./data/continent/continent')
continents <- spTransform(continents, CRS("+proj=cea +units=km"))
save(continents, file = './data/continent/continent_new')

# Enivironmental variables
# create a temperature raster
temp <- read.csv('./data/Environment/temp.csv')
temp$Meandepth <- rowMeans(temp[,3:87], na.rm = TRUE)
coordinates(temp) <- ~LONGITUDE + LATITUDE
proj4string(temp) <- "+proj=longlat +datum=WGS84"
temp <- spTransform(temp, CRS("+proj=cea +units=km"))
temp_raster <- rasterize(temp, oceans_raster, 'Meandepth')
save(temp_raster, file = './data/raster/temp_raster.Rdata')
load('./data/raster/temp_raster.Rdata')

# create a chlorophyll raster
chloro <- raster('./data/Environment/MY1DMM_CHLORA_2017-06-01_rgb_360x180.TIFF')
chloro_raster <- projectRaster(chloro, oceans_raster)
chloro_raster <- mask(chloro_raster, continents, inverse = T)
values(chloro_raster)[values(chloro_raster) == 255] <- NA
save(chloro_raster, file = './data/raster/chloro_raster.Rdata')
load('./data/raster/chloro_raster.Rdata')

# create a salinity raster
salinity <- read.csv('./data/Environment/salinity.csv')
salinity$Meandepth <- rowMeans(salinity[,3:86], na.rm = TRUE)
coordinates(salinity) <- ~ Longitude + Latitude
proj4string(salinity) <- "+proj=longlat +datum=WGS84"
salinity <- spTransform(salinity, CRS("+proj=cea +units=km"))
salinity_raster <- rasterize(salinity, oceans_raster, 'Meandepth')
save(salinity_raster, file = './data/raster/salinity_raster.Rdata')
load('./data/raster/salinity_raster.Rdata')

# create a bathymetry raster
bathy <- raster('./data/Environment/ETOPO1_Ice_c_gmt4.grd')
proj4string(bathy) = CRS("+proj=longlat +datum=WGS84")
bathy_raster <- projectRaster(bathy, oceans_raster)
values(bathy_raster)[values(bathy_raster) > 0] <- NA
save(bathy_raster, file = './data/raster/bathy_raster.Rdata')
load('./data/raster/bathy_raster.Rdata')

# create an area raster
area_raster <- rasterize(oceans, oceans_raster, getCover = T)
save(area_raster, file = './data/raster/area_raster.Rdata')
load('./data/raster/area_raster.Rdata')

# create a distance from the coast raster
rd <- setValues(oceans_raster, 0)
rd <- mask(rd, continents)
distance_raster <- distance(rd)
save(distance_raster, file = './data/raster/distance_raster.Rdata')
load('./data/raster/distance_raster.Rdata') 

# create a latitude raster
oceans_df <- coordinates(area_raster)
oceans_df <- data.frame(oceans_df)
new_field <- abs(oceans_df$y)
oceans_df$new_field <- new_field
coordinates(oceans_df) <- ~ x+y
proj4string(oceans_df) <- "+proj=cea +units=km"
latitude_raster <- rasterize(oceans_df, oceans_raster, 'new_field')
save(latitude_raster, file = './data/raster/latitude_raster.Rdata')
load('./data/raster/latitude_raster.Rdata')

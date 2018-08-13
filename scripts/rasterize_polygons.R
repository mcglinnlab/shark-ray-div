library(rgdal)
library(sp)
library(raster)
library(doParallel) 
library(foreach)
library(maps)
library(maptools)


# richness rasterize function
# w = directory where all the shape files are ('./data/polygon')
# y = directory to be created ('./data/raster/sp')
# output value is a new res_stack, name it and save it accordingly
richness_rasterize <- function(w, y) {
  load('./data/raster/oceans_raster.Rdata')
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
            temp_poly = spTransform(temp_poly, CRS("+proj=cea +units=km"))
            # add field that will provide values when rasterized
            temp_poly@data$occur = 1
            # rasterize
            sp_raster = rasterize(temp_poly, oceans_raster, field = 'occur')
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
    test <- rasterize(continents, richness_list[[i]], getCover = T)
    is.na(values(richness_list[[i]])) <- values(test) > 90
    plot(richness_list[[i]], 
         main=paste('resolution =', res(richness_list[[i]])))
    plot(continents, add = T, col = "black")
  }
  dev.off()
  return(richness_list)
}

# Taxonomic richness rasterize and plot
sp_res_stack <- richness_rasterize('./data/polygon', './data/raster/sp')
save(sp_res_stack, file = './data/raster/sp_res_stack.Rdata')
species_richness <- richness_plot(sp_res_stack, './figures/species_richness_maps.pdf')
save(species_richness, file = './data/raster/species_richness.Rdata')
load('./data/raster/species_richness.Rdata')

# IUCN richness rasterize and plot
iucn_res_stack <- richness_rasterize('./data/IUCN', './data/raster/iucn')
save(iucn_res_stack, file = './data/raster/iucn_res_stack.Rdata')
iucn_richness <- richness_plot(iucn_res_stack, './figures/IUCN_maps.pdf')
save(iucn_richness, file = './data/raster/iucn_richness.Rdata')
load('./data/raster/iucn_richness.Rdata')

# Carcharhiniformes rasterize and plot
car_res_stack <- richness_rasterize('./data/Carcharhiniformes', './data/raster/car')
save(car_res_stack, file = './data/raster/car_res_stack.Rdata')
car_richness <- richness_plot(car_res_stack, './figures/Carcharhiniforme_richness.pdf')
save(car_richness, file = './data/raster/car_richness.Rdata')
load('./data/raster/car_richness.Rdata')

# Lamniformes rasterize and plot
lam_res_stack <- richness_rasterize('./data/Lamniformes', './data/raster/lam')
save(lam_res_stack, file = './data/raster/lam_res_stack.Rdata')
lam_richness <- richness_plot(lam_res_stack, './figures/Lamniforme_richness.pdf')
save(lam_richness, file = './data/raster/lam_richness.Rdata')
load('./data/raster/lam_richness.Rdata')

# enviro_plot function to aggregate and plot environmental variables
# y = environmental raster from initial_cleanup
# w = file path for pdf
# output is raster list, name and save accordingly
enviro_plot <- function(y, w) {
  factor_val <- c(2, 4, 8, 16, 32)
  enviro_list <- lapply(factor_val, function (x)
  aggregate(y, fac = x, fun = mean))
  enviro_list <- c(y, enviro_list)
  return(enviro_list)
  pdf(w)
for (i in 1:6) {
  plot(enviro_list[[i]], main = paste('resolution =', res(enviro_list[[i]])))
  plot(continents, add = T, col = "black")
}
dev.off()
}

# temperature plot
temp_list <- enviro_plot(temp_raster, './figures/temperature.pdf')
save(temp_list, file = './data/raster/temp_list.Rdata')
load('./data/raster/temp_list.Rdata')

# chlorophyll plot
chloro_list <- enviro_plot(chloro_raster, './figures/chlorophyll.pdf')
save(chloro_list, file = './data/raster/chloro_list.Rdata')
load('./data/raster/chloro_list.Rdata')

# salinity plot
salinity_list <- enviro_plot(salinity_raster, './figures/salinity.pdf')
save(salinity_list, file = './data/raster/salinity_list.Rdata')
load('./data/raster/salinity_list.Rdata')

# bathymetry plot
bathy_list <- enviro_plot(bathy_raster, './figures/bathymetry.pdf')
save(bathy_list, file = './data/raster/bathy_list.Rdata')
load('./data/raster/bathy_list.Rdata')

# area plot
area_list <- enviro_plot(area_raster, './figures/area.pdf')
save(area_list, file = './data/raster/area_list.Rdata')
load('./data/raster/area_list.Rdata')

# distance from the coast plot
distance_list <- enviro_plot(distance_raster, './figures/distance_from_coast.pdf')
save(distance_list, file = './data/raster/distance_list.Rdata')
load('./data/raster/distance_list.Rdata')

# latitude plot
latitude_list <- enviro_plot(latitude_raster, './figures/latitude.pdf')
save(latitude_list, file = './data/raster/latitude_list.Rdata')
load('./data/raster/latitude_list.Rdata')

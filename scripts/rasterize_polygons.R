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
# z = res_stack to save (sp_res_stack)
# v = character string of res_stack file path ('./data/raster/sp_res_stack.Rdata')
richness_rasterize <- function(w, y, z, v) {
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
  z <- lapply(factor_val, function(x) 
    aggregate(sp_stack, fac = x, fun = sum) > 0)
  z <- c(sp_stack, z)
  
  save(z, file = v)
}


# richness_plot function  
# creating a species richness layer for each resolution
# z = res_stack created in above function (sp_res_stack)
# y = file path for loaded res_stack ('./data/raster/sp_res_stack.Rdata')
# v = output list name (species_richness)
# w = file path for pdf ('./figures/species_richness_maps.pdf')
# t = file path for raster list
richness_plot <- function(y, v, w, t) {
  load(y)
  v = lapply(z, function(x)
    calc(x, fun = sum, na.rm = T))
  
  save(v, file = t)
  
  pdf(w)
  for (i in 1:6) {
    test <- rasterize(continents, v[[i]], getCover = T)
    is.na(values(v[[i]])) <- values(test) > 90
    plot(v[[i]], 
         main=paste('resolution =', res(v[[i]])))
    plot(continents, add = T, col = "black")
  }
  dev.off()
}

# Taxonomic richness rasterize and plot
richness_rasterize('./data/use', './data/raster/sp', sp_res_stack, './data/raster/sp_res_stack.Rdata')
richness_plot('./data/raster/sp_res_stack.Rdata', species_richness, './figures/species_richness_maps.pdf', './data/raster/species_richness.Rdata')
load('./data/raster/species_richness.Rdata')

# IUCN richness rasterize and plot
richness_rasterize('./data/IUCN', './data/raster/iucn', iucn_res_stack, './data/raster/iucn_res_stack.Rdata')
richness_plot('./data/raster/iucn_res_stack.Rdata', iucn_richness, './figures/IUCN_maps.pdf', './data/raster/iucn_richness.Rdata')
load('./data/raster/iucn_richness.Rdata')

# Carcharhiniformes rasterize and plot

# Lamniformes rasterize and plot


# enviro_plot function to aggregate and plot environmental variables
# y = raster list for each variable
# v = environmental raster from initial_cleanup
# w = file path for pdf
# t = file path for raster list
enviro_plot <- function(y, v, w, t) {
  factor_val <- c(2, 4, 8, 16, 32)
  y <- lapply(factor_val, function (x)
  aggregate(v, fac = x, fun = mean))
  y <- c(v, y)
  pdf(w)
for (i in 1:6) {
  plot(y[[i]], main = paste('resolution =', res(y[[i]])))
  plot(continents, add = T, col = "black")
}
dev.off()
save(y, file = t)
}

# temperature plot
enviro_plot(temp_list, temp_raster, './figures/temperature.pdf', './data/raster/temp_list.Rdata')

# chlorophyll plot
enviro_plot(chloro_list, chloro_raster, './figures/chlorophyll.pdf', './data/raster/chloro_list.Rdata')

# salinity plot
enviro_plot(salinity_list, salinity_raster, './figures/salinity.pdf', './data/raster/salinity_list.Rdata')

# bathymetry plot
enviro_plot(bathy_list, bathy_raster, './figures/bathymetry.pdf', './data/raster/bathy_list.Rdata')

# area plot
enviro_plot(area_list, area_raster, './figures/area.pdf', './data/raster/area_list.Rdata')

# distance from the coast plot
enviro_plot(distance_list, distance_raster, './figures/distance_from_coast.pdf', './data/raster/distance_list.Rdata')

# latitude plot
enviro_plot(latitude_list, latitude_raster, './figures/latitude.pdf', './data/raster/latitude_list.Rdata')

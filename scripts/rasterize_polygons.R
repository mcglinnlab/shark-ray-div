rm(list = ls())

library(rgdal)
library(sp)
library(raster)
library(doParallel) 
library(foreach)
library(maps)
library(maptools)

source('./scripts/all_functions.R')

# This script reads in range polygons, transforms them into rasters, stacks them,
# and sums them to create richness plots. Environmental rasters are also plotted. 
# the richness_rasterize, enviro_plot, and richness_plot
# functions are coded in the all_functions.R script

# Read in continent shapefile for mapping purposes
load('./data/continent/continent_new.Rdata')

# Taxonomic richness rasterize and plot
sp_res_stack <- richness_rasterize('./data/polygon', './data/raster/sp', oceans_raster, 
                                   "+proj=cea +units=km")
save(sp_res_stack, file = './data/raster/sp_res_stack.Rdata')
species_richness <- richness_plot(sp_res_stack, './figures/species_richness_maps.pdf', continents)
save(species_richness, file = './data/raster/species_richness.Rdata')
load('./data/raster/species_richness.Rdata')

# IUCN richness rasterize and plot
iucn_res_stack <- richness_rasterize('./data/IUCN', './data/raster/iucn', oceans_raster, 
                                     "+proj=cea +units=km")
save(iucn_res_stack, file = './data/raster/iucn_res_stack.Rdata')
iucn_richness <- richness_plot(iucn_res_stack, './figures/IUCN_maps.pdf', continents)
save(iucn_richness, file = './data/raster/iucn_richness.Rdata')
load('./data/raster/iucn_richness.Rdata')

# Carcharhiniformes rasterize and plot
car_res_stack <- richness_rasterize('./data/Carcharhiniformes', './data/raster/car', oceans_raster, 
                                    "+proj=cea +units=km")
save(car_res_stack, file = './data/raster/car_res_stack.Rdata')
car_richness <- richness_plot(car_res_stack, './figures/Carcharhiniforme_richness.pdf', continents)
save(car_richness, file = './data/raster/car_richness.Rdata')
load('./data/raster/car_richness.Rdata')

# Lamniformes rasterize and plot
lam_res_stack <- richness_rasterize('./data/Lamniformes', './data/raster/lam', oceans_raster, 
                                    "+proj=cea +units=km")
save(lam_res_stack, file = './data/raster/lam_res_stack.Rdata')
lam_richness <- richness_plot(lam_res_stack, './figures/Lamniforme_richness.pdf', continents)
save(lam_richness, file = './data/raster/lam_richness.Rdata')
load('./data/raster/lam_richness.Rdata')

# Using enviro_plot function to rasterize environmental variables
# temperature plot
load('./data/raster/temp_raster.Rdata')
temp_list <- enviro_plot(temp_raster, './figures/temperature.pdf')
save(temp_list, file = './data/raster/temp_list.Rdata')
load('./data/raster/temp_list.Rdata')

# chlorophyll plot
load('./data/raster/chloro_raster.Rdata')
chloro_list <- enviro_plot(chloro_raster, './figures/chlorophyll.pdf')
save(chloro_list, file = './data/raster/chloro_list.Rdata')
load('./data/raster/chloro_list.Rdata')

# salinity plot
load('./data/raster/salinity_raster.Rdata')
salinity_list <- enviro_plot(salinity_raster, './figures/salinity.pdf')
save(salinity_list, file = './data/raster/salinity_list.Rdata')
load('./data/raster/salinity_list.Rdata')

# bathymetry plot
load('./data/raster/bathy_raster.Rdata')
bathy_list <- enviro_plot(bathy_raster, './figures/bathymetry.pdf')
save(bathy_list, file = './data/raster/bathy_list.Rdata')
load('./data/raster/bathy_list.Rdata')

# area plot
load('./data/raster/area_raster.Rdata')
area_list <- enviro_plot(area_raster, './figures/area.pdf')
save(area_list, file = './data/raster/area_list.Rdata')
load('./data/raster/area_list.Rdata')

# distance from the coast plot
load('./data/raster/distance_raster.Rdata')
distance_list <- enviro_plot(distance_raster, './figures/distance_from_coast.pdf')
save(distance_list, file = './data/raster/distance_list.Rdata')
load('./data/raster/distance_list.Rdata')

# latitude plot
load('./data/raster/latitude_raster.Rdata')
latitude_list <- enviro_plot(latitude_raster, './figures/latitude.pdf')
save(latitude_list, file = './data/raster/latitude_list.Rdata')
load('./data/raster/latitude_list.Rdata')

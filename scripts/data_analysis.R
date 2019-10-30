library(sp)
library(raster)
load('./data/raster/species_richness.Rdata')
load('./data/raster/car_richness.Rdata')
load('./data/raster/lam_richness.Rdata')
load('./data/raster/psv_raster_list.Rdata')
load('./data/raster/car_psv_list.Rdata')
load('./data/raster/lam_psv_list.Rdata')
load('./data/raster/mrd_raster_list.Rdata')
load('./data/raster/car_mrd_list.Rdata')
load('./data/raster/lam_mrd_list.Rdata')
load('./data/raster/beta_raster_list.Rdata')
load('./data/raster/temp_list.Rdata')
load('./data/raster/bathy_list.Rdata')
load('./data/raster/chloro_list.Rdata')
load('./data/raster/salinity_list.Rdata')
load('./data/raster/area_list.Rdata')
load('./data/raster/distance_list.Rdata')

# make an array with everything in it
data_list <- vector("list", length = 6)
for (i in 1:6) {
  coo <- sp::coordinates(species_richness[[i]])
  longitude <- coo[,1]
  latitude <- coo[,2]
  dat <- data.frame(latitude, longitude, values(species_richness[[i]]), 
                values(car_richness[[i]]), values(lam_richness[[i]]), 
                values(psv_raster_list[[i]]), values(car_psv_list[[i]]), 
                values(lam_psv_list[[i]]), values(mrd_raster_list[[i]]),
                values(car_mrd_list[[i]]), values(lam_mrd_list[[i]]),
                values(beta_raster_list[[i]]), values(temp_list[[i]]), 
                values(bathy_list[[i]]), values(chloro_list[[i]]),
                values(salinity_list[[i]]), values(area_list[[i]]),
                values(distance_list[[i]])) 
  colnames(dat) <- c("latitude", "longitude", "total_species_richness", 
                     "carcharhiniforme_species_richness", 
                     "lamniformes_species_richness", "total_psv", 
                     "carcharhiniforme_psv", "lamniforme_psv", "total_mrd",
                     "carcharhiniforme_mrd", "lamniforme_mrd", "total_beta",
                     "temperature", "bathymetry", "chlorophyll", "salinity", 
                     "area", "distance_from_coast")
  data_list[[i]] <- dat
}
save(data_list, file = './data/data_list.Rdata')
load('./data/data_list.Rdata')
  
# make_plot function plots two variables against one another
# x_var = x axis variable, must be chosen from data_list, in character form
# y_var = y axis variable, must be chosen from data_list, in character form
# x_lab = x axis label as character string
# y_lab = y axis label as character string
# figure_name = pdf file
# output are statistics summary, name and save accordingly
make_plot <- function(x_var, y_var, x_lab, y_lab, figure_name) {
   load('./data/data_list.Rdata')
   load('./data/raster/species_richness.Rdata')
   stats <- vector("list", length = 6)
   pdf(figure_name)
   for (i in 1:6) {
        column1 <- substitute(x_var)
        column2 <- substitute(y_var)
        lin <- lm(eval(column2, data_list[[i]]) ~ eval(column1, data_list[[i]])) 
        plot(eval(column1, data_list[[i]]), eval(column2, data_list[[i]]), 
             main = paste('resolution =', res(species_richness[[i]])), 
             xlab = x_lab, ylab = y_lab)
        abline(lin, col = 'red')
        stats[[i]] <- summary(lin)
    }
dev.off()
return(stats)
}


# temperature vs total richness
tempVrichness_stats <- make_plot(temperature, total_species_richness,
                                 "Temperature (°C)", "Species Richness",
                                 './figures/temperature_vs_richness.pdf')
save(tempVrichness_stats, file = './data/stats/tempVrichness_stats.Rdata')
load('./data/stats/tempVrichness_stats.Rdata')

# temperature vs total psv
tempVpsv_stats <- make_plot(temperature, total_psv,
                                 "Temperature (°C)", 
                                 "Phylogenetic Species Richness (PSV)",
                                 './figures/temperature_vs_psv.pdf')
save(tempVpsv_stats, file = './data/stats/tempVpsv_stats.Rdata')
load('./data/stats/tempVpsv_stats.Rdata')

# temperature vs total mrd
tempVmrd_stats <- make_plot(temperature, total_mrd,
                            "Temperature (°C)", 
                            "Mean Root Distance",
                            './figures/temperature_vs_mrd.pdf')
save(tempVmrd_stats, file = './data/stats/tempVmrd_stats.Rdata')
load('./data/stats/tempVmrd_stats.Rdata')

# temperature vs total beta
tempVbeta_stats <- make_plot(temperature, total_beta,
                            "Temperature (°C)", 
                            "Beta",
                            './figures/temperature_vs_beta.pdf')
save(tempVbeta_stats, file = './data/stats/tempVbeta_stats.Rdata')
load('./data/stats/tempVbeta_stats.Rdata')


# bathymetry vs total richness
bathVrichness_stats <- make_plot(bathymetry, total_species_richness,
                                 "Depth (m)", "Species Richness",
                                 './figures/bath_vs_richness.pdf')
save(bathVrichness_stats, file = './data/stats/bathVrichness_stats.Rdata')
load('./data/stats/bathVrichness_stats.Rdata')

# bathymetry vs total psv
bathVpsv_stats <- make_plot(bathymetry, total_psv,
                            "Depth (m)", 
                            "Phylogenetic Species Richness (PSV)",
                            './figures/bath_vs_psv.pdf')
save(bathVpsv_stats, file = './data/stats/bathVpsv_stats.Rdata')
load('./data/stats/bathVpsv_stats.Rdata')

# bathymetry vs total mrd
bathVmrd_stats <- make_plot(bathymetry, total_mrd,
                            "Depth (m)", 
                            "Mean Root Distance (MRD)",
                            './figures/bath_vs_mrd.pdf')
save(bathVmrd_stats, file = './data/stats/bathVmrd_stats.Rdata')
load('./data/stats/bathVmrd_stats.Rdata')

# bathymetry vs total beta
bathVbeta_stats <- make_plot(bathymetry, total_beta,
                            "Depth (m)", 
                            "Beta",
                            './figures/bath_vs_beta.pdf')
save(bathVbeta_stats, file = './data/stats/bathVbeta_stats.Rdata')
load('./data/stats/bathVbeta_stats.Rdata')


# chlorophyll vs total richness
chloroVrichness_stats <- make_plot(chlorophyll, total_species_richness,
                                 "Chlorophyll (mg/m3)", "Species Richness",
                                 './figures/chlorophyll_vs_richness.pdf')
save(chloroVrichness_stats, file = './data/stats/chloroVrichness_stats.Rdata')
load('./data/stats/chloroVrichness_stats.Rdata')

# chlorophyll vs total psv
chloroVpsv_stats <- make_plot(chlorophyll, total_psv,
                                   "Chlorophyll (mg/m3)", 
                                   "Phylogenetic Species Varianc (PSV)",
                                   './figures/chlorophyll_vs_psv.pdf')
save(chloroVpsv_stats, file = './data/stats/chloroVpsv_stats.Rdata')
load('./data/stats/chloroVpsv_stats.Rdata')

# chlorophyll vs total mrd
chloroVmrd_stats <- make_plot(chlorophyll, total_mrd,
                              "Chlorophyll (mg/m3)", 
                              "Mean Root Distance (MRD)",
                              './figures/chlorophyll_vs_mrd.pdf')
save(chloroVmrd_stats, file = './data/stats/chloroVmrd_stats.Rdata')
load('./data/stats/chloroVmrd_stats.Rdata')

# chlorophyll vs total beta
chloroVbeta_stats <- make_plot(chlorophyll, total_beta,
                              "Chlorophyll (mg/m3)", "Beta",
                              './figures/chlorophyll_vs_beta.pdf')
save(chloroVbeta_stats, file = './data/stats/chloroVbeta_stats.Rdata')
load('./data/stats/chloroVbeta_stats.Rdata')


# salinity vs total richness
salinityVrichness_stats <- make_plot(salinity, total_species_richness,
                                   "Salinity", "Species Richness",
                                   './figures/salinity_vs_richness.pdf')
save(salinityVrichness_stats, file = './data/stats/salinityVrichness_stats.Rdata')
load('./data/stats/salinityVrichness_stats.Rdata')

# salinity vs total psv
salinityVpsv_stats <- make_plot(salinity, total_psv,
                                     "Salinity", 
                                     "Phylogenetic Species Variance (PSV)",
                                     './figures/salinity_vs_psv.pdf')
save(salinityVpsv_stats, file = './data/stats/salinityVpsv_stats.Rdata')
load('./data/stats/salinityVpsv_stats.Rdata')

# salinity vs total mrd
salinityVmrd_stats <- make_plot(salinity, total_mrd,
                                "Salinity", 
                                "Mean Root Distance (MRD)",
                                './figures/salinity_vs_mrd.pdf')
save(salinityVmrd_stats, file = './data/stats/salinityVmrd_stats.Rdata')
load('./data/stats/salinityVmrd_stats.Rdata')

# salinity vs total beta
salinityVbeta_stats <- make_plot(salinity, total_beta,
                                "Salinity", "Beta",
                                './figures/salinity_vs_beta.pdf')
save(salinityVbeta_stats, file = './data/stats/salinityVbeta_stats.Rdata')
load('./data/stats/salinityVbeta_stats.Rdata')


# latitude vs total richness
latVrichness_stats <- make_plot(abs(latitude), total_species_richness,
                                   "Latitude", "Species Richness",
                                   './figures/latitude_vs_richness.pdf')
save(latVrichness_stats, file = './data/stats/latVrichness_stats.Rdata')
load('./data/stats/latVrichness_stats.Rdata')

# latitude vs total psv
latVpsv_stats <- make_plot(latitude, total_psv,
                                "Latitude", "Phylogentic Species Variance (PSV)",
                                './figures/latitude_vs_psv.pdf')
save(latVpsv_stats, file = './data/stats/latVpsv_stats.Rdata')
load('./data/stats/latVpsv_stats.Rdata')

# latitude vs total mrd
latVmrd_stats <- make_plot(latitude, total_mrd,
                                "Latitude", "Mean Root Distance (MRD)",
                                './figures/latitude_vs_mrd.pdf')
save(latVmrd_stats, file = './data/stats/latVmrd_stats.Rdata')
load('./data/stats/latVmrd_stats.Rdata')

# latitude vs total beta
latVbeta_stats <- make_plot(latitude, total_beta,
                                "Latitude", "Beta",
                                './figures/latitude_vs_beta.pdf')
save(latVbeta_stats, file = './data/stats/latVbeta_stats.Rdata')
load('./data/stats/latVbeta_stats.Rdata')


# distance vs total richness
distanceVrichness_stats <- make_plot(log(distance_from_coast+1), total_species_richness,
                                   "Distance From Coast", "Species Richness",
                                   './figures/distance_vs_richness.pdf')
save(distanceVrichness_stats, file = './data/stats/distanceVrichness_stats.Rdata')
load('./data/stats/distanceVrichness_stats.Rdata')

# distance vs total psv
distanceVpsv_stats <- make_plot(distance_from_coast, total_psv,
                                     "Distance From Coast", 
                                     "Phylogenetic Species Variance (PSV)",
                                     './figures/distance_vs_psv.pdf')
save(distanceVpsv_stats, file = './data/stats/distanceVpsv_stats.Rdata')
load('./data/stats/distanceVpsv_stats.Rdata')

# distance vs total mrd
distanceVmrd_stats <- make_plot(distance_from_coast, total_mrd,
                                "Distance From Coast", 
                                "Mean Root Distance (MRD)",
                                './figures/distance_vs_mrd.pdf')
save(distanceVmrd_stats, file = './data/stats/distanceVmrd_stats.Rdata')
load('./data/stats/distanceVmrd_stats.Rdata')

# distance vs total beta
distanceVbeta_stats <- make_plot(distance_from_coast, total_beta,
                                "Distance From Coast", 
                                "Beta",
                                './figures/distance_vs_beta.pdf')
save(distanceVbeta_stats, file = './data/stats/distanceVbeta_stats.Rdata')
load('./data/stats/distanceVbeta_stats.Rdata')


# area vs total richness
areaVrichness_stats <- make_plot(area, total_species_richness,
                                   "Area", "Species Richness",
                                   './figures/area_vs_richness.pdf')
save(areaVrichness_stats, file = './data/stats/areaVrichness_stats.Rdata')
load('./data/stats/areaVrichness_stats.Rdata')

# area vs total psv
areaVpsv_stats <- make_plot(area, total_psv,
                                 "Area", "Phylogenetic Species Variance (PSV)",
                                 './figures/area_vs_psv.pdf')
save(areaVpsv_stats, file = './data/stats/areaVpsv_stats.Rdata')
load('./data/stats/areaVpsv_stats.Rdata')

# area vs total mrd
areaVmrd_stats <- make_plot(area, total_mrd,
                            "Area", "Mean Root Distance (MRD)",
                            './figures/area_vs_mrd.pdf')
save(areaVmrd_stats, file = './data/stats/areaVmrd_stats.Rdata')
load('./data/stats/areaVmrd_stats.Rdata')

# area vs total beta
areaVbeta_stats <- make_plot(area, total_beta,
                            "Area", "Beta",
                            './figures/area_vs_beta.pdf')
save(areaVbeta_stats, file = './data/stats/areaVbeta_stats.Rdata')
load('./data/stats/areaVbeta_stats.Rdata')


# total richness vs total mrd
mrdVrichness_stats <- make_plot(total_mrd, total_species_richness,
                                     "Mean Root Distance", "Species Richness",
                                     './figures/mrd_vs_richness.pdf')
save(mrdVrichness_stats, file = './data/stats/mrdVrichness_stats.Rdata')
load('./data/stats/mrdVrichness_stats.Rdata')

# total mrd vs total psv
psvVmrd_stats <- make_plot(total_psv, total_mrd,
                                "Phylogenetic Species Variance (PSV)", 
                                "Mean Root Distance (MRD)",
                                './figures/psv_vs_mrd.pdf')
save(psvVmrd_stats, file = './data/stats/psvVmrd_stats.Rdata')
load('./data/stats/psvVmrd_stats.Rdata')

# total psv vs total richness
psvVrichness_stats <- make_plot(total_psv, total_species_richness,
                                "Phylogenetic Species Variance (PSV)", 
                                "Species Richness",
                                './figures/psv_vs_richness.pdf')
save(psvVrichness_stats, file = './data/stats/psvVrichness_stats.Rdata')
load('./data/stats/psvVrichness_stats.Rdata')



# all for Carcharhiniformes
# temperature vs richness
tempVcarrichness_stats <- make_plot(temperature, carcharhiniforme_species_richness,
                                 "Temperature (°C)", 
                                 "Carcharhiniformes Species Richness",
                                 './figures/temperature_vs_Car_richness.pdf')
save(tempVcarrichness_stats, file = './data/stats/tempVcarrichness_stats.Rdata')
load('./data/stats/tempVcarrichness_stats.Rdata')

# temperature vs psv
tempVcarpsv_stats <- make_plot(temperature, carcharhiniforme_psv,
                            "Temperature (°C)", 
                            "Phylogenetic Species Richness (PSV) Carcharhiniformes",
                            './figures/temperature_vs_Car_psv.pdf')
save(tempVcarpsv_stats, file = './data/stats/tempVcarpsv_stats.Rdata')
load('./data/stats/tempVcarpsv_stats.Rdata')

# temperature vs mrd
tempVcarmrd_stats <- make_plot(temperature, carcharhiniforme_mrd,
                            "Temperature (°C)", 
                            "Mean Root Distance (MRD) Carcharhiniformes",
                            './figures/temperature_vs_Car_mrd.pdf')
save(tempVcarmrd_stats, file = './data/stats/tempVcarmrd_stats.Rdata')
load('./data/stats/tempVcarmrd_stats.Rdata')


# bathymetry vs richness
bathVcarrichness_stats <- make_plot(bathymetry, carcharhiniforme_species_richness,
                                 "Depth (m)", 
                                 "Carcharhiniformes Species Richness",
                                 './figures/bath_vs_Car_richness.pdf')
save(bathVcarrichness_stats, file = './data/stats/bathVcarrichness_stats.Rdata')
load('./data/stats/bathVcarrichness_stats.Rdata')

# bathymetry vs psv
bathVcarpsv_stats <- make_plot(bathymetry, carcharhiniforme_psv,
                            "Depth (m)", 
                            "Phylogenetic Species Richness (PSV) carcharhiniformes",
                            './figures/bath_vs_Car_psv.pdf')
save(bathVcarpsv_stats, file = './data/stats/bathVcarpsv_stats.Rdata')
load('./data/stats/bathVcarpsv_stats.Rdata')

# bathymetry vs mrd
bathVcarmrd_stats <- make_plot(bathymetry, carcharhiniforme_mrd,
                            "Depth (m)", 
                            "Mean Root Distance (MRD) Carcharhiniformes",
                            './figures/bath_vs_Car_mrd.pdf')
save(bathVcarmrd_stats, file = './data/stats/bathVcarmrd_stats.Rdata')
load('./data/stats/bathVcarmrd_stats.Rdata')


# distance vs richness
distanceVcarrichness_stats <- make_plot(distance_from_coast, 
                                        carcharhiniforme_species_richness,
                                        "Distance (km)", 
                                        "Carcharhiniformes Species Richness",
                                    './figures/distance_vs_Car_richness.pdf')
save(distanceVcarrichness_stats, file = './data/stats/distanceVcarrichness_stats.Rdata')
load('./data/stats/distanceVcarrichness_stats.Rdata')

# distance vs psv
distanceVcarpsv_stats <- make_plot(distance_from_coast, carcharhiniforme_psv,
                               "Distance (km)", 
                               "Phylogenetic Species Richness (PSV) Carcharhiniformes",
                               './figures/distance_vs_Car_psv.pdf')
save(distanceVcarpsv_stats, file = './data/stats/distanceVcarpsv_stats.Rdata')
load('./data/stats/distanceVcarpsv_stats.Rdata')

# distance vs mrd
distanceVcarmrd_stats <- make_plot(distance_from_coast, carcharhiniforme_mrd,
                               "Distance (km)", 
                               "Mean Root Distance (MRD) Carcharhiniformes",
                               './figures/distance_vs_Car_mrd.pdf')
save(distanceVcarmrd_stats, file = './data/stats/distanceVcarmrd_stats.Rdata')
load('./data/stats/distanceVcarmrd_stats.Rdata')


# richness vs mrd
carmrdVcarrichness_stats <- make_plot(carcharhiniforme_mrd, 
                                      carcharhiniforme_species_richness,
                                "Mean Root Distance (MRD) Carcharhiniformes", 
                                "Carcharhiniformes Species Richness",
                                './figures/Car_mrd_vs_Car_richness.pdf')
save(carmrdVcarrichness_stats, file = './data/stats/carmrdVcarrichness_stats.Rdata')
load('./data/stats/carmrdVcarrichness_stats.Rdata')

# mrd vs psv
carpsvVcarmrd_stats <- make_plot(carcharhiniforme_psv, carcharhiniforme_mrd,
                           "Phylogenetic Species Variance (PSV) Carcharhiniformes", 
                           "Mean Root Distance (MRD) Carcharhiniformes",
                           './figures/Car_psv_vs_Car_mrd.pdf')
save(carpsvVcarmrd_stats, file = './data/stats/carpsvVcarmrd_stats.Rdata')
load('./data/stats/carpsvVcarmrd_stats.Rdata')

# psv vs richness
carpsvVcarrichness_stats <- make_plot(carcharhiniforme_psv, 
                                      carcharhiniforme_species_richness,
                                "Phylogenetic Species Variance (PSV) Carcharhiniformes", 
                                "Carcharhiniformes Species Richness",
                                './figures/Car_psv_vs_Car_richness.pdf')
save(carpsvVcarrichness_stats, file = './data/stats/carpsvVcarrichness_stats.Rdata')
load('./data/stats/carpsvVcarrichness_stats.Rdata')



# All with Lamniformes
# temperature vs richness
tempVlamrichness_stats <- make_plot(temperature, lamniformes_species_richness,
                                    "Temperature (°C)", 
                                    "Lamniformes Species Richness",
                                    './figures/temperature_vs_Lam_richness.pdf')
save(tempVlamrichness_stats, file = './data/stats/tempVlamrichness_stats.Rdata')
load('./data/stats/tempVlamrichness_stats.Rdata')

# temperature vs psv
tempVlampsv_stats <- make_plot(temperature, lamniforme_psv,
                               "Temperature (°C)", 
                               "Phylogenetic Species Richness (PSV) Lamniformes",
                               './figures/temperature_vs_Lam_psv.pdf')
save(tempVlampsv_stats, file = './data/stats/tempVlampsv_stats.Rdata')
load('./data/stats/tempVlampsv_stats.Rdata')

# temperature vs mrd
tempVlammrd_stats <- make_plot(temperature, lamniforme_mrd,
                               "Temperature (°C)", 
                               "Mean Root Distance (MRD) Lamniformes",
                               './figures/temperature_vs_Lam_mrd.pdf')
save(tempVlammrd_stats, file = './data/stats/tempVlammrd_stats.Rdata')
load('./data/stats/tempVlammrd_stats.Rdata')


# bathymetry vs richness
bathVlamrichness_stats <- make_plot(bathymetry, lamniformes_species_richness,
                                    "Depth (m)", 
                                    "Lamniformes Species Richness",
                                    './figures/bath_vs_Lam_richness.pdf')
save(bathVlamrichness_stats, file = './data/stats/bathVlamrichness_stats.Rdata')
load('./data/stats/bathVlamrichness_stats.Rdata')

# bathymetry vs psv
bathVlampsv_stats <- make_plot(bathymetry, lamniforme_psv,
                               "Depth (m)", 
                               "Phylogenetic Species Richness (PSV) Lamniformes",
                               './figures/bath_vs_Lam_psv.pdf')
save(bathVlampsv_stats, file = './data/stats/bathVlampsv_stats.Rdata')
load('./data/stats/bathVlampsv_stats.Rdata')

# bathymetry vs mrd
bathVlammrd_stats <- make_plot(bathymetry, lamniforme_mrd,
                               "Depth (m)", 
                               "Mean Root Distance (MRD) Lamniformes",
                               './figures/bath_vs_Lam_mrd.pdf')
save(bathVlammrd_stats, file = './data/stats/bathVlammrd_stats.Rdata')
load('./data/stats/bathVlammrd_stats.Rdata')


# distance vs richness
distanceVlamrichness_stats <- make_plot(distance_from_coast, 
                                        lamniformes_species_richness,
                                        "Distance (km)", 
                                        "Lamniformes Species Richness",
                                        './figures/distance_vs_Lam_richness.pdf')
save(distanceVlamrichness_stats, file = './data/stats/distanceVlamrichness_stats.Rdata')
load('./data/stats/distanceVlamrichness_stats.Rdata')

# distance vs psv
distanceVlampsv_stats <- make_plot(distance_from_coast, lamniforme_psv,
                                   "Distance (km)", 
                                   "Phylogenetic Species Richness (PSV) Lamniformes",
                                   './figures/distance_vs_Lam_psv.pdf')
save(distanceVlampsv_stats, file = './data/stats/distanceVlampsv_stats.Rdata')
load('./data/stats/distanceVlampsv_stats.Rdata')

# distance vs mrd
distanceVlammrd_stats <- make_plot(distance_from_coast, lamniforme_mrd,
                                   "Distance (km)", 
                                   "Mean Root Distance (MRD) Lamniformes",
                                   './figures/distance_vs_Lam_mrd.pdf')
save(distanceVlammrd_stats, file = './data/stats/distanceVlammrd_stats.Rdata')
load('./data/stats/distanceVlammrd_stats.Rdata')


# richness vs mrd
lammrdVlamrichness_stats <- make_plot(lamniforme_mrd, 
                                      lamniformes_species_richness,
                                      "Mean Root Distance (MRD) Lamniformes", 
                                      "Lamniformes Species Richness",
                                      './figures/Lam_mrd_vs_Lam_richness.pdf')
save(lammrdVlamrichness_stats, file = './data/stats/lammrdVlamrichness_stats.Rdata')
load('./data/stats/lammrdVlamrichness_stats.Rdata')

# mrd vs psv
lampsvVlammrd_stats <- make_plot(lamniforme_psv, lamniforme_mrd,
                                 "Phylogenetic Species Variance (PSV) Lamniformes", 
                                 "Mean Root Distance (MRD) Lamniformes",
                                 './figures/Lam_psv_vs_Lam_mrd.pdf')
save(lampsvVlammrd_stats, file = './data/stats/lampsvVlammrd_stats.Rdata')
load('./data/stats/lampsvVlammrd_stats.Rdata')

# psv vs richness
lampsvVlamrichness_stats <- make_plot(lamniforme_psv, 
                                      lamniformes_species_richness,
                                      "Phylogenetic Species Variance (PSV) Lamniformes", 
                                      "Lamniformes Species Richness",
                                      './figures/Lam_psv_vs_Lam_richness.pdf')
save(lampsvVlamrichness_stats, file = './data/stats/lampsvVlamrichness_stats.Rdata')
load('./data/stats/lampsvVlamrichness_stats.Rdata')


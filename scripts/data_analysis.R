library(spatialEco)
library(raster)

# linear regression temp vs richness
richnessVtemp_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(species_richness_list[[i]]), values(temp_list[[i]]), 
          xlab = "Shark Richness", ylab = "Temperature (°C)")
     abline(lm(values(temp_list[[i]]) ~ values(species_richness_list[[i]])), 
            col = 'red')
     richnessVtemp <- lm(values(temp_list[[i]]) ~ 
                           values(species_richness_list[[i]]))
     richnessVtemp_list[[i]] <- richnessVtemp
     print(summary(richnessVtemp_list[[i]]))
}

# linear regression chlorophyll vs richness
richnessVchloro_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(species_richness_list[[i]]), values(chloro_list[[i]]), 
       xlab = "Shark Richness", ylab = "Chlorophyll (mg/m^3)")
     abline(lm(values(chloro_list[[i]]) ~ values(species_richness_list[[i]])), 
         col = 'red')
     richnessVchloro <- lm(values(chloro_list[[i]]) ~ 
                        values(species_richness_list[[i]]))
     richnessVchloro_list[[i]] <- richnessVchloro
     print(summary(richnessVchloro_list[[i]]))
}

# linear regression normal vs threatened
normalVthreatened_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(species_richness_list[[i]]), values(IUCN_richness_list[[i]]), 
          xlab = "Normal Shark Richness", ylab = "Threatened Shark Richness")
     abline(lm(values(IUCN_richness_list[[i]]) ~ 
                 values(species_richness_list[[i]])), col = 'red')
     normalVthreatened <- lm(values(IUCN_richness_list[[i]]) ~ 
                               values(species_richness_list[[i]]))
     normalVthreatened_list[[i]] <- normalVthreatened
     print(normalVthreatened_list[[i]])
}

# latitude vs richness
latVrichness_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(abs(latitude_list[[i]]), values(species_richness_list[[i]]), 
          xlab = "Latitude", ylab = "Shark Richness")
     abline(lm(values(species_richness_list[[i]]) ~ abs(latitude_list[[i]])), 
            col = 'red')
     latVrichness <- lm(values(species_richness_list[[i]]) ~ 
                          abs(latitude_list[[i]]))
     latVrichness_list[[i]] <- latVrichness
     print(latVrichness_list[[i]])
}

# salinity vs richness
salinityVrichness_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(species_richness_list[[i]]), values(salinity_list[[i]]), 
          xlab = "Shark Richness", ylab = "Salinity")
     abline(lm(values(salinity_list[[i]]) ~ 
                 values(species_richness_list[[i]])), col = 'red')
     salinityVrichness <- lm(values(salinity_list[[i]]) ~ 
                               values(species_richness_list[[i]]))
     salinityVrichness_list[[i]] <- salinityVrichness
     print(salinityVrichness_list[[i]])
}

# experimentation
stat_stack <- stack(species_richness_list[[1]], temp_list[[1]])
richnessVtemp <- layerStats(stat_stack, 'pearson', na.rm = T)
richnessVtemp

richnessVtemp <- rasterCorrelation(species_richness_list[[1]], temp_list[[1]], 
                                   s = 3, type = "spearman")
plot(richnessVtemp)

# basic correlation
richness_temp_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     richnessVtemp_c <- cor.test(values(species_richness_list[[i]]), 
                            values(temp_list[[i]]), alternative = "two.sided", 
                            method = "pearson")
     richness_temp_list[[i]] <- richnessVtemp_c
}

richness_chloro_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     richnessVchloro_c <- cor.test(values(species_richness_list[[i]]), 
                              values(chloro_list[[i]]), 
                              alternative = "two.sided", method = "pearson")
     richness_chloro_list[[i]] <- richnessVchloro_c
}

normal_threatened_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     normalVthreatened_c <- cor.test(values(species_richness_list[[i]]), 
                                values(IUCN_richness_list[[i]]), 
                                alternative = "two.sided", method = "pearson")
     normal_threatened_list[[i]] <- normalVthreatened_c
}

richness_latitude_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     richnessVlat_c <- cor.test(values(species_richness_list[[i]]), 
                           abs(latitude_list[[i]]), alternative = "two.sided", 
                           method = "pearson")
     richness_latitude_list[[i]] <- richnessVlat_c
}

richness_salinity_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2 ,3 ,4 , 5, 6, 7))) {
     richnessVsalinity_c <- cor.test(values(species_richness_list[[i]]), 
                                     values(salinity_list[[i]]), 
                                     alternative = "two.sided", 
                                     method = "pearson")
     richness_salinity_list[[i]] <- richnessVsalinity_c
} 
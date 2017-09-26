library(spatialEco)
library(raster)

# linear regression temp vs richness
pdf('./figures/temperature_vs_richness.pdf')
richnessVtemp_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(temp_list[[i]]), values(species_richness_list[[i]]),
          main = paste('resolution =', res(res_list[[i]])), 
          xlab = "Temperature (°C)", ylab = "Shark Richness")
     abline(lm(values(species_richness_list[[i]]) ~ values(temp_list[[i]])), 
            col = 'red')
     richnessVtemp <- lm(values(species_richness_list[[i]]) ~ 
                           values(temp_list[[i]]))
     richnessVtemp_list[[i]] <- richnessVtemp
     print(summary(richnessVtemp_list[[i]]))
}
dev.off()

# linear regression chlorophyll vs richness
pdf('./figures/chlorophyll_vs_richness.pdf')
richnessVchloro_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(chloro_list[[i]]), values(species_richness_list[[i]]), 
          main = paste('resolution =', res(res_list[[i]])), 
          xlab = "Chlorophyll (mg/m3)", ylab = "Shark Richness")
     abline(lm(values(species_richness_list[[i]]) ~ values(chloro_list[[i]])), 
            col = 'red')
     richnessVchloro <- lm(values(species_richness_list[[i]]) ~ 
                           values(chloro_list[[i]]))
     richnessVchloro_list[[i]] <- richnessVchloro
     print(summary(richnessVchloro_list[[i]]))
}
dev.off()

# linear regression normal vs threatened
pdf('./figures/IUSN_vs_Normal.pdf')
normalVthreatened_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(species_richness_list[[i]]), values(IUCN_richness_list[[i]]), 
          main = paste('resolution =', res(res_list[[i]])), 
          xlab = "Normal Shark Richness", ylab = "Threatened Shark Richness")
     abline(lm(values(IUCN_richness_list[[i]]) ~ 
            values(species_richness_list[[i]])), col = 'red')
     normalVthreatened <- lm(values(IUCN_richness_list[[i]]) ~ 
                             values(species_richness_list[[i]]))
     normalVthreatened_list[[i]] <- normalVthreatened
     print(normalVthreatened_list[[i]])
}
dev.off()

# latitude vs richness
pdf('./figures/latitude_vs_richness.pdf')
latVrichness_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(abs(latitude_list[[i]]), values(species_richness_list[[i]]), 
          main = paste('resolution =', res(res_list[[i]])), xlab = "Latitude", 
          ylab = "Shark Richness")
     abline(lm(values(species_richness_list[[i]]) ~ abs(latitude_list[[i]])), 
            col = 'red')
     latVrichness <- lm(values(species_richness_list[[i]]) ~ 
                        abs(latitude_list[[i]]))
     latVrichness_list[[i]] <- latVrichness
     print(latVrichness_list[[i]])
}
dev.off()

# salinity vs richness
pdf('./figures/salinity_vs_richness.pdf')
salinityVrichness_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(salinity_list[[i]]), values(species_richness_list[[i]]), 
          main = paste('resolution =', res(res_list[[i]])), xlab = "Salinity", 
          ylab = "Shark Richness")
     abline(lm(values(species_richness_list[[i]]) ~ 
                 values(salinity_list[[i]])), col = 'red')
     salinityVrichness <- lm(values(species_richness_list[[i]]) ~ 
                               values(salinity_list[[i]]))
     salinityVrichness_list[[i]] <- salinityVrichness
     print(salinityVrichness_list[[i]])
}
dev.off()

# distance from coast vs richness normal
pdf('./figures/coast_vs_richness_norm.pdf')
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(coast_distance_list[[i]]), values(species_richness_list[[i]]), 
          main = paste('resolution =', res(res_list[[i]])), xlab = "Distance", 
          ylab = "Shark Richness")
}
dev.off()

# distance from coast vs richness log
pdf('./figures/coast_vs_richness_log.pdf')
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
  plot(log(values(coast_distance_list[[i]])), values(species_richness_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), xlab = "Distance", 
       ylab = "Shark Richness")
  abline(lm(values(species_richness_list[[i]]) ~ 
              log(values(coast_distance_list[[i]]))), col = 'red')
}
dev.off()

# taxonomic richness vs faith's pd
pdf('./figures/richness_vs_pd.pdf')
for (i in 1:7) {
  plot(values(pd_raster_list[[i]]), values(species_richness_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), 
       xlab = "Faith's Phylogentic Diversity", 
       ylab = "Taxonomic Richness")
  abline(lm(values(species_richness_list[[i]]) ~ 
              values(pd_raster_list[[i]])), col = 'red')
}
dev.off()

# taxonomic richness vs psv
pdf('./figures/richness_vs_psv.pdf')
for (i in 1:7) {
  plot(values(psv_raster_list[[i]]), values(species_richness_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), 
       xlab = "PSV", 
       ylab = "Taxonomic Richness")
  abline(lm(values(species_richness_list[[i]]) ~ 
              values(psv_raster_list[[i]])), col = 'red')
}
dev.off()

# taxonomic richness vs psr
pdf('./figures/richness_vs_psr.pdf')
for (i in 1:7) {
  plot(values(psr_raster_list[[i]]), values(species_richness_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), 
       xlab = "PSR", 
       ylab = "Taxonomic Richness")
  abline(lm(values(species_richness_list[[i]]) ~ 
              values(psr_raster_list[[i]])), col = 'red')
}
dev.off()

# taxonomic richness vs random pd
pdf('./figures/richness_vs_randompd.pdf')
for (i in 1:7) {
  plot(values(randompd_raster_list[[i]]), values(species_richness_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), 
       xlab = "Random PD", 
       ylab = "Taxonomic Richness")
  abline(lm(values(species_richness_list[[i]]) ~ 
              values(randompd_raster_list[[i]])), col = 'red')
}
dev.off()

# normal psv vs random psv
pdf('./figures/normalpsv_vs_randompsv.pdf')
for (i in 1:7) {
  plot(values(psv_raster_list[[i]]), values(randompsv_raster_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), 
       xlab = "Normal PSV", 
       ylab = "Random PSV")
  abline(lm(values(randompsv_raster_list[[i]]) ~ 
              values(psv_raster_list[[i]])), col = 'red')
}
dev.off()

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




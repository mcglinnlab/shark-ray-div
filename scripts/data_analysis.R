library(spatialEco)
library(raster)

# necessary variables
load('./data/raster/species_richness.Rdata')
load('./data/raster/temp_list.Rdata')
load('./data/raster/chloro_list.Rdata')
load('./data/raster/IUCN_richness_list.Rdata')
load('./data/raster/salinity_list')
load('./data/raster/mrd_raster_list.Rdata')
load('./data/raster/pd_raster_list.Rdata')
load('./data/raster/psv_raster_list.Rdata')
load('./data/raster/psr_raster_list.Rdata')
load('./data/raster/area_list.Rdata')
load('./data/raster/beta_raster_list.Rdata')

# linear regression temp and area vs richness
temp_stats <- vector("list", length = 6)
pdf('./figures/temperature_and_area_vs_richness.pdf')
for (i in 1:6) {
     temp <- lm(values(species_richness[[i]]) ~ values(temp_list[[i]]) + 
                  values(area_list[[i]]))
     plot(values(temp_list[[i]]) + values(area_list[[i]]), 
          values(species_richness[[i]]),
          main = paste('resolution =', res(species_richness[[i]])), 
          xlab = "Temperature (°C) vs Area", ylab = "Shark Richness")
     temp_stats[[i]] <- summary(temp)
}
dev.off()
save(temp_stats, file = './data/stats/temp_stats.Rdata')
load('./data/stats/temp_stats.Rdata')

# temp and area vs psv
tempVpsv_stats <- vector("list", length = 6)
pdf('./figures/temperature_and_area_vs_psv.pdf')
for (i in 1:6) {
  temp <- lm(values(psv_raster_list[[i]]) ~ values(temp_list[[i]]) + 
               values(area_list[[i]]))
  plot(values(temp_list[[i]]) + values(area_list[[i]]), 
       values(psv_raster_list[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C) vs Area", ylab = "Phylogenetic Species Variance (PSV)")
  tempVpsv_stats[[i]] <- summary(temp)
}
dev.off()
save(tempVpsv_stats, file = './data/stats/tempVpsv_stats.Rdata')
load('./data/stats/tempVpsv_stats.Rdata')

# temp and area vs mrd
tempVmrd_stats <- vector("list", length = 6)
pdf('./figures/temperature_and_area_vs_mrd.pdf')
for (i in 1:6) {
  temp <- lm(values(mrd_raster_list[[i]]) ~ values(temp_list[[i]]) + 
               values(area_list[[i]]))
  plot(values(temp_list[[i]]) + area_list[[i]], values(mrd_raster_list[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C) vs Area", ylab = "mrd")
  tempVmrd_stats <- temp
}
dev.off()
save(tempVmrd_stats, file = './data/stats/tempVmrd_stats.Rdata')
load('./data/stats/tempVmrd_stats.Rdata')

# temp and area vs beta
tempVbeta_stats <- vector("list", length = 6)
pdf('./figures/temperature_and_area_vs_beta.pdf')
for (i in 1:6) {
  temp <- lm(values(beta_raster_list[[i]]) ~ values(temp_list[[i]]) + 
               values(area_list[[i]]))
  plot(values(temp_list[[i]]) + area_list[[i]], values(beta_raster_list[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C) vs Area", ylab = "beta")
  tempVbeta_stats[[i]] <- temp
}
dev.off()
save(tempVbeta_stats, file = './data/stats/tempVbeta_stats.Rdata')
load('./data/stats/tempVbeta_stats.Rdata')

# linear regression temp vs richness
tempVrichness_stats <- vector("list", length = 6)
pdf('./figures/temperature_vs_richness.pdf')
for (i in 1:6) {
  temp <- lm(values(species_richness[[i]]) ~ values(temp_list[[i]])) 
  plot(values(temp_list[[i]]), 
       values(species_richness[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C)", ylab = "Shark Richness")
  abline(lm(values(species_richness[[i]]) ~ values(temp_list[[i]])), col = 'red')
  tempVrichness_stats[[i]] <- summary(temp)
}
dev.off()
save(tempVrichness_stats, file = './data/stats/tempVrichness_stats.Rdata')
load('./data/stats/tempVrichness_stats.Rdata')

# temp vs psv
tempVpsv_stats_alone <- vector("list", length = 6)
pdf('./figures/temperature_vs_psv.pdf')
for (i in 1:6) {
  temp <- lm(values(psv_raster_list[[i]]) ~ values(temp_list[[i]]))
  plot(values(temp_list[[i]]), 
       values(psv_raster_list[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C) vs Area", ylab = "Phylogenetic Species Variance (PSV)")
  abline(lm(values(psv_raster_list[[i]]) ~ values(temp_list[[i]])), col = 'red')
  tempVpsv_stats_alone[[i]] <- summary(temp)
}
dev.off()
save(tempVpsv_stats_alone, file = './data/stats/tempVpsv_stats_alone.Rdata')
load('./data/stats/tempVpsv_stats_alone.Rdata')

# temp vs mrd
tempVmrd_stats_alone <- vector("list", length = 6)
pdf('./figures/temperature_vs_mrd.pdf')
for (i in 1:6) {
  temp <- lm(values(mrd_raster_list[[i]]) ~ values(temp_list[[i]]))
  plot(values(temp_list[[i]]), values(mrd_raster_list[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C) vs Area", ylab = "mrd")
  abline(lm(values(mrd_raster_list[[i]]) ~ values(temp_list[[i]])), col = 'red')
  tempVmrd_stats_alone[[i]] <- temp
}
dev.off()
save(tempVmrd_stats_alone, file = './data/stats/tempVmrd_stats_alone.Rdata')
load('./data/stats/tempVmrd_stats_alone.Rdata')

# temp vs beta
tempVbeta_stats_alone <- vector("list", length = 6)
pdf('./figures/temperature_vs_beta.pdf')
for (i in 1:6) {
  temp <- lm(values(beta_raster_list[[i]]) ~ values(temp_list[[i]]))
  plot(values(temp_list[[i]]) + area_list[[i]], values(beta_raster_list[[i]]),
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Temperature (°C) vs Area", ylab = "beta")
  abline(lm(values(beta_raster_list[[i]]) ~ values(temp_list[[i]])), col = 'red')
  tempVbeta_stats_alone[[i]] <- temp 
}
dev.off()
save(tempVbeta_stats_alone, file = './data/stats/tempVbeta_stats_alone.Rdata')
load('./data/stats/tempVbeta_stats_alone.Rdata')

# linear regression chlorophyll vs richness
chloro_stats <- vector("list", length = 6)
pdf('./figures/chlorophyll_vs_richness.pdf')
for (i in 1:6) {
     chlor <- lm(values(species_richness[[i]]) ~ values(chloro_list[[i]])) 
     plot(values(chloro_list[[i]]), values(species_richness[[i]]), 
          main = paste('resolution =', res(species_richness[[i]])), 
          xlab = "Chlorophyll (mg/m3)", ylab = "Shark Richness")
     abline(lm(values(species_richness[[i]]) ~ values(chloro_list[[i]])), 
            col = 'red')
     chloro_stats[[i]] <- chlor
}
dev.off()
save(chloro_stats, file = './data/stats/chloro_stats.Rdata')
load('./data/stats/chloro_stats.Rdata')

# linear regression chlorophyll vs psv
chloroVpsv_stats <- vector("list", length = 6)
pdf('./figures/chlorophyll_vs_psv.pdf')
for (i in 1:6) {
  chlor <- lm(values(psv_raster_list[[i]]) ~ values(chloro_list[[i]])) 
  plot(values(chloro_list[[i]]), values(psv_raster_list[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Chlorophyll (mg/m3)", ylab = "PSV")
  abline(lm(values(psv_raster_list[[i]]) ~ values(chloro_list[[i]])), 
         col = 'red')
  chloroVpsv_stats[[i]] <- chlor
}
dev.off()
save(chloroVpsv_stats, file = './data/stats/chloroVpsv_stats.Rdata')
load('./data/stats/chloroVpsv_stats.Rdata')

# linear regression chlorophyll vs mrd
chloroVmrd_stats <- vector("list", length = 6)
pdf('./figures/chlorophyll_vs_mrd.pdf')
for (i in 1:6) {
  chlor <- lm(values(mrd_raster_list[[i]]) ~ values(chloro_list[[i]])) 
  plot(values(chloro_list[[i]]), values(mrd_raster_list[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Chlorophyll (mg/m3)", ylab = "MRD")
  abline(lm(values(mrd_raster_list[[i]]) ~ values(chloro_list[[i]])), 
         col = 'red')
  chloroVmrd_stats[[i]] <- chlor
}
dev.off()
save(chloroVmrd_stats, file = './data/stats/chloroVmrd_stats.Rdata')
load('./data/stats/chloroVmrd_stats.Rdata')

# linear regression chlorophyll vs beta 
chloroVbeta_stats <- vector("list", length = 6)
pdf('./figures/chlorophyll_vs_beta.pdf')
for (i in 1:6) {
  chlor <- lm(values(beta_raster_list[[i]]) ~ values(chloro_list[[i]])) 
  plot(values(chloro_list[[i]]), values(beta_raster_list[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Chlorophyll (mg/m3)", ylab = "Beta")
  abline(lm(values(beta_raster_list[[i]]) ~ values(chloro_list[[i]])), 
         col = 'red')
  chloroVbeta_stats[[i]] <- chlor
}
dev.off()
save(chloroVbeta_stats, file = './data/stats/chloroVbeta_stats.Rdata')
load('./data/stats/chloroVbeta_stats.Rdata')


# linear regression normal vs threatened
iucn_stats <- vector("list", length = 6)
pdf('./figures/IUCN_vs_Normal.pdf')
for (i in 1:6) {
     iucn <- lm(values(iucn_richness[[i]]) ~ values(species_richness[[i]]))
     plot(values(species_richness[[i]]), values(iucn_richness[[i]]), 
          main = paste('resolution =', res(species_richness[[i]])), 
          xlab = "Normal Shark Richness", ylab = "Threatened Shark Richness")
     abline(lm(values(iucn_richness[[i]]) ~ 
            values(species_richness[[i]])), col = 'red')
     iucn_stats[[i]] <- iucn
}
dev.off()
save(iucn_stats, file = './data/stats/iucn_stats.Rdata')
load('./data/stats/iucn_stats.Rdata')

# latitude vs richness
pdf('./figures/latitude_vs_richness.pdf')
for (i in 1:6) {
     plot(abs(latitude_list[[i]]), values(species_richness[[i]]), 
          main = paste('resolution =', res(species_richness[[i]])), xlab = "Latitude", 
          ylab = "Shark Richness")
     abline(lm(values(species_richness[[i]]) ~ abs(latitude_list[[i]])), 
            col = 'red')
}
dev.off()

# salinity vs richness
salinity_stats <- vector("list", length = 6)
pdf('./figures/salinity_vs_richness.pdf')
for (i in 1:6) {
     salinity <- lm(values(species_richness[[i]]) ~ values(salinity_list[[i]]))  
     plot(values(salinity_list[[i]]), values(species_richness[[i]]), 
          main = paste('resolution =', res(species_richness[[i]])), xlab = "Salinity", 
          ylab = "Shark Richness")
     abline(lm(values(species_richness[[i]]) ~ 
                 values(salinity_list[[i]])), col = 'red')
     salinity_stats[[i]] <- salinity
}
dev.off()
save(salinity_stats, file = './data/stats/salinity_stats.Rdata')
load('./data/stats/salinity_stats.Rdata')

# salinity vs psv
salinityVpsv_stats <- vector("list", length = 6)
pdf('./figures/salinity_vs_psv.pdf')
for (i in 1:6) {
  salinity <- lm(values(psv_raster_list[[i]]) ~ values(salinity_list[[i]]))  
  plot(values(salinity_list[[i]]), values(psv_raster_list[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), xlab = "Salinity", 
       ylab = "PSV")
  abline(lm(values(psv_raster_list[[i]]) ~ 
              values(salinity_list[[i]])), col = 'red')
  salinityVpsv_stats[[i]] <- salinity
}
dev.off()
save(salinityVpsv_stats, file = './data/stats/salinityVpsv_stats.Rdata')
load('./data/stats/salinityVpsv_stats.Rdata')

# salinity vs mrd
salinityVmrd_stats <- vector("list", length = 6)
pdf('./figures/salinity_vs_mrd.pdf')
for (i in 1:6) {
  salinity <- lm(values(mrd_raster_list[[i]]) ~ values(salinity_list[[i]]))  
  plot(values(salinity_list[[i]]), values(mrd_raster_list[[i]]), 
       main = paste('resolution =', res(mrd_raster_list[[i]])), xlab = "Salinity", 
       ylab = "MRD")
  abline(lm(values(mrd_raster_list[[i]]) ~ 
              values(salinity_list[[i]])), col = 'red')
  salinityVmrd_stats[[i]] <- salinity
}
dev.off()
save(salinityVmrd_stats, file = './data/stats/salinityVmrd_stats.Rdata')
load('./data/stats/salinityVmrd_stats.Rdata')

# salinity vs beta
salinityVbeta_stats <- vector("list", length = 6)
pdf('./figures/salinity_vs_beta.pdf')
for (i in 1:6) {
  salinity <- lm(values(beta_raster_list[[i]]) ~ values(salinity_list[[i]]))  
  plot(values(salinity_list[[i]]), values(beta_raster_list[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), xlab = "Salinity", 
       ylab = "Beta")
  abline(lm(values(beta_raster_list[[i]]) ~ 
              values(salinity_list[[i]])), col = 'red')
  salinityVbeta_stats[[i]] <- salinity
}
dev.off()
save(salinityVbeta_stats, file = './data/stats/salinityVbeta_stats.Rdata')
load('./data/stats/salinityVbeta_stats.Rdata')

# distance from coast vs richness normal
pdf('./figures/coast_vs_richness_norm.pdf')
for (i in 1:9) {
     plot(values(coast_distance_list[[i]]), values(species_richness_list[[i]]), 
          main = paste('resolution =', res(res_list[[i]])), xlab = "Distance", 
          ylab = "Shark Richness")
}
dev.off()

# distance from coast vs richness log
pdf('./figures/coast_vs_richness_log.pdf')
for (i in 1:9) {
  plot(log(values(coast_distance_list[[i]])), values(species_richness_list[[i]]), 
       main = paste('resolution =', res(res_list[[i]])), xlab = "Distance", 
       ylab = "Shark Richness")
  abline(lm(values(species_richness_list[[i]]) ~ 
              log(values(coast_distance_list[[i]]))), col = 'red')
}
dev.off()

# taxonomic richness vs faith's pd
richnessVpd_stats <- vector("list", length = 6)
pdf('./figures/richness_vs_pd.pdf')
for (i in 1:6) {
     pd <- lm(values(species_richness[[i]]) ~ 
             values(pd_raster_list[[i]]))
     plot(values(pd_raster_list[[i]]), values(species_richness[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Faith's Phylogentic Diversity", 
       ylab = "Taxonomic Richness")
     abline(lm(values(species_richness[[i]]) ~ 
              values(pd_raster_list[[i]])), col = 'red')
     richnessVpd_stats[[i]] <- pd
}
dev.off()
save(richnessVpd_stats, file = './data/stats/richnessVpd_stats.Rdata')
load('./data/stats/richnessVpd_stats.Rdata')

# taxonomic richness vs psv
richnessVpsv_stats <- vector("list", length = 6)
pdf('./figures/richness_vs_psv.pdf')
for (i in 1:6) {
     psv <- lm(values(species_richness[[i]]) ~ 
              values(psv_raster_list[[i]]))
     plot(values(psv_raster_list[[i]]), values(species_richness[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "PSV", 
       ylab = "Taxonomic Richness")
     abline(lm(values(species_richness[[i]]) ~ 
              values(psv_raster_list[[i]])), col = 'red')
     richnessVpsv_stats[[i]] <- psv
}
dev.off()
save(richnessVpsv_stats, file = './data/stats/richnessVpsv_stats.Rdata')
load('./data/stats/richnessVpsv_stats.Rdata')

# taxonomic richness vs psr
richnessVpsr_stats <- vector("list", length = 6)
pdf('./figures/richness_vs_psr.pdf')
for (i in 1:6) {
  psr <- lm(values(species_richness[[i]]) ~ 
              values(psr_raster_list[[i]]))
  plot(values(psr_raster_list[[i]]), values(species_richness[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "PSR", 
       ylab = "Taxonomic Richness")
  abline(lm(values(species_richness[[i]]) ~ 
              values(psr_raster_list[[i]])), col = 'red')
  richnessVpsr_stats[[i]] <- psr
}
dev.off()
save(richnessVpsr_stats, file = './data/stats/richnessVpsr_stats.Rdata')
load('./data/stats/richnessVpsr_stats.Rdata')

# MRD vs psv
mrdVpsv_stats <- vector("list", length = 6)
pdf('./figures/MRD_vs_psv.pdf')
for (i in 1:6) {
  mrd <- lm(values(mrd_raster_list[[i]]) ~ 
              values(psv_raster_list[[i]]))
  plot(values(psv_raster_list[[i]]), values(mrd_raster_list[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Phylogenetic Species Variance (PSV)", 
       ylab = "Mean Root Distance (MRD)")
  abline(lm(values(mrd_raster_list[[i]]) ~ 
              values(psv_raster_list[[i]])), col = 'red')
  mrdVpsv_stats[[i]] <- mrd
}
dev.off()
save(mrdVpsv_stats, file = './data/stats/mrdVpsv_stats.Rdata')
load('./data/stats/mrdVpsv_stats.Rdata')

# MRD vs taxonomic richness
richnessVmrd_stats <- vector("list", length = 6)
pdf('./figures/richness_vs_MRD.pdf')
for (i in 1:6) {
  values(species_richness[[i]])[values(species_richness[[i]]) == 0] <- NA
  values(mrd_raster_list[[i]])[values(mrd_raster_list[[i]]) == 0] <- NA
  rich <- lm(values(species_richness[[i]]) ~ 
               values(mrd_raster_list[[i]]))
  plot(values(mrd_raster_list[[i]]), values(species_richness[[i]]), 
       main = paste('resolution =', res(species_richness[[i]])), 
       xlab = "Mean Root Distance (MRD)", 
       ylab = "Shark Richness")
  abline(lm(values(species_richness[[i]]) ~ 
              values(mrd_raster_list[[i]])), col = 'red')
  richnessVmrd_stats[[i]] <- rich
}
dev.off()
save(richnessVmrd_stats, file = './data/stats/richnessVmrd_stats.Rdata')
load('./data/stats/richnessVmrd_stats.Rdata')



 




library(spatialEco)
library(raster)

# linear regression temp vs richness
richnessVtemp_list <- vector("list", length = 7)
for (i in seq_along(c(1, 2, 3, 4, 5, 6, 7))) {
     plot(values(species_richness_list[[i]]), values(temp_list[[i]]), 
          xlab = "richness", ylab = "temperature (C)")
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
       xlab = "richness", ylab = "chlorophyll (mg/m^3)")
     abline(lm(values(chloro_list[[i]]) ~ values(species_richness_list[[i]])), 
         col = 'red')
     richnessVchloro <- lm(values(chloro_list[[i]]) ~ 
                        values(species_richness_list[[i]]))
     richnessVchloro_list[[i]] <- richnessVchloro
     print(summary(richnessVchloro_list[[i]]))
}

stat_stack <- stack(species_richness_list[[1]], temp_list[[1]])
richnessVtemp <- layerStats(stat_stack, 'pearson', na.rm = T)
richnessVtemp

richnessVtemp <- cor(values(species_richness_list[[1]]), values(temp_list[[1]]), use = "na.or.complete")
richnessVtemp

richnessVtemp <- rasterCorrelation(species_richness_list[[1]], temp_list[[1]], s = 3, type = "spearman")
plot(richnessVtemp)

richnessVchloro <- cor(values(species_richness_list[[1]]), values(chloro), use = "na.or.complete")
richnessVtemp

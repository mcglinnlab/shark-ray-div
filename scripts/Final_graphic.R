library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

# Read in statistics
load('./data/stats/mrdVrichness_stats.Rdata')
load('./data/stats/tempVmrd_stats.Rdata')
load('./data/stats/tempVrichness_stats.Rdata')
load('./data/stats/carmrdVcarrichness_stats.Rdata')
load('./data/stats/tempVcarmrd_stats.Rdata')
load('./data/stats/tempVcarrichness_stats.Rdata')
load('./data/stats/lammrdVlamrichness_stats.Rdata')
load('./data/stats/tempVlammrd_stats.Rdata')
load('./data/stats/tempVlamrichness_stats.Rdata')
load('./data/stats/trop_atlantic_richnessVmrd.Rdata')
load('./data/stats/trop_atlantic_tempVmrd.Rdata')
load('./data/stats/trop_atlantic_richnessVtemp.Rdata')
load('./data/stats/indopacific_richnessVmrd.Rdata')
load('./data/stats/indopacific_tempVmrd.Rdata')
load('./data/stats/indopacific_richnessVtemp.Rdata')
load('./data/stats/sims.Rdata')

sims$bet[1] <- 0.4
sims$bet[2] <- 0.4


# for correlation coefficient
indopacific_richnessVmrd[[1]][[4]]

# create a data table where the columns are the analyses
table_list <- vector("list", length = 6)
for (i in 1:6) {
  cn <- c("Metric", "Niche_Conservatism_(Temperate)", "Niche_Conservatism_(Tropical)", 
          "Ecological_Limits_(Temperate)", "Ecological_Limits_(Tropical)", 
          "Global_Analysis", "Carcharhiniformes", "Lamniformes", 
          "Tropical_Atlantic", "Central_Indo-Pacific")
  rn <- c("Richness vs Temperature", "Richness vs MRD", 
          "Energy Gradient vs MRD", "Beta")
  temp_table <- matrix(ncol = 9, nrow = 4)
  temp_table[1, 1:4] <- sims$richvtemp
  temp_table[2, 1:4] <- sims$richvmrd
  temp_table[3, 1:4] <- sims$energyvmrd
  temp_table[4, 1:4] <- sims$bet
  temp_table[,5] <- c(tempVrichness_stats[[i]][[4]], 
                      mrdVrichness_stats[[i]][[4]], 
                      tempVmrd_stats[[i]][[4]], -0.8868239)
  temp_table[,6] <- c(tempVcarrichness_stats[[i]][[4]], 
                      carmrdVcarrichness_stats[[i]][[4]], 
                      tempVcarmrd_stats[[i]][[4]], -0.9441983)
  temp_table[,7] <- c(tempVlamrichness_stats[[i]][[4]], 
                      lammrdVlamrichness_stats[[i]][[4]], 
                      tempVlammrd_stats[[i]][[4]], -0.6364018)
  temp_table[,8] <- c(trop_atlantic_richnessVtemp[[i]][[4]], 
                      trop_atlantic_richnessVmrd[[i]][[4]], 
                      trop_atlantic_tempVmrd[[i]][[4]], -0.9621481)
  temp_table[,9] <- c(indopacific_richnessVtemp[[i]][[4]], 
                      indopacific_richnessVmrd[[i]][[4]], 
                      indopacific_tempVmrd[[i]][[4]], 
                      -0.8929302)
  temp_table <- data.frame(temp_table)
  temp_table <- cbind(rn, temp_table)
  colnames(temp_table) <- cn
  table_list[[i]] <- temp_table
}

save(table_list, file = './data/stats/table_list.Rdata')

# for error bars
# [3] for correlation coefficient: upper
confint(mod, "x", level=0.95)[2]

# lower
confint(mod, "x", level=0.95)[1]

error_list <- vector("list", length = 6)
for (i in 1:6) {
  error_df <- matrix(ncol = 4, nrow = 20)
  colnames(error_df) <- c("rel", "analysis", "upper_error", "lower_error")
  error_df[,1] <- c(rep("Beta", each = 5), rep("Energy Gradient vs MRD", each = 5), 
                    rep("Richness vs MRD", each = 5), rep("Richness vs Temperature", each = 5))
  error_df[,2] <- rep(c("GA", "Car", "Lam", "Trop_Atl", "Indo"), times = 4)
  error_df[,3] <- c(-0.621985, -0.5738432, 6.558514, -0.5981778, -0.5615923, 
                      confint(tempVmrd_stats[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(tempVcarmrd_stats[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(tempVlammrd_stats[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(trop_atlantic_tempVmrd[[i]][[3]], "x_std", level = 0.95)[2],
                      confint(indopacific_tempVmrd[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(mrdVrichness_stats[[i]][[3]], "x", level = 0.95)[2],
                      confint(carmrdVcarrichness_stats[[i]][[3]], "x", level = 0.95)[2], 
                      confint(lammrdVlamrichness_stats[[i]][[3]], "x", level = 0.95)[2],
                      confint(trop_atlantic_richnessVmrd[[i]][[3]], "x", level = 0.95)[2], 
                      confint(indopacific_richnessVmrd[[i]][[3]], "x", level = 0.95)[2],
                      confint(tempVrichness_stats[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(tempVcarrichness_stats[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(tempVlamrichness_stats[[i]][[3]], "x_std", level = 0.95)[2], 
                      confint(trop_atlantic_richnessVtemp[[i]][[3]], "x_std", level = 0.95)[2],
                      confint(indopacific_richnessVtemp[[i]][[3]], "x_std", level = 0.95)[2])
  
  error_df[,4] <- c(-1.1103004, -1.2345835, -2, -1.2485749, -1.1603735, 
                      confint(tempVmrd_stats[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(tempVcarmrd_stats[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(tempVlammrd_stats[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(trop_atlantic_tempVmrd[[i]][[3]], "x_std", level = 0.95)[1],
                      confint(indopacific_tempVmrd[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(mrdVrichness_stats[[i]][[3]], "x", level = 0.95)[1],
                      confint(carmrdVcarrichness_stats[[i]][[3]], "x", level = 0.95)[1], 
                      confint(lammrdVlamrichness_stats[[i]][[3]], "x", level = 0.95)[1],
                      confint(trop_atlantic_richnessVmrd[[i]][[3]], "x", level = 0.95)[1], 
                      confint(indopacific_richnessVmrd[[i]][[3]], "x", level = 0.95)[1],
                      confint(tempVrichness_stats[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(tempVcarrichness_stats[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(tempVlamrichness_stats[[i]][[3]], "x_std", level = 0.95)[1], 
                      confint(trop_atlantic_richnessVtemp[[i]][[3]], "x_std", level = 0.95)[1],
                      confint(indopacific_richnessVtemp[[i]][[3]], "x_std", level = 0.95)[1])
  error_list[[i]] <- error_df
}
save(error_list, file = './data/stats/error_list.Rdata')

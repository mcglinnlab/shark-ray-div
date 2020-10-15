library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

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
load('./data/stats/sim_values.Rdata')

sim_values$Beta[1] <- 0.4
sim_values$Beta[2] <- 0.4

get_r = function(mod_sum) { 
    r = sqrt(mod_sum$r.squared)
    if (coef(mod_sum)[2, 1] < 0) # if slope is negative
       r = -r
    r
}

# for slope
coef(mrdVrichness_stats[[1]])[2 , 1]

table_list <- vector("list", length = 6)
for (i in 1:6) {
rn <- c("Niche Conservatism (Temperate)", "Niche Conservatism (Tropical)", 
        "Ecological Limits (Temperate)", "Ecological Limits (Tropical)", 
        "Global Analysis", "Carcharhiniformes", "Lamniformes", 
        "Tropical Atlantic", "Central Indo-Pacific")
cn <- c("Analyses", "Richness vs Temperature", "Richness vs MRD", 
        "Energy Gradient vs MRD", "Beta")
temp_table <- matrix(nrow = 9, ncol = 4)
temp_table[1:4,1] <- sim_values$`Richness vs Temperature`
temp_table[1:4,2] <- sim_values$`Richness vs MRD`
temp_table[1:4,3] <- sim_values$`Energy Gradient vs MRD`
temp_table[1:4,4] <- sim_values$Beta
temp_table[5,] <- c(get_r(tempVrichness_stats[[i]]), 
                    coef(mrdVrichness_stats[[i]])[2,1], 
                    get_r(tempVmrd_stats[[i]]), -0.88)
temp_table[6,] <- c(get_r(tempVcarrichness_stats[[i]]), 
                    coef(carmrdVcarrichness_stats[[i]])[2,1], 
                    get_r(tempVcarmrd_stats[[i]]), -0.91)
temp_table[7,] <- c(get_r(tempVlamrichness_stats[[i]]), 
                    coef(lammrdVlamrichness_stats[[i]])[2,1], 
                    get_r(tempVlammrd_stats[[i]]), NA)
temp_table[8,] <- c(get_r(trop_atlantic_richnessVtemp[[i]]), 
                    coef(trop_atlantic_richnessVmrd[[i]])[2,1], 
                    get_r(trop_atlantic_tempVmrd[[i]]), -0.97)
temp_table[9,] <- c(get_r(indopacific_richnessVtemp[[i]]), 
                    coef(indopacific_richnessVmrd[[i]])[2,1], 
                    get_r(indopacific_tempVmrd[[i]]), 
                    -0.79)
temp_table <- data.frame(temp_table)
temp_table <- cbind(rn, temp_table)
colnames(temp_table) <- cn
table_list[[i]] <- temp_table
}

save(table_list, file = './data/stats/table_list.Rdata')

pdf('./figures/gradient_table.pdf', width = 7 * 1.25)
for (i in 1:6) {
    tableau.m <- reshape::melt(table_list[[i]])
    tableau.m$value <- ifelse(is.na(tableau.m$value), 0, tableau.m$value)
    tableau.m = tableau.m %>% 
        group_by(variable) %>% 
        mutate(rescale = value)

tableau.m$Analyses = factor(tableau.m$Analyses,
                            levels = rev(unique(tableau.m$Analyses)))


p <- ggplot(tableau.m, aes(variable, Analyses)) + 
    geom_tile(aes(fill = rescale), colour = "white") + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
p

base_size <- 9
p + theme_grey(base_size = base_size) + 
  labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none", axis.ticks = element_blank(), 
        axis.text.x = element_text(size = base_size * 0.5, angle = 0, 
                                   hjust = 0, colour = "gray50"))
print(p)
}
dev.off()

# testing different methods
tes <- table_list[[4]]
colnames(tes) <- c("Analyses", "Richness_vs_Temperature", "Richness_vs_MRD", 
                   "Energy_Gradient_vs_MRD","Beta")

tes1 <- tes[1, 2:5]
tes2 <- tes[5, 2:5]
tesr <- data.frame(tes1, tes2)

d1 <- ggplot() + geom_col(data = tes, aes(x = Analyses, y = Richness_vs_Temperature, fill = Analyses)) +
  scale_x_discrete(limits = c("Niche Conservatism (Temperate)", "Niche Conservatism (Tropical)", 
                              "Ecological Limits (Temperate)", "Ecological Limits (Tropical)", 
                              "Global Analysis", "Carcharhiniformes", "Lamniformes", 
                              "Tropical Atlantic", "Central Indo-Pacific")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

d2 <- ggplot() + geom_col(data = tes, aes(x = Analyses, y = Richness_vs_MRD, fill = Analyses)) +
  scale_x_discrete(limits = c("Niche Conservatism (Temperate)", "Niche Conservatism (Tropical)", 
                              "Ecological Limits (Temperate)", "Ecological Limits (Tropical)", 
                              "Global Analysis", "Carcharhiniformes", "Lamniformes", 
                              "Tropical Atlantic", "Central Indo-Pacific")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +

d3 <- ggplot() + geom_col(data = tes, aes(x = Analyses, y = Energy_Gradient_vs_MRD, fill = Analyses)) +
  scale_x_discrete(limits = c("Niche Conservatism (Temperate)", "Niche Conservatism (Tropical)", 
                              "Ecological Limits (Temperate)", "Ecological Limits (Tropical)", 
                              "Global Analysis", "Carcharhiniformes", "Lamniformes", 
                              "Tropical Atlantic", "Central Indo-Pacific")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

d4 <- ggplot() + geom_col(data = tes, aes(x = Analyses, y = Beta, fill = Analyses)) +
  scale_x_discrete(limits = c("Niche Conservatism (Temperate)", "Niche Conservatism (Tropical)", 
                              "Ecological Limits (Temperate)", "Ecological Limits (Tropical)", 
                              "Global Analysis", "Carcharhiniformes", "Lamniformes", 
                              "Tropical Atlantic", "Central Indo-Pacific")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf('./figures/test_graphic.pdf')
grid.arrange(d1, d2, d3, d4)
dev.off()
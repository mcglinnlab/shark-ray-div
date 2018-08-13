library(dplyr)
library(ggplot2)
library(scales)

load('./data/stats/richnessVmrd_stats.Rdata')
load('./data/stats/tempVmrd_stats.Rdata')
load('./data/stats/carmrdVcarrichness_stats.Rdata')
load('./data/stats/tempVcarmrd_stats.Rdata')
load('./data/stats/lammrdVlamrichness_stats.Rdata')
load('./data/stats/tempVlammrd_stats.Rdata')
load('./data/stats/trop_atlantic_richnessVmrd.Rdata')
load('./data/stats/trop_atlantic_tempVmrd.Rdata')
load('./data/stats/indopacific_richnessVmrd.Rdata')
load('./data/stats/indopacific_tempVmrd.Rdata')


table_list <- vector("list", length = 6)
for (i in 1:6) {
rn <- c("Niche Conservatism", "Ecological Limits (Temperate)", "Ecological Limits (Tropical)", 
        "Global Analysis", "Carcharhiniformes", "Lamniformes", "Tropical Atlantic", "Central Indo-Pacific")
cn <- c("Analyses", "Richness vs MRD", "Energy Gradient vs MRD", "Beta")
temp_table <- matrix(nrow = 8, ncol = 3)
temp_table[1,] <- c(0, 0, 1)
temp_table[2,] <- c(1, 1, -1)
temp_table[3,] <- c(-1, -1, -1)
temp_table[4,] <- c(richnessVmrd_stats[[i]]$coefficients[2], tempVmrd_stats[[i]]$coefficients[2], -0.88)
temp_table[5,] <- c(carmrdVcarrichness_stats[[i]]$coefficients[2], 
                    tempVcarmrd_stats[[i]]$coefficients[2], -0.91)
temp_table[6,] <- c(lammrdVlamrichness_stats[[i]]$coefficients[2], tempVlammrd_stats[[i]]$coefficients[2], NA)
temp_table[7,] <- c(trop_atlantic_richnessVmrd[[i]]$coefficients[2], 
                    trop_atlantic_tempVmrd[[i]]$coefficients[2], -0.97)
temp_table[8,] <- c(indopacific_richnessVmrd[[i]]$coefficients[2], indopacific_tempVmrd[[i]]$coefficients[2], 
                    -0.79)
temp_table <- data.frame(temp_table)
temp_table <- cbind(rn, temp_table)
colnames(temp_table) <- cn
table_list[[i]] <- temp_table
}

save(table_list, file = './data/stats/table_list.Rdata')

pdf('./figures/gradient_table.pdf')
for (i in 1:6) {
tableau.m <- reshape::melt(table_list[[i]])
tableau.m$value = as.numeric(as.character(tableau.m$value))
tableau.m = tableau.m %>% 
  group_by(variable) %>% 
  mutate(rescale = rescale(value))

(p <- ggplot(tableau.m, aes(variable, Analyses)) + 
    geom_tile(aes(fill = rescale), colour = "white") + 
    scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0))

base_size <- 9
p + theme_grey(base_size = base_size) + 
  labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none", axis.ticks = element_blank(), 
        axis.text.x = element_text(size = base_size * 0.5, angle = 0, 
                                   hjust = 0, colour = "gray50"))
}
dev.off()


tableau <- read.table(
  text = 
    "Net    B   C   D   E.(e)   F.(f)
a   1.88    0.15    0.60    10.00   90.00
b   2.05    0.23    0.51    55.00   80.00
c   2.09    0.29    0.40    58.00   88.00
d   2.07    0.52    0.36    80.00   84.00
e   2.13    0.30    0.27    7.00    90.00",
  header = TRUE)

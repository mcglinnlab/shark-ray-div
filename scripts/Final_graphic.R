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
load('./data/stats/sim_values.Rdata')

sim_values$Beta[1] <- 0.4
sim_values$Beta[2] <- 0.4

# function to get r squared value
get_r = function(mod_sum) { 
  r = sqrt(mod_sum$r.squared)
  if (coef(mod_sum)[2, 1] < 0) # if slope is negative
    r = -r
  r
}

# for slope
coef(mrdVrichness_stats[[1]])[2 , 1]

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
  temp_table[1, 1:4] <- sim_values$`Richness vs Temperature`
  temp_table[2, 1:4] <- sim_values$`Richness vs MRD`
  temp_table[3, 1:4] <- sim_values$`Energy Gradient vs MRD`
  temp_table[4, 1:4] <- sim_values$Beta
  temp_table[,5] <- c(get_r(tempVrichness_stats[[i]]), 
                      coef(mrdVrichness_stats[[i]])[2,1], 
                      get_r(tempVmrd_stats[[i]]), -0.88)
  temp_table[,6] <- c(get_r(tempVcarrichness_stats[[i]]), 
                      coef(carmrdVcarrichness_stats[[i]])[2,1], 
                      get_r(tempVcarmrd_stats[[i]]), -0.91)
  temp_table[,7] <- c(get_r(tempVlamrichness_stats[[i]]), 
                      coef(lammrdVlamrichness_stats[[i]])[2,1], 
                      get_r(tempVlammrd_stats[[i]]), NA)
  temp_table[,8] <- c(get_r(trop_atlantic_richnessVtemp[[i]]), 
                      coef(trop_atlantic_richnessVmrd[[i]])[2,1], 
                      get_r(trop_atlantic_tempVmrd[[i]]), -0.97)
  temp_table[,9] <- c(get_r(indopacific_richnessVtemp[[i]]), 
                      coef(indopacific_richnessVmrd[[i]])[2,1], 
                      get_r(indopacific_tempVmrd[[i]]), 
                      -0.79)
  temp_table <- data.frame(temp_table)
  temp_table <- cbind(rn, temp_table)
  colnames(temp_table) <- cn
  table_list[[i]] <- temp_table
}

save(table_list, file = './data/stats/table_list.Rdata')

# function to get goodness of fit value
# where pred is the expected model
# where obs is the observed result
gof <- function(pred, obs) {
  pred_sub <- substitute(pred)
  obs_sub <- substitute(obs)
res_result <- ((table_list[[4]][eval(y_sub)] - 
                  table_list[[4]][eval(x_sub)])^2)
res_result <- sum(res_result, na.rm = T)
print(res_result)
}

# Function to get r squared between model and data. Function has been checked by hand
gof <- function(pred, obs, na.rm = F) {
  pred_sub <- substitute(pred)
  obs_sub <- substitute(obs)
  observed <- table_list[[4]][eval(obs_sub)]
  predicted <- table_list[[4]][eval(pred_sub)]
    if (na.rm) {
      true = !(is.na(observed) | is.na(predicted))
      observed = observed[true]
      predicted = predicted[true]
      observed = data.frame(observed)
    }
  SSerr = sum((observed - predicted)^2)
  SStot = sum((observed - mean(observed[,1]))^2)
  R2 = 1 - SSerr / SStot
  return(R2)
}

nctempXg <- gof('Niche_Conservatism_(Temperate)', 'Global_Analysis')
nctempXc <- gof('Niche_Conservatism_(Temperate)', 'Carcharhiniformes')
nctempXl <- gof('Niche_Conservatism_(Temperate)', 'Lamniformes', na.rm = T)
nctempXta <- gof('Niche_Conservatism_(Temperate)', 'Tropical_Atlantic')
nctempXip <- gof('Niche_Conservatism_(Temperate)', 'Central_Indo-Pacific')

nctropXg <- gof('Niche_Conservatism_(Tropical)', 'Global_Analysis')
nctropXc <- gof('Niche_Conservatism_(Tropical)', 'Carcharhiniformes')
nctropXl <- gof('Niche_Conservatism_(Tropical)', 'Lamniformes', na.rm = T)
nctropXta <- gof('Niche_Conservatism_(Tropical)', 'Tropical_Atlantic')
nctropXip <- gof('Niche_Conservatism_(Tropical)', 'Central_Indo-Pacific')

eltempXg <- gof('Ecological_Limits_(Temperate)', 'Global_Analysis')
eltempXc <- gof('Ecological_Limits_(Temperate)', 'Carcharhiniformes')
eltempXl <- gof('Ecological_Limits_(Temperate)', 'Lamniformes', na.rm = T)
eltempXta <- gof('Ecological_Limits_(Temperate)', 'Tropical_Atlantic')
eltempXip <- gof('Ecological_Limits_(Temperate)', 'Central_Indo-Pacific')

eltropXg <- gof('Ecological_Limits_(Tropical)', 'Global_Analysis')
eltropXc <- gof('Ecological_Limits_(Tropical)', 'Carcharhiniformes')
eltropXl <- gof('Ecological_Limits_(Tropical)', 'Lamniformes', na.rm = T)
eltropXta <- gof('Ecological_Limits_(Tropical)', 'Tropical_Atlantic')
eltropXip <- gof('Ecological_Limits_(Tropical)', 'Central_Indo-Pacific')


analysis <- c("Global Analysis", "Carcharhiniformes", "Lamniformes", 
              "Tropical Atlantic", "Central Indo-Pacific")
analysis <- rep(analysis, times = 4)
gof_value <- c(nctempXg, nctempXc, nctempXl, nctempXta, nctempXip, 
               nctropXg, nctropXc, nctropXl, nctropXta, nctropXip,
               eltempXg, eltempXc, eltempXl, eltempXta, eltempXip,
               eltropXg, eltropXc, eltropXl, eltropXta, eltropXip)


model <- c(rep("Niche Conservatism (Temperate)", times = 5), 
           rep("Niche Conservatism (Tropical)", times = 5),
           rep("Ecological limits (Temperate)", times = 5),
           rep("Ecological limits (Tropical)", times = 5))
new_df <- data.frame(analysis, gof_value, model)
colnames(new_df) <- c("analysis", "gof_value", "model")
new_df$analysis <- factor(new_df$analysis, levels=unique(new_df$analysis))

# untransformed figure
pdf('./figures/final_figure_no_transformation.pdf', width = 7*1.5)
ggplot() + geom_col(data = new_df, aes(x = model, y = gof_value, fill = model)) +
  facet_wrap(~analysis, nrow=1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_discrete(name = "Model", labels = c("Ecological Limits (Temperate)", 
                                                 "Ecological Limits (Tropical)",
                                                 "Niche Conservatism (Temperate)",
                                                 "Niche Conservatism (Tropical)")) +
  xlab("Model") +
  ylab("R Squared Value")
dev.off()

# transformed figure
new_df_trans <- data.frame(analysis, abs(gof_value), model)
colnames(new_df_trans) <- c("analysis", "gof_value", "model")
new_df_trans$analysis <- factor(new_df_trans$analysis, levels=unique(new_df_trans$analysis))

pdf('./figures/final_graphic_transformed.pdf')
ggplot() + geom_col(data = new_df_trans, aes(x = model, y = gof_value, fill = model)) +
  facet_wrap(~analysis, nrow=1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_discrete(name = "Model", labels = c("Ecological Limits (Temperate)", 
                                                 "Ecological Limits (Tropical)",
                                                 "Niche Conservatism (Temperate)",
                                                 "Niche Conservatism (Tropical)")) +
  xlab("Model") +
  ylab("Goodness of Fit Value")
dev.off()





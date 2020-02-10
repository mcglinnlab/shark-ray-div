library(tidyr)
library(ggplot2)
library(purrr)

load('./figures/table_list_for_plotting.Rdata')
ls()
dim(table_list[[1]])

out <- do.call("rbind", table_list)
out <- data.frame(scale = rep(1:6, each = 4), out)

# fix column names
names(out)
names(out) = c('scale', 'rel', 'NC_temp', 'NC_trop', 'EC_temp', 
               'EC_trop', 'GA', 'Car', 'Lam', 'Trop_Atl', 'Indo')

out %>%
    ggplot(aes(GA, EC_temp)) + 
    geom_point(aes(color = scale, pch = metric)) +
    geom_abline(intercept = 0, slope = 1)

out %>%
    ggplot(aes(GA, EC_trop)) + 
    geom_point(aes(color = scale, pch = metric)) +
    geom_abline(intercept = 0, slope = 1)

out %>%
    ggplot(aes(GA, NC_trop)) + 
    geom_point(aes(color = scale, pch = metric)) +
    geom_abline(intercept = 0, slope = 1)

# restructure from wide to long format %>%
pred_long <- subset(out, select = -c(GA, Car, Lam, Trop_Atl, Indo)) %>%
    pivot_longer(-c(scale, rel), names_to = "model", values_to = "pred")

obs_long <- subset(out, select = -c(NC_temp, NC_trop, EC_temp, EC_trop)) %>%
    pivot_longer(-c(scale, rel), names_to = "analysis", values_to = "obs")

# merge back observed metrics into long format expected matrix

out_long <- merge(pred_long, obs_long)

# now make plots

out_long %>%
    subset(subset = scale == 3) %>%
    ggplot(aes(pred, obs)) +
    geom_point(aes(col = model, pch = rel, cex = 2)) + 
    geom_abline(intercept = 0, slope = 1) +
    ylim(-1, 1) + 
    xlim(-1, 1) + 
    facet_wrap(. ~ analysis)

ggsave('./figures/pred_obs_metrics.pdf', width = 7 * 2, height = 7 * 1)


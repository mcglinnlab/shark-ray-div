library(tidyr)
library(ggplot2)
library(purrr)
library(dplyr)
library(wesanderson)

load('./data/stats/table_list.Rdata')
load('./data/stats/error_list.Rdata')
dim(table_list[[1]])

out <- do.call("rbind", table_list)
out <- data.frame(scale = rep(1:6, each = 4), out)
error_out <- do.call("rbind", error_list)
error_out <- data.frame(scale = rep(1:6, each = 20), error_out)

# fix column names
names(out)
names(out) = c('scale', 'rel', 'NC_temp', 'NC_trop', 'EC_temp', 
               'EC_trop', 'GA', 'Car', 'Lam', 'Trop_Atl', 'Indo')

out %>%
    ggplot(aes(GA, EC_temp)) + 
    geom_point(aes(color = scale, pch = rel)) +
    geom_abline(intercept = 0, slope = 1)

out %>%
    ggplot(aes(GA, EC_trop)) + 
    geom_point(aes(color = scale, pch = rel)) +
    geom_abline(intercept = 0, slope = 1)

out %>%
    ggplot(aes(GA, NC_trop)) + 
    geom_point(aes(color = scale, pch = rel)) +
    geom_abline(intercept = 0, slope = 1)

# restructure from wide to long format %>%
pred_long <- subset(out, select = -c(GA, Car, Lam, Trop_Atl, Indo)) %>%
    pivot_longer(-c(scale, rel), names_to = "model", values_to = "pred")

obs_long <- subset(out, select = -c(NC_temp, NC_trop, EC_temp, EC_trop)) %>%
    pivot_longer(-c(scale, rel), names_to = "analysis", values_to = "obs")

# merge error with observed metrics
obs_long <- merge(obs_long, error_out)

# merge back observed metrics into long format expected matrix

out_long <- merge(pred_long, obs_long)

# compute R2 of 1:1 lines between pred and obs

get_SSerr = function(obs, pred, na.rm=TRUE) {
  if (na.rm) {
    true = !(is.na(obs) | is.na(pred))
    obs = obs[true]
    pred = pred[true]
  }
  SSerr = mean((obs - pred)^2)
  return(SSerr)
}

SSerr <- out_long %>%
    group_by(scale, analysis, model) %>%
    summarize(SSerr = get_SSerr(obs, pred))

SSerr$sSSerr <- SSerr$SSerr / 2
SSerr
write.csv(SSerr, file = './data/SSerr.csv', row.names = F) 
# now make plots
SSerr <- read.csv('./data/SSerr.csv')

# fix this
out_long$upper_error <- as.numeric(as.character(out_long$upper_error))
out_long$lower_error <- as.numeric(as.character(out_long$lower_error))

out_long2 <- out_long
out_long2$analysis <- gsub(pattern = "GA", replacement = "Global Analysis", 
                           x = out_long2$analysis)
out_long2$analysis <- gsub(pattern = "Car", replacement = "Requiem Sharks", 
                           x = out_long2$analysis)
out_long2$analysis <- gsub(pattern = "Lam", replacement = "Mackerel Sharks", 
                           x = out_long2$analysis)
out_long2$analysis <- gsub(pattern = "Trop_Atl", replacement = "Tropical Atlantic", 
                           x = out_long2$analysis)
out_long2$analysis <- gsub(pattern = "Indo", replacement = "Central Indo-Pacific", 
                           x = out_long2$analysis)
out_long2$model <- gsub(pattern = "EC_temp", replacement = "ELH (temperate origin)", 
                           x = out_long2$model)
out_long2$model <- gsub(pattern = "EC_trop", replacement = "ELH (tropical origin)", 
                        x = out_long2$model)
out_long2$model <- gsub(pattern = "NC_temp", replacement = "NCH (temperate origin)", 
                        x = out_long2$model)
out_long2$model <- gsub(pattern = "NC_trop", replacement = "NCH (tropical origin)", 
                        x = out_long2$model)
neworder <- c("Global Analysis", "Requiem Sharks", "Mackerel Sharks", 
              "Tropical Atlantic", "Central Indo-Pacific")
out_long2 <- arrange(transform(out_long2,
                               analysis=factor(analysis,levels=neworder)), analysis)

# obs pred for all
out_long2 %>%
    subset(subset = scale == 4) %>%
    ggplot(aes(pred, obs)) +
    geom_point(aes(col = model, pch = rel, cex = 2), alpha = 0.9) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_errorbar(aes(ymin=lower_error, ymax=upper_error), width=.05, alpha = 0.5) +
    labs(col = "Model", pch = "Metric") +
    scale_color_manual(values = rev(wes_palette("GrandBudapest2"))) +
    ylim(-2, 6.6) + 
    xlim(-2, 2) + 
    xlab("Predicted Values") +
    ylab("Observed Values") +
    theme(strip.text.x = element_text(size = 12), axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10), legend.text = element_text(size = 9),
        legend.title = element_text(size = 10), aspect.ratio = 1.1) +
    facet_wrap(. ~ analysis, scales = "free", nrow = 1)

ggsave('./figures/pred_obs_metrics.png', width = 7 * 2, height = 7 * 2)


# jitter x-position
tmp <- out_long2
rand <- runif(nrow(tmp), min = -.05, max = 0.05)
tmp$pred <- tmp$pred + rand
# split the variable model up into hypo and orig
tmp$hypo <- ifelse(grepl('ELH', tmp$model), 'ELH', 'NCH')    
tmp$orig <- ifelse(grepl('temperate', tmp$model), 'temperate', 'tropical')
# fix labels on the variable rel
tmp$rel <- gsub('Energy Gradient vs MRD', 'temperature vs MRD', tmp$rel)
tmp$rel <- gsub('Richness vs MRD', 'richness vs MRD', tmp$rel)
tmp$rel <- gsub('Richness vs Temperature', 'temperature vs richness', tmp$rel)
# fix labels on the hypo
hypo_labs <- c("Energy Limits Hypothesis", "Niche Conservatism Hypothesis") 
names(hypo_labs) <- c("ELH", "NCH")

tmp %>%
  subset(subset = scale == 4) %>%
  filter(analysis == "Global Analysis") %>%
  ggplot(aes(pred, obs)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 0, lty = 2, col = 'grey') + 
  geom_vline(xintercept = 0, lty = 2, col = 'grey') +     
  geom_errorbar(aes(ymin=lower_error, ymax=upper_error), width=.05, alpha = 0.5) +
  geom_point(aes(col = orig, pch = rel), cex = 2, stroke = 2) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(col = "Center of Origin", pch = "Metric") +
  scale_color_manual(values = c('blue', 'red')) +
  ylim(-1.5, 1.5) + 
  xlim(-1.5, 1.5) + 
  xlab("Predicted Values") +
  ylab("Observed Values") +
  theme_classic() + 
  theme(strip.background = element_blank(), strip.text = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), aspect.ratio = 1.1, 
        panel.spacing = unit(1, "lines")) +
  facet_wrap(. ~ hypo, labeller = labeller(hypo = hypo_labs),
             scales = "fixed", nrow = 1)
ggsave('./figures/pred_obs_metrics_global_sepearted.png')



# obs pred for remainder without global
out_long2 %>%
  subset(subset = scale == 4) %>%
  filter(analysis != "Global Analysis") %>%
  ggplot(aes(pred, obs)) +
  geom_point(aes(col = model, pch = rel, cex = 6), alpha = 0.9) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_errorbar(aes(ymin=lower_error, ymax=upper_error), width=.05, alpha = 0.5) +
  scale_color_manual(values = rev(wes_palette("GrandBudapest2"))) +
  labs(col = "Model", pch = "Metric") +
  xlab("Predicted Values") +
  ylab("Observed Values") +
  theme(strip.text.x = element_text(size = 12), axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), legend.text = element_text(size = 9),
        legend.title = element_text(size = 10), aspect.ratio = 1.1) +
  facet_wrap(. ~ analysis, scales = "free", nrow = 1)

ggsave('./figures/pred_obs_metrics_rest.pdf')

# Plot for final bar graph
SSerr_wanted_scale <- filter(SSerr, scale == 4)
SSerr_wanted_scale$analysis <- gsub(pattern = "GA", replacement = "Global Analysis", 
                           x = SSerr_wanted_scale$analysis)
SSerr_wanted_scale$analysis <- gsub(pattern = "Car", replacement = "Requiem Sharks", 
                           x = SSerr_wanted_scale$analysis)
SSerr_wanted_scale$analysis <- gsub(pattern = "Lam", replacement = "Mackerel Sharks", 
                           x = SSerr_wanted_scale$analysis)
SSerr_wanted_scale$analysis <- gsub(pattern = "Trop_Atl", replacement = "Tropical Atlantic", 
                           x = SSerr_wanted_scale$analysis)
SSerr_wanted_scale$analysis <- gsub(pattern = "Indo", replacement = "Central Indo-Pacific", 
                           x = SSerr_wanted_scale$analysis)
neworder <- c("Global Analysis", "Requiem Sharks", "Mackerel Sharks", 
              "Tropical Atlantic", "Central Indo-Pacific")
SSerr_wanted_scale <- arrange(transform(SSerr_wanted_scale,
                                        analysis=factor(analysis,levels=neworder)), analysis)

# bar graph all
pdf('./figures/final_figure_bar_plot.pdf', width = 7*1.5)
SSerr_wanted_scale %>% 
  ggplot(aes(x = model, y = SSerr, fill = model)) +
  geom_col() +
  facet_wrap(~analysis, nrow=1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size = 15), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_fill_manual(values = rev(wes_palette("GrandBudapest2")), 
                    name = "Model", labels = c("ELH (Temperate)","ELH (Tropical)",
                                               "NCH (Temperate)","NCH (Tropical)")) +
  xlab("Model") +
  ylab("Mean Sum Squared Error (MSE)")
dev.off()

# bar graph global only
pdf('./figures/final_figure_bar_plot_global.pdf', width = 7*1.5)
SSerr_wanted_scale %>% 
  filter(analysis == "Global Analysis") %>%
  ggplot(aes(x = model, y = SSerr, fill = model)) +
  geom_col() +
  facet_wrap(~analysis) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size = 30), axis.title.x = element_text(size = 19),
        axis.title.y = element_text(size = 19), legend.text = element_text(size = 18),
        legend.title = element_text(size = 19)) +
  scale_fill_manual(values = rev(wes_palette("GrandBudapest2")), 
                    name = "Model", labels = c("ELH (Temperate)","ELH (Tropical)",
                                              "NCH (Temperate)","NCH (Tropical)")) +
  xlab("Model") +
  ylab("Mean Sum of Squares Error (MSE)")
dev.off()

# bar graph for rest
pdf('./figures/final_figure_bar_plot_rest.pdf', width = 7*1.5)
SSerr_wanted_scale %>% 
  filter(analysis != "Global Analysis") %>%
  ggplot(aes(x = model, y = SSerr, fill = model)) +
  geom_col() +
  facet_wrap(~analysis, nrow=1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
  scale_fill_manual(values = rev(wes_palette("GrandBudapest2")), 
                    name = "Model", labels = c("ELH (Temperate)","ELH (Tropical)",
                                               "NCH (Temperate)","NCH (Tropical)")) +
  xlab("Model") +
  ylab("Mean Sum of Squares Error (MSE)")
dev.off()

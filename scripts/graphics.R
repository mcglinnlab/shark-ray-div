library(ggplot2)
library(BiocManager)
library(ggtree)
library(egg)
library(apTreeshape)

# This script creates graphics for the simulation results and expectations given each model.
# the graphics produced here comprise Figure 1. Correlation coefficients for 
# richness vs temperature, richness vs MRD, and MRD vs temperature, and the beta
# metric of tree imbalance are computed here for every hypothesis. These are to be called
# later when computing MSE

#!/usr/bin/env Rscript

#sim = commandArgs();
#sim = as.numeric(sim[length(sim)]);

# Choose number of time slices per simulation to analyze
num.of.time.slices = 1; # use -999 if you want to define specific time slices
# which.time.slices is (apparently) only for specifying particular, unevenly spaced time slices;
# if not being used it should be set to -999
which.time.slices = 30000;
# time.sequence is (apparently) for when the time slices occur for a regular interval; set to -999 if not being used
# Note that due to the slow calculation of tree imbalance (beta) for large trees, it may be best to specify only ~20 time slices
time.sequence = -999
#time.sequence = seq(2,300,by=2); # for time scenario sims
#time.sequence = seq(1000,100000,length=100); # for energy gradient sims

# choose root only or all clades
root.only = 1 # 0 means all clades, 1 means just the root

# Set minimum number of species in a clade needed to proceed with analysis
min.num.spp = 8;

partial.analysis = 1; # toggle to determine whether we're looking at all sims or just some
setwd("/home/sheahane/shark-ray-div")
sim_dir = "./data/raw_sim_output"
analysis_dir = "./data/analysis_output"
already_unzipped = 1
# Simulation workflow

#(2) load simulation and analysis functions
library(ape);
library(permute);
library(nlme);
library(vegan);
library(picante);
library(mvtnorm);
library(caper);
library(paleotree);
library(plyr);
library(phytools);
library(foreach);
library(doParallel);
library(apTreeshape);

package.vector = c('ape','permute','nlme','vegan','picante','mvtnorm','caper','paleotree','plyr','phytools','apTreeshape');

source('./scripts/reg_calc_and_analysis.r');
source('./scripts/make.phylo.jimmy.fun.r');
source('./scripts/lat.grad.time.plot.r');
source('./scripts/clade.origin.corr.plot.r');
source('./scripts/clade.exmpl.figs.r');
source('./scripts/extinct.calc.r');
source('./scripts/unzipping_files.r');

cl = makeCluster(20);
registerDoParallel(cl);

#(3) read in master simulation matrix with chosen parameter combinations;
# then add fields for storing output summary
sim.matrix = as.data.frame(read.csv("./data/SENC_Master_Simulation_Matrix.csv",header=T));
sim.matrix$n.regions = NA
sim.matrix$extant.S = NA
sim.matrix$extinct.S = NA
sim.matrix$skipped.clades = NA
sim.matrix$skipped.times = NA
#sim.matrix$BK.reg = NA
#sim.matrix$BK.env = NA

# Our analysis only concerns the energy limits and niche conservatism under different
# temperature origin scenarios
partial.analysis <- 1

# This function must run for each set of sim ID's representing each model
# function to return reg.summary2 for each needed sim
# wanted_sim = the numeric sim id
regsum_fun <- function(wanted_sim) {
  sim <- wanted_sim
  .packages <- package.vector
  .combine <- 'rbind'

  # (5) read in simulation results for specified simulation from the output zip file 
  all.populations = read.csv(paste(sim_dir,'/sim',sim,'_out/SENC_all.pops_sim',sim,'.csv',sep=''), header=T)
  time.richness = read.csv(paste(sim_dir,'/sim', sim, '_out/SENC_time.rich_sim',sim,'.csv',sep=''), header=T)
  phylo.out = read.tree(paste(sim_dir,'/sim', sim, '_out/SENC_phylo_sim',sim,'.tre',sep=''))
  params.out = read.csv(paste(sim_dir,'/sim', sim, '_out/SENC_params.out_sim',sim,'.csv',sep=''), header=T)

  max.time.actual = max(time.richness$time);
  # If just a single timeslice, then use the end of the simulation or a designated time, otherwise space them equally (which.time.slices == -999)
  # or use specified vector in which.time.slices
  if (num.of.time.slices == 1) { 
    timeslices = max.time.actual 
  } else {
   if (which.time.slices != -999 & num.of.time.slices == - 999) { timeslices = which.time.slices };
   if (which.time.slices == -999 & num.of.time.slices > 1) {timeslices = as.integer(round(seq(max(time.richness$time)/num.of.time.slices,max(time.richness$time),length=num.of.time.slices),digits=0))};
   if (time.sequence[1] != -999) {timeslices = subset(time.sequence, time.sequence <= max.time.actual)};
  }

  skipped.clades = 0
  skipped.times = ""
  for (t in timeslices) {
    # vector of species in existence at time t
    sub.species = as.character(unique(subset(all.populations,time.of.sp.origin <= t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);
  
    # Some species may be extant globally (extant==1) but in our boundary regions (0,11) only;
    # we need to eliminate species that are not extant within regions 1-10 (which is all that is
    # reflected in the all.populations dataframe)
    time.slice.populations = all.populations;
    time.slice.populations$extant = 0;
    time.slice.populations$extant[time.slice.populations$time.of.origin <= t & time.slice.populations$time.of.extinction > t] = 1
    extant.ornot = aggregate(time.slice.populations$extant,by=list(time.slice.populations$spp.name),sum)
    extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
  
    # FIXME:
    # Add more explanatory comments justifying why we don't need to consider species that existed
    # at time t but went extinct before the present.
    # In some cases (e.g. sim 1 or 2, t=6000), tips.to.drop includes all tips and so sub.phylo is empty.
    # Does it make sense for this to ever happen? If not, fix it.
    # If so, need to provide an if-else error catch both in the creation of sub.phylo,
    # and of sub.clade.phylo inside the clade loop. (Sim 3, t = 156 bonks at that point)
    # NOTE: code runs for sim==5 currently as a test case
    sub.species2 = sub.species[!sub.species %in% extinct.species]
    tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
  
    # check to see if there are at least min.num.spp species for continuing with the analysis; if not store the skipped timeslice
    if ( (length(phylo.out$tip.label) - length(tips.to.drop)) < min.num.spp) {
      skipped.times = paste(skipped.times, t) # keep track of the timeslices that were skipped in a text string
    } else {
    
      sub.phylo = drop.tip(phylo.out,tips.to.drop);
      temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1]); temp.root.time;
      most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1])))]; most.recent.spp;
      extinct.time.most.recent = unique(all.populations$time.of.sp.extinction[all.populations$spp.name==most.recent.spp]); extinct.time.most.recent;
      sub.phylo$root.time = temp.root.time + max(c(0,max.time.actual - extinct.time.most.recent)); sub.phylo$root.time;
      sub.phylo = collapse.singles(timeSliceTree(sub.phylo,sliceTime=(max.time.actual - t),plot=F,drop.extinct = T));
      num.of.spp = length(sub.phylo$tip.label);
    
      if (root.only == 1) { sub.clade.loop.end = (num.of.spp+1) }
      if (root.only == 0) { sub.clade.loop.end = max(sub.phylo$edge) }
      for (c in (num.of.spp+1):sub.clade.loop.end) {
      
        #pull out list of species names belonging to each subclade
        sub.clade = clade.members(c, sub.phylo, tip.labels=T)
        subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
      
        #sub.populations is the subset of populations specific to a particular clade and timeslice
        sub.populations = subset(subset.populations, time.of.origin <= t & time.of.extinction > t)
      
        #sub.clade.phylo is a specific simulation clade pulled from the phylogeny that was sliced at timeslice t
        tips.to.drop2 = as.character(sub.phylo$tip.label[which(is.element(sub.phylo$tip.label,as.character(sub.populations$spp.name))==F)]);
      
        # check to see if there are at least min.num.spp species for continuing with the analysis; if not increment skipped.clades
        if((length(sub.phylo$tip.label) - length(tips.to.drop2)) < min.num.spp) {
          skipped.clades = skipped.clades + 1
        } else {
        
          sub.clade.phylo = drop.tip(sub.phylo,tips.to.drop2);
          sub.clade.phylo$root.time = max(dist.nodes(sub.clade.phylo)[1:Ntip(sub.clade.phylo),Ntip(sub.clade.phylo) + 1]); sub.clade.phylo$root.time;
          sub.clade.phylo$origin.time = t - sub.clade.phylo$root.time; sub.clade.phylo$origin.time;
        
          if (identical(sort(as.integer(unique(sub.populations$spp.name))) ,
                        sort(as.integer(sub.clade.phylo$tip.label))) == F ) {
            print(c(c,t,'Error: trimmed phylogeny does not contain the correct species')); break} else{};
        
          reg.summary = regional.calc(sub.populations[,c('region','spp.name','time.of.origin','reg.env','extant')], sub.clade.phylo, as.integer(t));
        
          #Note that extinction calculation must be done on subset.populations, not sub.populations
          extinction = extinct.calc(subset.populations, timeslice=t)
          reg.summary2 = merge(reg.summary,extinction[,c('region','extinction.rate')],by='region')
      }
    }
  }
  }
  return(reg.summary2)
}

# testing regsum_fun function
tes <- regsum_fun(3465)

# for Niche conservatism tropical origin (c(3465:3474)) and temperate origin (c(3565:3574))
NCtrop <- c(3465:3474)
regsum_list_NCtrop <- vector("list", length = length(NCtrop))
for (i in seq_along(NCtrop)) {
  tes <- regsum_fun(NCtrop[[i]])
  regsum_list_NCtrop[[i]] <- tes
}
regframe_NCtrop <- rbind(regsum_list_NCtrop[[1]], regsum_list_NCtrop[[2]], regsum_list_NCtrop[[3]], 
                         regsum_list_NCtrop[[4]], regsum_list_NCtrop[[5]], regsum_list_NCtrop[[6]],
                         regsum_list_NCtrop[[7]], regsum_list_NCtrop[[8]], regsum_list_NCtrop[[9]],
                         regsum_list_NCtrop[[10]])
origin <- rep("Tropical", 56)
regframe_NC <- cbind(regframe_NCtrop, origin)


NCtemp <- c(3565:3574)
regsum_list_NCtemp <- vector("list", length = length(NCtemp))
for (i in seq_along(NCtemp)) {
  tes <- regsum_fun(NCtemp[[i]])
  regsum_list_NCtemp[[i]] <- tes
}
regframe_NCtemp <- rbind(regsum_list_NCtemp[[1]], regsum_list_NCtemp[[2]], regsum_list_NCtemp[[3]], 
                         regsum_list_NCtemp[[4]], regsum_list_NCtemp[[5]], regsum_list_NCtemp[[6]],
                         regsum_list_NCtemp[[7]], regsum_list_NCtemp[[8]], regsum_list_NCtemp[[9]],
                         regsum_list_NCtemp[[10]])
origin <- rep("Temperate", 58)
regframe_NCtemp <- cbind(regframe_NCtemp, origin)

regframe_NC <- rbind(regframe_NC, regframe_NCtemp)
save(regframe_NC, file = './data/regframe_NC.Rdata')

# for ecological limits of a tropical origin (c(4065:4074)) and a temperate origin (c(4075:4084))
ELtrop <- c(4065:4074)
regsum_list_ELtrop <- vector("list", length = length(ELtrop))
for (i in seq_along(ELtrop)) {
  tes <- regsum_fun(ELtrop[[i]])
  regsum_list_ELtrop[[i]] <- tes
}
regframe_ELtrop <- rbind(regsum_list_ELtrop[[1]], regsum_list_ELtrop[[2]], regsum_list_ELtrop[[3]], 
                         regsum_list_ELtrop[[4]], regsum_list_ELtrop[[5]], regsum_list_ELtrop[[6]],
                         regsum_list_ELtrop[[7]], regsum_list_ELtrop[[8]], regsum_list_ELtrop[[9]],
                         regsum_list_ELtrop[[10]])
origin <- rep("Tropical", 100)
regframe_EL <- cbind(regframe_ELtrop, origin)


ELtemp <- c(4075:4084)
regsum_list_ELtemp <- vector("list", length = length(ELtemp))
for (i in seq_along(ELtemp)) {
  tes <- regsum_fun(ELtemp[[i]])
  regsum_list_ELtemp[[i]] <- tes
}
regframe_ELtemp <- rbind(regsum_list_ELtemp[[1]], regsum_list_ELtemp[[2]], regsum_list_ELtemp[[3]], 
                         regsum_list_ELtemp[[4]], regsum_list_ELtemp[[5]], regsum_list_ELtemp[[6]],
                         regsum_list_ELtemp[[7]], regsum_list_ELtemp[[8]], regsum_list_ELtemp[[9]],
                         regsum_list_ELtemp[[10]])
origin <- rep("Temperate", 100)
regframe_ELtemp <- cbind(regframe_ELtemp, origin)

regframe_EL <- rbind(regframe_EL, regframe_ELtemp)
save(regframe_EL, file = './data/regframe_El.Rdata')

# Creating linear plots for each metric for Niche conservatism and ecological limits colored by origin
# richness vs temperature
pdf('./figures/niche_conservatism_richness_vs_temperature.pdf', width = 7*1.5)
NC_rvt <- ggplot(regframe_NC, aes(x = reg.env, y = richness, color = origin)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = origin), alpha = 0.2) +
  labs(x = "Temperature (°C)", y = "Species Richness", color = "Origin") +
  scale_color_manual(values = c('#FF0000', '#0000FF')) +
  scale_fill_manual(values = c('#FF0000', '#0000FF')) +
  guides(fill = F) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), 
        legend.position = "top", 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0, 1.5, 1.5, 1.5), "lines"),
        aspect.ratio = 1)
NC_rvt
dev.off()

grob1 <- ggplotGrob(NC_rvt)

pdf('./figures/energy_limits_richness_vs_temperature.pdf', width = 7*1.5)
EL_rvt <- ggplot(regframe_EL, aes(x = reg.env, y = richness, color = origin)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = origin), alpha = 0.2) +
  labs(x = "Temperature (°C)", y = "Species Richness", color = "Origin") +
  scale_color_manual(values = c('#FF0000', '#0000FF')) +
  scale_fill_manual(values = c('#FF0000', '#0000FF')) +
  guides(fill = F) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), 
        legend.position = "top", 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0, 1.5, 1.5, 1.5), "lines"),
        aspect.ratio = 1)
EL_rvt
dev.off()

grob2 <- ggplotGrob(EL_rvt)

# richness vs MRD
pdf('./figures/niche_conservatism_richness_vs_mrd.pdf', width = 7*1.5)
NC_rvm <- ggplot(regframe_NC, aes(x = richness, y = MRD, color = origin)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = origin), alpha = 0.2) +
  labs(x = "Species Richness", y = "Mean Root Distance (MRD)", color = "Climatic Origin") +
  scale_color_manual(values = c('#FF0000', '#0000FF')) +
  scale_fill_manual(values = c('#FF0000', '#0000FF')) +
  guides(fill = F) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), 
        legend.position = "none",
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        aspect.ratio = 1)
NC_rvm
dev.off()

grob3 <- ggplotGrob(NC_rvm)

pdf('./figures/energy_limits_richness_vs_mrd.pdf', width = 7*1.5)
EL_rvm <- ggplot(regframe_EL, aes(x = richness, y = MRD, color = origin)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = origin), alpha = 0.2) +
  labs(x = "Species Richness", y = "Mean Root Distance (MRD)", color = "Climatic Origin") +
  scale_color_manual(values = c('#FF0000', '#0000FF')) +
  scale_fill_manual(values = c('#FF0000', '#0000FF')) +
  guides(fill = F) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), 
        legend.position = "none", 
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        aspect.ratio = 1)
EL_rvm
dev.off()

grob4 <- ggplotGrob(EL_rvm)

# MRD vs temperature
pdf('./figures/niche_conservatism_temperature_vs_mrd.pdf', width = 7*1.5)
NC_mvt <- ggplot(regframe_NC, aes(x = reg.env, y = MRD, color = origin)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = origin), alpha = 0.2) +
  labs(x = "Temperature (°C)", y = "Mean Root Distance (MRD)", color = "Climatic Origin") +
  scale_color_manual(values = c('#FF0000', '#0000FF')) +
  scale_fill_manual(values = c('#FF0000', '#0000FF')) +
  guides(fill = F) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), 
        legend.position = "none", 
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        aspect.ratio = 1)
NC_mvt
dev.off()

grob5 <- ggplotGrob(NC_mvt)

pdf('./figures/energy_limits_temperature_vs_mrd.pdf', width = 7*1.5)
EL_mvt <- ggplot(regframe_EL, aes(x = reg.env, y = MRD, color = origin)) +
  geom_point() +
  geom_smooth(method = lm, aes(fill = origin), alpha = 0.2) +
  labs(x = "Temperature (°C)", y = "Mean Root Distance (MRD)", color = "Climatic Origin") +
  scale_color_manual(values = c('#FF0000', '#0000FF')) +
  scale_fill_manual(values = c('#FF0000', '#0000FF')) +
  guides(fill = F) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), 
        legend.position = "none", 
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
        aspect.ratio = 1)
EL_mvt
dev.off()

grob6 <- ggplotGrob(EL_mvt)

# pulling phylogenetic trees for each model
phy_nc <- read.tree(paste("./data/raw_sim_output",'/sim', 3565, '_out/SENC_phylo_sim', 
                            3565,'.tre',sep=''))
te <- sample(phy_nc$tip.label, 50, replace = F)
nc_drop <- phy_nc$tip.label[!(phy_nc$tip.label %in% te)]
phy_nc_clean <- drop.tip(phy = phy_nc, tip = nc_drop)

phy_el <- read.tree(paste("./data/raw_sim_output",'/sim', 4065, '_out/SENC_phylo_sim', 
                               4065,'.tre',sep=''))
ce <- sample(phy_el$tip.label, 50, replace = F)
drop_el <- phy_el$tip.label[!(phy_el$tip.label %in% ce)]
phy_el_clean <- drop.tip(phy = phy_el, tip = drop_el)

# creating ggtree objects so they can be included in grid.arrange
pdf('./figures/phy_sim_nc.pdf')
phy_nc_g <- ggplot(phy_nc_clean) +
  geom_tree() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 0, 1.5), "lines"))
dev.off()

grob7 <- ggplotGrob(phy_nc_g)

pdf('./figures/phy_sim_el.pdf')
phy_el_g <- ggplot(phy_el_clean) +
  geom_tree() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(1.5, 1.5, 0, 1.5), "lines"))
dev.off()

grob8 <- ggplotGrob(phy_el_g)

# calculating beta for these trees
phy_nc_b <- maxlik.betasplit(phy_nc, confidence.interval = "profile")
save(phy_nc_b, file = './data/stats/phy_nc_b.Rdata')

phy_el_b <- maxlik.betasplit(phy_el, confidence.interval = "profile")
save(phy_el_b, file = './data/stats/phy_el_b.Rdata')

# arranging plots for final graphic
pdf('./figures/figure1.pdf', height = 15, width = 10)
ggarrange(EL_rvt, NC_rvt, EL_rvm, NC_rvm, EL_mvt, NC_mvt, phy_el_g, phy_nc_g, 
          ncol = 2,
          labels = c("a)", "e)",
                     "b)", "f)", 
                     "c)", "g)",
                     "d)", "h)"))
dev.off()

# checking statistics
NCtemp_env_v_richness_stats <- cor.test(regframe_NCtemp$reg.env, regframe_NCtemp$richness, 
                                        method = "pearson")

NCtrop_env_v_richness_stats <- cor.test(regframe_NCtrop$reg.env, regframe_NCtrop$richness, 
                                        method = "pearson")

ELtemp_env_v_richness_stats <- cor.test(regframe_ELtemp$reg.env, regframe_ELtemp$richness, 
                                        method = "pearson")

ELtrop_env_v_richness_stats <- cor.test(regframe_ELtrop$reg.env, regframe_ELtrop$richness, 
                                        method = "pearson")

NCtemp_mrd_v_richness_stats <- cor.test(regframe_NCtemp$richness, regframe_NCtemp$MRD, 
                                        method = "pearson")

NCtrop_mrd_v_richness_stats <- cor.test(regframe_NCtrop$richness, regframe_NCtrop$MRD, 
                                        method = "pearson")

ELtemp_mrd_v_richness_stats <- cor.test(regframe_ELtemp$richness, regframe_ELtemp$MRD, 
                                        method = "pearson")

ELtrop_mrd_v_richness_stats <- cor.test(regframe_ELtrop$richness, regframe_ELtrop$MRD, 
                                        method = "pearson")

NCtemp_mrd_v_env_stats <- cor.test(regframe_NCtemp$reg.env, regframe_NCtemp$MRD, 
                                        method = "pearson")

NCtrop_mrd_v_env_stats <- cor.test(regframe_NCtrop$reg.env, regframe_NCtrop$MRD, 
                                   method = "pearson")

ELtemp_mrd_v_env_stats <- cor.test(regframe_ELtemp$reg.env, regframe_ELtemp$MRD, 
                                   method = "pearson")

ELtrop_mrd_v_env_stats <- cor.test(regframe_ELtrop$reg.env, regframe_ELtrop$MRD, 
                                   method = "pearson")

# converting stats into data frames
NCtempr <- c(NCtemp_env_v_richness_stats$estimate, NCtemp_mrd_v_richness_stats$estimate,
             NCtemp_mrd_v_env_stats$estimate)
NCtropr <- c(NCtrop_env_v_richness_stats$estimate, NCtrop_mrd_v_richness_stats$estimate,
             NCtrop_mrd_v_env_stats$estimate)
ELtempr <- c(ELtemp_env_v_richness_stats$estimate, ELtemp_mrd_v_richness_stats$estimate,
             ELtemp_mrd_v_env_stats$estimate)
ELtropr <- c(ELtrop_env_v_richness_stats$estimate, ELtrop_mrd_v_richness_stats$estimate,
             ELtrop_mrd_v_env_stats$estimate)

cor_df_sims <- rbind(NCtropr, NCtempr, ELtropr, ELtempr)
colnames(cor_df_sims) <- c("temp_v_rich", "rich_v_mrd", "mrd_v_temp")

# pulling beta values for statistics
# two trees have been dropped from niche conservatism tropical origin, 3 trees
# from Ecological Limits Tropical Origin, and because
# they're so large that maxlik.betasplit fails to calculate beta for them
ncto <- c(3466, 3467, 3468, 3470, 3471, 3472, 3473, 3474)
ncte <- 3565:3574
elto <- c(4065, 4066, 4068, 4069, 4070, 4072, 4073)
elte <- c(4076, 4078, 4079, 4080, 4082, 4083, 4084)

# beta_maker is a function which will calculate the beta value for each tree
# created by each sim repetition. sims is a vector of sim ids. it returns a list
# of beta values
beta_maker <- function(sims) {
  beta_vals <- vector("list", length = length(sims))
  for (i in seq_along(sims)) {
    simtree <- read.tree(paste("./data/raw_sim_output",'/sim', sims[[i]], '_out/SENC_phylo_sim', 
                               sims[[i]],'.tre',sep=''))
    simtree <- multi2di.phylo(simtree)
    simtree <- na.omit(simtree)
    simb <- maxlik.betasplit(simtree, confidence.interval = 'profile')
    beta_vals[[i]] <- simb
  }
  return(beta_vals)
}

# beta for Niche Conservatism Tropical Origin
ncto_b <- beta_maker(ncto)
ncto_b_mean <- mean(ncto_b[[1]]$max_lik, ncto_b[[2]]$max_lik, 
                    ncto_b[[3]]$max_lik, ncto_b[[4]]$max_lik,
                    ncto_b[[5]]$max_lik, ncto_b[[6]]$max_lik,
                    ncto_b[[7]]$max_lik, ncto_b[[8]]$max_lik)

# beta for Niche Conservatism Temperate Origin
ncte_b <- beta_maker(ncte)
ncte_b_mean <- mean(ncte_b[[1]]$max_lik, ncte_b[[2]]$max_lik, 
                    ncte_b[[3]]$max_lik, ncte_b[[4]]$max_lik,
                    ncte_b[[5]]$max_lik, ncte_b[[6]]$max_lik,
                    ncte_b[[7]]$max_lik, ncte_b[[8]]$max_lik,
                    ncte_b[[9]]$max_lik, ncte_b[[10]]$max_lik)

# beta for Ecological Limits Tropical Origin
elto_b <- beta_maker(elto)
elto_b_mean <- mean(elto_b[[1]]$max_lik, elto_b[[2]]$max_lik, 
                    elto_b[[3]]$max_lik, elto_b[[4]]$max_lik,
                    elto_b[[5]]$max_lik, elto_b[[6]]$max_lik,
                    elto_b[[7]]$max_lik)

# beta for Ecological Limits Temperate Origin
elte_b <- beta_maker(elte)
elte_b_mean <- mean(elte_b[[1]]$max_lik, elte_b[[2]]$max_lik, 
                    elte_b[[3]]$max_lik, elte_b[[4]]$max_lik,
                    elte_b[[5]]$max_lik, elte_b[[6]]$max_lik,
                    elte_b[[7]]$max_lik)

# Adding beta to the data frame with the correlation coefficients for the other
# metrics
beta_stat <- c(ncto_b_mean, ncte_b_mean, elto_b_mean, elte_b_mean)
cor_df_sims <- cbind(cor_df_sims, beta_stat)
save(cor_df_sims, file = './data/stats/cor_df_sims.Rdata')

# creating data frame for p-values
NCtempp <- c(NCtemp_env_v_richness_stats$p.value, NCtemp_mrd_v_richness_stats$p.value,
             NCtemp_mrd_v_env_stats$p.value)
NCtropp <- c(NCtrop_env_v_richness_stats$p.value, NCtrop_mrd_v_richness_stats$p.value,
             NCtrop_mrd_v_env_stats$p.value)
ELtempp <- c(ELtemp_env_v_richness_stats$p.value, ELtemp_mrd_v_richness_stats$p.value,
             ELtemp_mrd_v_env_stats$p.value)
ELtropp <- c(ELtrop_env_v_richness_stats$p.value, ELtrop_mrd_v_richness_stats$p.value,
             ELtrop_mrd_v_env_stats$p.value)

p_df_sims <- rbind(NCtropp, NCtempp, ELtropp, ELtempp)
colnames(p_df_sims) <- c("temp_v_rich", "rich_v_mrd", "mrd_v_temp")
save(p_df_sims, file = './data/stats/p_df_sims.Rdata')

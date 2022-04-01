#!/usr/bin/env Rscript
# This is the script where the original simulation is run

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

#(4) start analyses based on value of 'sim' which draws parameter values from sim.matrix
if (partial.analysis == 0) {which.sims = 1:max(sim.matrix$sim.id)};
#if (partial.analysis == 1) {which.sims = c(sim.matrix$sim.id[sim.matrix$carry.cap == 'on' & sim.matrix$energy.gradient == 'on' & sim.matrix$sim.id > 3464])}; # which.sims = c(read.csv(paste(analysis_dir,"/sims.to.analyze.csv",sep=""))$x)
if (partial.analysis == 1) {which.sims = c(3465:3474, 3565:3574, 4065:4074, 4075:4084)};

foo = foreach(sim=which.sims,.packages = package.vector,.combine='rbind') %dopar% {

  rm(list=c('all.populations', 'time.richness', 'phylo.out', 'params.out', 'output', 'sim.results'))
  output = numeric();
  
  # (5) read in simulation results for specified simulation from the output zip file 
  # -- (the final version should delete the code in the first if statement and assume a zip file exists)
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
        
          MRD.range = max(reg.summary$MRD,na.rm = T) - min(reg.summary$MRD,na.rm = T)
          MRD.mean = mean(reg.summary$MRD,na.rm = T)
          MRD.var = var(reg.summary$MRD,na.rm = T)
          MRD.rich.slope = lm(reg.summary$MRD ~ reg.summary$richness)$coefficients[2]
          MRD.env.slope = lm(reg.summary$MRD ~ reg.summary$reg.env)$coefficients[2]
          PSV.range = max(reg.summary$PSV,na.rm = T) - min(reg.summary$PSV,na.rm = T)
          PSV.mean = mean(reg.summary$PSV,na.rm = T)
          PSV.var = var(reg.summary$PSV,na.rm = T)
          PSV.rich.slope = lm(reg.summary$PSV ~ reg.summary$richness)$coefficients[2]
          PSV.env.slope = lm(reg.summary$PSV ~ reg.summary$reg.env)$coefficients[2]
          n.div.regions = length(reg.summary$region[reg.summary$richness > 1])
        
          PSV.MRD.slope = lm(reg.summary$PSV ~ reg.summary$MRD)$coefficients[2]
          
          corr.results = cbind(xregion.analysis(reg.summary2),MRD.range,MRD.mean,MRD.var,MRD.rich.slope,MRD.env.slope,
                               PSV.range,PSV.mean,PSV.var,PSV.rich.slope,PSV.env.slope,PSV.MRD.slope,n.div.regions)
        
          #Pybus & Harvey (2000)'s gamma statistic
          Gamma.stat = gammaStat(sub.clade.phylo)
        
          #Calculate Blum & Francois (2006)'s Beta metric of tree imbalance using apTreeshape package
          # --seems to bonk on very large phylogenies, so only try calculating for fewer than 6000 species
          #if(length(sub.phylo$tip.label) < 6000) {
            tree.beta.out = maxlik.betasplit(sub.clade.phylo)
            tree.beta = tree.beta.out$max_lik
          #} else {
            #tree.beta = NA
          #}
        
          #Calculate Blomberg's K for two traits: environmental optimum, and mean region of occurrence
          #spp.traits = aggregate(sub.populations$region, by = list(sub.populations$spp.name, sub.populations$env.opt),
          #                       function(x) mean(x, na.rm=T))
          #names(spp.traits) = c('spp.name','env.opt','region')
          
          #spp.env = spp.traits$env.opt
          #names(spp.env) = spp.traits$spp.name
          #BK.env = phylosig(sub.clade.phylo, spp.env[sub.clade.phylo$tip.label], method="K")
          
          #spp.reg = spp.traits$region
          #names(spp.reg) = spp.traits$spp.name
          #BK.reg = phylosig(sub.clade.phylo, spp.reg[sub.clade.phylo$tip.label], method="K")
        
          output = rbind(output, cbind(sim=sim,clade.id = c, time = t, corr.results, gamma.stat = Gamma.stat,
                                       clade.richness = length(unique(sub.populations$spp.name)), 
                                       #BK.env = BK.env , BK.reg = BK.reg, 
                                       tree.beta = tree.beta))
          print(paste(sim,sub.clade.loop.end,c,t,date(),length(sub.clade.phylo$tip.label),sep="   "));
          flush.console()
        } # end third else
      } # end sub clade for loop
    } # end second else
  } # end timeslice loop

  #write all of this output to files
  if (num.of.time.slices == 1) {write.csv(output,paste(analysis_dir,"/reanalysis/",sim,"_time",t,".csv",sep=""),quote=F,row.names=F)};
  if (num.of.time.slices > 1 & root.only == 0) {write.csv(output,paste(analysis_dir,"/NEW_Stats_sim",sim,"_mult_times_all_clades.csv",sep=""),quote=F,row.names=F)};
  if (num.of.time.slices > 1 & root.only == 1) {write.csv(output,paste(analysis_dir,"/NEW_Stats_sim",sim,"_mult_times_root_only.csv",sep=""),quote=F,row.names=F)};
  if (which.time.slices != -999 & root.only == 1) {write.csv(output,paste(analysis_dir,"/NEW_Stats_sim",sim,"_specific_times_root_only.csv",sep=""),quote=F,row.names=F)};
  if (which.time.slices != -999 & root.only == 0) {write.csv(output,paste(analysis_dir,"/NEW_Stats_sim",sim,"_specific_times_all_clades.csv",sep=""),quote=F,row.names=F)};
  if (time.sequence[1] != -999 & root.only == 1) {write.csv(output,paste(analysis_dir,"/NEW_Stats_sim",sim,"_time_seq_root_only.csv",sep=""),quote=F,row.names=F)};
  if (time.sequence[1] != -999 & root.only == 0) {write.csv(output,paste(analysis_dir,"/NEW_Stats_sim",sim,"_time_seq_all_clades.csv",sep=""),quote=F,row.names=F)};
    
  #FIXME: store these warnings to a file, along with sim.id? Or is this being done in the shell?
  #print(c(warnings(),sim.start,sim.end,analysis.end));
  
  # Add overall summary info
  sim.matrix[sim.matrix$sim.id==sim,'n.regions'] = length(unique(all.populations$region))
  sim.matrix[sim.matrix$sim.id==sim,'extant.S'] = nrow(extant.ornot[extant.ornot$x>0,])
  sim.matrix[sim.matrix$sim.id==sim,'extinct.S'] = length(extinct.species)
  sim.matrix[sim.matrix$sim.id==sim,'skipped.clades'] = skipped.clades # number of clades skipped over for analysis, summed over timeslices
  sim.matrix[sim.matrix$sim.id==sim,'skipped.times'] = skipped.times # number of time slices skipped over for analysis
  #sim.matrix[sim.matrix$sim.id==sim,'BK.reg'] = BK.reg # blomberg's K based on region
  #sim.matrix[sim.matrix$sim.id==sim,'BK.env'] = BK.env # blomberg's K based on environment

  write.csv(sim.matrix[sim.matrix$sim.id==sim,],paste(analysis_dir,"/reanalysis/sim.matrix.output.",sim,"_time",t,".csv",sep=""),quote=F,row.names=F);
  sim.matrix[sim.matrix$sim.id==sim,]
} # end sim loop

write.csv(foo,paste(analysis_dir,'/reanalysis/sim.matrix.output_',Sys.Date(),'.csv',sep=''),row.names=F)

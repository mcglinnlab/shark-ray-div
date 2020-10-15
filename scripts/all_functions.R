library(ggplot2)

# This script contains every function used in this procedure and should be sourced for:
# rasterize_polygons.R
# Phylogenetic_metrics.R
# data_analysis.R
# ecoregions.R
# Final_graphic.R

# functions for rasterization:
# richness rasterize function
# poly_files = directory where all the shape files are ('./data/polygon')
# raster_files = directory to be created ('./data/raster/sp')
# base_raster = raster upon which data will be added
# wanted_crs = projection the raster will be in
# output value is a new res_stack, name it and save it accordingly
load('./data/raster/oceans_raster.Rdata')
richness_rasterize <- function(poly_files, raster_files, base_raster, wanted_crs) {
  sp_poly_files <- dir(poly_files)
  sp_raster_files <- sub("json", "grd", sp_poly_files)
  sp_raster_files <- sub(" ", "_", sp_raster_files)
  
  dir.create(raster_files)
  
  cl <- makeCluster(24) 
  registerDoParallel(cl) 
  #clusterExport(cl, list('sp_poly_files', 'sp_raster_files'))
  
  foreach(i = seq_along(sp_poly_files),
          .packages = c("raster", "rgdal")) %dopar% {
            # read in 
            temp_poly = readOGR(dsn = paste0(poly_files, "/", sp_poly_files[i]),
                                layer = sub(".json", "", sp_poly_files[i]))
            # convert projection to cea
            temp_poly = spTransform(temp_poly, CRS(wanted_crs))
            # add field that will provide values when rasterized
            temp_poly@data$occur = 1
            # rasterize
            sp_raster = rasterize(temp_poly, base_raster, field = 'occur')
            writeRaster(sp_raster, 
                        filename = paste0(raster_files, "/", sp_raster_files[i]),
                        datatype = "LOG1S", overwrite = TRUE)
          }
  stopCluster(cl)
  
  sp_stack = stack(sapply(sp_raster_files, function(x) 
    paste0(raster_files, "/", x)))
  names(sp_stack) = sub('.grd', ' ', names(sp_stack))
  
  # create a list of stack at each resultion
  factor_val <- c(2, 4, 8, 16, 32)
  res_stack <- lapply(factor_val, function(x) 
    aggregate(sp_stack, fac = x, fun = sum) > 0)
  res_stack <- c(sp_stack, res_stack)
  return(res_stack)
}

# richness_plot function  
# creating a species richness layer for each resolution
# res_stack = res_stack created in above function (sp_res_stack)
# figure_name = file path for pdf ('./figures/species_richness_maps.pdf')
# mask = continents shapefile
# output is the richness list, name and save accordingly
richness_plot <- function(res_stack, figure_name, mask) {
  richness_list = lapply(res_stack, function(x)
    calc(x, fun = sum, na.rm = T))
  
  pdf(figure_name)
  for (i in 1:6) {
    test <- rasterize(mask, richness_list[[i]], getCover = T)
    richness_list[[i]][values(test) > 0.9] <- NA
    plot(richness_list[[i]],
         main=paste('resolution =', res(richness_list[[i]])))
    plot(mask, add = T, col = "black")
  }
  dev.off()
  return(richness_list)
}

# enviro_plot function to aggregate and plot environmental variables
# enviro_raster = environmental raster from initial_cleanup
# figure_name = file path for pdf
#mask = continents mask
# output is raster list, name and save accordingly
enviro_plot <- function(enviro_raster, figure_name, mask) {
  factor_val <- c(2, 4, 8, 16, 32)
  enviro_list <- lapply(factor_val, function (x)
    aggregate(enviro_raster, fac = x, fun = mean))
  enviro_list <- c(enviro_raster, enviro_list)
  return(enviro_list)
  pdf(figure_name)
  for (i in 1:6) {
    test <- rasterize(mask, enviro_list[[i]], getCover = T)
    enviro_list[[i]][values(test) > 0.9] <- NA
    plot(enviro_list[[i]], main = paste('resolution =', res(enviro_list[[i]])))
    plot(continents, add = T, col = "black")
  }
  dev.off()
}

# functions for phylogenetic metrics:
# mini_tree function to make subtrees based on clades
# poly_dir = file directory to species polygons
mini_tree <- function(poly_dir) {
  names_poly <- dir(poly_dir)
  names_poly <- sub(pattern = '.json', "", names_poly)
  keep <- shark_tree_clean$tip.label %in% names_poly
  drop <- shark_tree_clean$tip.label[!keep]
  new_tree <- drop.tip(shark_tree_clean, drop)
  return(new_tree)
}

# com_mat function makes community matrix for phylogenetic diversity
# poly_dir = file directory of species
# needed_tree = clade tree
# res_stack = res_stack of clade
# output value is community matrix for specified clade, name accordingly
com_mat <- function(poly_dir, needed_tree, res_stack) {
  names_poly <- dir(poly_dir)
  names_poly <- sub(pattern = '.json', "", names_poly)
  test_species <- names_poly %in% needed_tree$tip.label
  i <- 1
  mat_list <- vector("list", length = length(res_stack))
  for (j in 1:6) {
    mat_list[[j]] = matrix(NA, ncol = sum(test_species), 
                           nrow=length(res_stack[[j]][[i]]))
    colnames(mat_list[[j]]) = names_poly[test_species]
    icol = 1
    for (i in which(test_species)) {
      mat_list[[j]][ , icol] = raster::extract(res_stack[[j]][[i]], 
                                               1:nrow(mat_list[[j]]))
      icol = icol + 1
    }
    mat_list[[j]] = ifelse(is.na(mat_list[[j]]), 0, mat_list[[j]])
  }
  return(mat_list)
}

# mrd_runplot function for mean root distance
# needed_tree = clade tree
# needed_mat = community matrix for clade
# mask = continents
# figure_name = file name for pdf
# raster_base = empty raster to be filled
# output value is the mrd raster list, name and save accordingly
load('./data/raster/species_richness.Rdata')
load('./data/continent/continent_new.Rdata')
mrd_runplot <- function(needed_tree, needed_mat, mask, figure_name, raster_base) {
  phylo_bl1 <- compute.brlen(needed_tree, 1)
  all_dist <- dist.nodes(phylo_bl1)
  root_dist <- all_dist[length(needed_tree$tip.label) + 1, 
                        1:length(needed_tree$tip.label)]
  tips_to_root <- data.frame(spp.name=needed_tree$tip.label, root_dist)
  
  mrd_test_list <- vector("list", length = length(needed_mat)) 
  for (i in 1:length(needed_mat)) {
    allspecies = colnames(needed_mat[[i]])
    mrd_test_list[[i]] = rep(0, nrow(needed_mat[[i]]))
    for(j in 1:nrow(needed_mat[[i]])) {
      sp_list = data.frame(spp.name = allspecies[needed_mat[[i]][j, ] == 1])
      if (nrow(sp_list) > 0) {
        root_dist_tot <- merge(sp_list, tips_to_root, sort = F)
        mrd_test_list[[i]][j] <- mean(root_dist_tot$root_dist)
      }
    }
  }
  
  raster_list <- vector("list", length = 6)
  pdf(figure_name)
  for (i in 1:6) {
    mrd_raster <- raster_base[[i]]
    mrd_raster@data@values <- mrd_test_list[[i]]
    test <- rasterize(mask, mrd_raster, getCover = T)
    mrd_raster[values(test) > 0.9] <- NA
    plot(mrd_raster)
    plot(mask, add = T, col = "black")
    raster_list[[i]] <- mrd_raster
  }
  dev.off()
  return(raster_list)
}

# psv_runplot function for phylogenetic species diversity 
# needed_tree = clade tree
# needed_mat = community matrix for clade
# mask = continents mask
# figure_name = pdf file
# raster_base = empty raster to be filled
# output is raster list, name and save accordingly
psv_runplot <- function(needed_tree, needed_mat, mask, figure_name, raster_base) {
  psv_test_list <- vector("list", length = 6)
  for (i in 1:6) {
    psv_test <- psv(needed_mat[[i]], needed_tree)
    psv_test_list[[i]] <- psv_test
  }
  
  pdf(figure_name)
  raster_list <- vector("list", length = 6)
  for (i in 1:6) {
    psv_raster <- raster_base[[i]]
    psv_raster@data@values <- psv_test_list[[i]]$PSVs
    test <- rasterize(mask, psv_raster, getCover = T)
    psv_raster[values(test) > 0.9] <- NA
    plot(psv_raster)
    plot(mask, add = T, col = "black")
    raster_list[[i]] <- psv_raster
  }
  dev.off()
  return(raster_list)
}

# beta_runplot
# can only be applied to total and carcharhiniformes, lamniformes are too unstable
# needed_tree = clade tree
# needed_mat = community matrix
# mask = continents mask
# figure_name = pdf file
# raster_base = empty raster to be filled
# output is beta raster list, name and save accordingly
beta_runplot <- function(needed_tree, needed_mat, mask, figure_name, raster_base) {
  beta_list <- vector("list", length = length(needed_mat))
  for (i in 1:length(needed_mat)) {
    allspecies = colnames(needed_mat[[i]])
    beta_list[[i]] = rep(NA, nrow(needed_mat[[i]]))
    for(j in 1:nrow(needed_mat[[i]])) {
      sp_list = data.frame(spp.name = allspecies[needed_mat[[i]][j, ] == 1])
      if (nrow(sp_list) > 1) {
        drop_these <- needed_tree$tip.label[!(needed_tree$tip.label %in% sp_list$spp.name)]
        temp_tree <- drop.tip(needed_tree, drop_these)
        temp_tree <- multi2di(temp_tree)
        b <- maxlik.betasplit(temp_tree)
        beta_list[[i]][j] <- b$max_lik
      }
    }
  }
  
  raster_list <- vector("list", length = 6)
  pdf(figure_name)
  for (i in 1:6) {
    beta_raster <- raster_base[[i]]
    beta_raster@data@values <- beta_list[[i]]
    test <- rasterize(mask, beta_raster, getCover = T)
    beta_raster[values(test) > 0.9] <- NA
    plot(beta_raster)
    plot(mask, add = T, col = "black")
    raster_list[[i]] <- beta_raster
  }
  dev.off()
  return(raster_list)
}


#' functions for data analysis:
#' make_plot function plots two variables against one another
#' @param dat = data frame from which variables are pulled
#' @param x_var = x axis variable, in character form
#' @param y_var = y axis variable, in character form
#' @param x_lab = x axis label as character string
#' @param y_lab = y axis label as character string
#' @param figure_name = pdf file name
#' output is a list of vectors for each scale. The first element is the model, 
#' the second is the slope, the third is the standardized model, 
#' the fourth is the correlation coefficient, name and save accordingly
make_plot <- function(dat, x_var, y_var, x_lab, y_lab, figure_path) {
  stats_list <- vector("list", length = 6)
  pdf(figure_path, onefile = TRUE)
  for (i in 1:6) {
    tmp <- subset(dat, scale == i)
    x = eval(parse(text = paste0('tmp$', x_var)))
    y = eval(parse(text = paste0('tmp$', y_var)))
    p <- ggplot(tmp, mapping = aes(x,y)) +
      geom_point() +
      geom_smooth(method = 'lm') +
      xlab(x_lab) +
      ylab(y_lab)
    print(p)
    lin <- lm(data = tmp, y ~ x) 
    stats <- summary(lin)
    slope <- coef(stats)[2,1]
    x_std <- scale(x)
    y_std <- scale(y)
    mod_std <- lm(y_std ~ x_std)
    r_value <- coef(mod_std)[2]
    nest_list <- vector("list", length = 4)
    nest_list[[1]] <- lin
    nest_list[[2]] <- slope
    nest_list[[3]] <- mod_std
    nest_list[[4]] <- r_value
    stats_list[[i]] <- nest_list
  }
  dev.off()
  return(stats_list)
}

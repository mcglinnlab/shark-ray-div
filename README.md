# shark-ray-div
Global patterns of shark and ray diversity

All raw polygon data for this project can be downloaded from zenodo here: https://zenodo.org/record/6321610#.YkEfVefMK3A

Naylor, Gavin. (2022). Global Dataset of Shark and Chimaera Ranges [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6321610

Chlorophyll data was retrived from NASA Earth Observatory: https://neo.sci.gsfc.nasa.gov/view.php?datasetId=MY1DMM_CHLORA

All other environmental variables were retrieved from NOAA's World Ocean Database: https://www.nodc.noaa.gov/OC5/WOD/pr_wod.html

## Instructions for recreating analysis

1. Simulation
    * Code: `./scripts/Analysis_workflow_local.R`
        - Produces needed simulation results and adapted from Hurlbert and Stegen (2014). We specifically modified this to remove
        an if else which prevented the beta value for large trees from being computed. The sim IDs listed in `which.sims` corresponds
        to the following models:
            + 3465:3474 = Niche Conservatism Tropical Origin
            + 3565:3574 = Niche Conservatism Temperate Origin
            + 4065:4074 = Ecological Limits Tropical Origin
            + 4075:4084 = Ecological Limits Temperate Origin
    * Code: `./scripts/Graphics.R`
        - Creates the graphics for figure 1. It also creates two data frames, `data/stats/cor_df_sims.Rdata` and `data/stats/p_df_sims.Rdata`,
        which are used for calculating MSE.
2. Shark analysis
    * Code: `./scripts/Initial_cleanup.R`
        - Where all environmental rasters are read in, created, and cleaned
    * Code: `./scripts/Rasterize_polygons.R`
        - Where range map polygons are read in, transformed into rasters, stacked, and summed to create richness maps and rasters of species richness
        for the global analysis, the subclade analysis, and an IUCN test case
        - Environmental rasters are also plotted here
    * Code: `./scripts/Phylogenetic_metrics.R`
        - The whole tree is read in and cleaned, so duplicates and polytomies are removed
        - The `mini_tree` function creates subset trees for sub-analyses
        - The `com_mat` function creates community matrices to be used to calculate MRD
        - Rasters and maps are created for MRD
        - The beta statistic is calculated using the `maxlik.betasplit` function for the total analysis and the two subclade analyses
        - A ggtree object is created for the total shark tree to be used in the figure 2 graphic
   * Code: `./scripts/Data_analysis.R`
       - All raster data is pulled into a single data frame, and the `make_plot` function creates regression plots for all variables
       - Figure 2 is generated here using the `make_plot_fig2` function and `grid.arrange`
   * Code: `./scripts/Ecoregions.R`
        - Realm polygons are pulled in and cleaned
        - All components of the analysis (rasterization, phylogenetic metrics, and data analysis/graphing) are repeated in this script for the Tropical
        Atlantic and Central Indo-Pacific realms
   * Code: `./scripts/Final_graphic.R`
       - All correlation coefficients and beta values for the global analysis, subclade analyses, and ecoregion analyses are pulled into a list of data
       frames corresponding to scale resolution, along with the same statistics from the simulation
       - Data frames are similarly created for the upper and lower confidence interval of each statistic
   * Code: `Obs_pred_plot.R`
       - MSE is calculated for each analysis and hypothesis combination using the `table_list` of correlation coefficients and beta values and the
       `error_list` of confidence intervals
       - Figure 3 is created here along with other data visualizations describing how close each hypothesis is to fitting the real world data from each analysis

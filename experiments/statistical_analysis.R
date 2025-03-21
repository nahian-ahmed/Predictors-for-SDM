#################################################################################################
# Statistical analysis comparing environmental predictors across species. Run this after 
# generating all_metrics.RData with evaluation.R
#
# Laurel Hopkins 
# 10/2/2020
#################################################################################################

library(pROC)  
library(plotrix)
library(PRROC)
library(ggplot2)

Working.directory <-"C:\\Users\\Laurel\\Documents\\Oregon State\\Research\\ICB" 
setwd(Working.directory)


#########################################################
# 1a. Prepare data for statistical analysis - SUMMER MEANS
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean")
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
#                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
#                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
#                  "Yellow-breasted Chat")

data_df <- all_metrics[all_metrics$model %in% models & all_metrics$species %in% species_list,]
data_df <- data_df[!(data_df$feature_set=="texture"),]


#calculate percent difference in AUC from the mean of all AUCs per species
for(species in species_list) {
  raw_band_val <- data_df[data_df$species==species & data_df$feature_set=="raw_bands", "AUC"]
  # get the mean AUC for the species
  species_mean <- mean(data_df[data_df$species==species, "AUC"])
  data_df$species_mean[data_df$species==species] <- species_mean
  # calculate the difference of each index from the species mean
  data_df$diff_from_mean[data_df$species==species] <- data_df[data_df$species==species, "AUC"] - species_mean
  # calculate the percent difference of each index from the species mean
  data_df$percent_diff_from_mean[data_df$species==species] <- (data_df[data_df$species==species, "AUC"] - 
                                                                species_mean) / species_mean
  # calculate the difference od each index from the raw band value
  data_df$diff_from_raw_band[data_df$species==species] <- raw_band_val - data_df[data_df$species==species, "AUC"]
}


#########################################################
# 1b. Visualize data - SUMMER MEANS
#########################################################
library(ggplot2)

# box plots of percent difference from species means
qplot(feature_set, percent_diff_from_mean, data = data_df, geom = "boxplot") + 
   labs(x = "Spectral predictors", y = "Percent change from species mean AUC", 
        title = "Percent change in AUC from species mean AUC") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(limits = c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi", "raw_bands", 
                              "tasseled_cap"), 
                   labels = c("evi"="EVI", "msavi"="MSAVI", "savi"="SAVI", "nbr"="NBR", "nbr_2"="NBR2", 
                              "ndmi"="NDMI", "ndsi"="NDSI", "ndvi"="NDVI", "raw_bands"="Raw bands", 
                              "tasseled_cap"="Tasseled cap"))

# RESULTS: Equal variance assumption is violated; need a non-parametric test.


# box plots of mean AUCs
qplot(feature_set, AUC, data = data_df, geom = "boxplot") + 
  labs(x = "Spectral predictors", y = "Mean AUC across species", title = "AUC averaged across all species") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(limits = c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi", "raw_bands", 
                              "tasseled_cap"), 
                   labels = c("evi"="EVI", "msavi"="MSAVI", "savi"="SAVI", "nbr"="NBR", "nbr_2"="NBR2", 
                              "ndmi"="NDMI", "ndsi"="NDSI", "ndvi"="NDVI", "raw_bands"="Raw bands", 
                              "tasseled_cap"="Tasseled cap"))


#########################################################
# 1c. Summary statistics - SUMMER MEANS
#########################################################
# means
raw_bands_mean <- mean(data_df[data_df$feature_set=="raw_bands","AUC"])
tasseled_cap_mean <- mean(data_df[data_df$feature_set=="tasseled_cap","AUC"])
ndvi_mean <- mean(data_df[data_df$feature_set=="ndvi","AUC"])
single_indices <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi")
single_indices_mean <- mean(data_df[data_df$feature_set %in% single_indices,"AUC"])
single_indices_wo_ndvi <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi")
single_indices_wo_ndvi_mean <- mean(data_df[data_df$feature_set %in% single_indices_wo_ndvi,"AUC"])

# difference from raw bands mean
tasseled_cap_diff <- (tasseled_cap_mean - raw_bands_mean);tasseled_cap_diff
ndvi_diff <- (ndvi_mean - raw_bands_mean);ndvi_diff
single_indices_diff <- (single_indices_mean - raw_bands_mean);single_indices_diff
single_indices_wo_ndvi_diff <- (single_indices_wo_ndvi_mean - raw_bands_mean);single_indices_wo_ndvi_diff

# species mean AUC across summer means spectral predictors
for(species in species_list) {
  species_df <- data_df[data_df$species==species & data_df$model=="mean",]
  mean_AUC <- mean(species_df$AUC)
  print(paste0(species, ": ", mean_AUC))
}

# mean raw bands summer mean AUC across species
raw_bands_df <- data_df[data_df$feature_set=="raw_bands" & data_df$model=="mean",]
raw_bands_mean_AUC <- mean(raw_bands_df$AUC)
print(paste0("raw bands mean AUC: ", raw_bands_mean_AUC))

# highest summer means spectral predictor per species
for(species in species_list) {
  species_df <- data_df[data_df$species==species & data_df$model=="mean",]
  max_feature_set <- species_df[which.max(species_df$AUC),"feature_set"]
  print(paste0(species, ": ", max_feature_set))
  
}

# highest summer means single-index spectral predictor per species
for(species in species_list) {
  species_df <- data_df[data_df$species==species & data_df$feature_set %in% single_indices 
                        & data_df$model=="mean",]
  max_feature_set <- species_df[which.max(species_df$AUC),"feature_set"]
  print(paste0(species, ": ", max_feature_set))
  
}

# mean AUC for single-index spectral predictors across species
for(index in single_indices) {
  species_df <- data_df[data_df$feature_set==index & data_df$model=="mean",]
  mean_AUC <- mean(species_df$AUC)
  print(paste0(index, ": ", mean_AUC))
  
}
#########################################################
# 1d. Friedma's test for non-parametric one-way repeated 
#    measures ANOVA - SUMMER MEANS
#########################################################
# Need repeasted measures test since subjects (i.e. species) are repeatedly 
# measured across all groups (violates independence assumption)
#
# Need a non-parametric test since equal variance assumption is also 
# violated

# select only relevant data
fried_df <- data_df[, c("species", "feature_set", "percent_diff_from_mean"),]
# reshape for friedman
reshpaed_fried_df <- reshape(fried_df, idvar = "species", timevar = "feature_set", direction ="wide")
# convert to matrix and drop species column
reshpaed_fried_df <- as.matrix(reshpaed_fried_df[c(-1)])
# run test
friedman.test(reshpaed_fried_df)


#########################################################
# 2a. Prepare data for statistical analysis - ALL SEASON MEANS
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "mean_allSeasons")
species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
                  "Yellow-breasted Chat")
all_seasons_data_df <- all_metrics[all_metrics$model %in% models & all_metrics$species %in% species_list,]
# Remove NDSI
#all_seasons_data_df <- all_seasons_data_df[!(all_seasons_data_df$feature_set=="ndsi"),]

#calculate absolute difference in AUC from the summer mean and spring, summer, fall mean per species
spectral_predictors <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi", "raw_bands", 
                         "tasseled_cap")
for(species in species_list) {
  for(predictor in spectral_predictors) {
    summer_mean <- all_seasons_data_df[all_seasons_data_df$model == "mean" & 
                                         all_seasons_data_df$species == species & 
                                         all_seasons_data_df$feature_set == predictor,"AUC"]
    allSeasons_mean <- all_seasons_data_df[all_seasons_data_df$model == "mean_allSeasons" & 
                                         all_seasons_data_df$species == species & 
                                         all_seasons_data_df$feature_set == predictor,"AUC"]
    all_seasons_data_df[all_seasons_data_df$model == "mean_allSeasons" & 
                                           all_seasons_data_df$species == species & 
                                           all_seasons_data_df$feature_set == predictor, "diff_from_summer"] <- 
      allSeasons_mean - summer_mean
  }
}


#########################################################
# 2b. Visualize data - ALL SEASON MEANS
#########################################################
library(ggplot2)

# box plots of AUCs excluding all indices expect for ndvi (highest performing index)
to_plot <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndvi", "raw_bands", "tasseled_cap")
qplot(feature_set, diff_from_summer, data = all_seasons_data_df[all_seasons_data_df$feature_set %in% to_plot & 
                                                                  all_seasons_data_df$model == "mean_allSeasons",], 
      geom = "boxplot") + 
  labs(x = "Spectral predictors", y = "Increase in AUC", title = "Increase in AUC by adding spring and fall means to summer means") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(limits = to_plot, 
                   labels = c("evi"="EVI", "msavi"="MSAVI", "savi"="SAVI", "nbr"="NBR", "nbr_2"="NBR2", 
                              "ndmi"="NDMI", "ndsi"="NDSI", "ndvi"="NDVI", "raw_bands"="Raw bands", 
                              "tasseled_cap"="Tasseled cap", "texture"="Texture"))
  
  
#########################################################
# 2c. Summary statistics - ALL SEASON MEANS
#########################################################
# mean incresae in AUC by adding spring and fall (all predictors)
mean(all_seasons_data_df[,"diff_from_summer"], na.rm=TRUE)

# mean incresae in AUC by adding spring and fall (raw bands and Tasseled Cap)
top_predictors <- c("raw_bands", "tasseled_cap")
top_mean_diff <- mean(all_seasons_data_df[all_seasons_data_df$feature_set %in% top_predictors, 
                                          "diff_from_summer"], na.rm=TRUE);top_mean_diff

# mean incresae in AUC by adding spring and fall (single-index models)
index_predictors <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi")
index_mean_diff <- mean(all_seasons_data_df[all_seasons_data_df$feature_set %in% index_predictors, 
                                            "diff_from_summer"], na.rm=TRUE);index_mean_diff


#########################################################
# 3a. Prepare data for statistical analysis - SUMMER MEAN & 
#     STD. DEV
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "mean_stdDev")
species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
                  "Yellow-breasted Chat")
mean_stdDev_data_df <- all_metrics[all_metrics$model %in% models & all_metrics$species %in% species_list,]
# Remove NDSI
# mean_stdDev_data_df <- mean_stdDev_data_df[!(mean_stdDev_data_df$feature_set=="ndsi"),]

# calculate absolute difference in AUC from the summer mean and spring, summer, fall mean per species
spectral_predictors <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi", "raw_bands", 
                         "tasseled_cap")
for(species in species_list) {
  for(predictor in spectral_predictors) {
    summer_mean <- mean_stdDev_data_df[mean_stdDev_data_df$model == "mean" & 
                                         mean_stdDev_data_df$species == species & 
                                         mean_stdDev_data_df$feature_set == predictor,"AUC"]
    summer_mean_stdDev <- mean_stdDev_data_df[mean_stdDev_data_df$model == "mean_stdDev" & 
                                                mean_stdDev_data_df$species == species & 
                                                mean_stdDev_data_df$feature_set == predictor,"AUC"]
    mean_stdDev_data_df[mean_stdDev_data_df$model == "mean_stdDev" & 
                          mean_stdDev_data_df$species == species & 
                          mean_stdDev_data_df$feature_set == predictor, "diff_from_mean"] <- 
      summer_mean_stdDev - summer_mean
  }
}


#########################################################
# 3b. Visualize data - SUMMER MEAN & STD. DEV
#########################################################
library(ggplot2)

# box plots of AUCs excluding all indices expect for ndvi (highest performing index)
to_plot <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndvi", "raw_bands", "tasseled_cap")
qplot(feature_set, diff_from_mean, data = mean_stdDev_data_df[mean_stdDev_data_df$feature_set %in% to_plot & 
                                                                  mean_stdDev_data_df$model == "mean_stdDev",], 
      geom = "boxplot") + 
  labs(x = "Spectral predictors", y = "Increase in AUC", title = "Increase in AUC by adding standard deviation to summer means") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(limits = to_plot, labels = c("evi"="EVI", "msavi"="MSAVI", "savi"="SAVI", "nbr"="NBR", 
                                                "nbr_2"="NBR2", "ndmi"="NDMI", "ndvi"="NDVI", 
                                                "raw_bands"="Raw bands", "tasseled_cap"="Tasseled cap"))


#########################################################
# 3c. Summary statistics - SUMMER MEAN & STD. DEV
#########################################################
# mean incresae in AUC by adding standard deviation (all predictors)
mean(mean_stdDev_data_df[,"diff_from_mean"], na.rm=TRUE)

# mean incresae in AUC by adding standard deviation (raw bands and Tasseled Cap)
top_predictors <- c("raw_bands", "tasseled_cap")
top_mean_diff <- mean(mean_stdDev_data_df[mean_stdDev_data_df$feature_set %in% top_predictors, 
                                          "diff_from_mean"], na.rm=TRUE);top_mean_diff

# mean incresae in AUC by adding spring and fall (single-index models)
index_predictors <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi")
index_mean_diff <- mean(mean_stdDev_data_df[mean_stdDev_data_df$feature_set %in% index_predictors, 
                                            "diff_from_mean"], na.rm=TRUE);index_mean_diff

# calculate increase in performance for riparian species 
riparian_species <- c("Yellow Warbler", "Yellow-breasted Chat")
riparian_diff <- mean(mean_stdDev_data_df[mean_stdDev_data_df$model == "mean_stdDev" & 
                                            mean_stdDev_data_df$species %in% riparian_species & 
                                            mean_stdDev_data_df$feature_set %in% all_predictors, 
                                            "diff_from_mean"]);riparian_diff


#########################################################
# 4a. Prepare data for statistical analysis - SUMMER MEAN & 
#     TEXTURE
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "mean_texture")
species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
                  "Yellow-breasted Chat")
mean_texture_data_df <- all_metrics[all_metrics$model %in% models & all_metrics$species %in% species_list,]
# Remove NDSI
#mean_texture_data_df <- mean_texture_data_df[!(mean_texture_data_df$feature_set=="ndsi"),]

#calculate absolute difference in AUC from the summer mean and spring, summer, fall mean per species
spectral_predictors <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi", "raw_bands", 
                         "tasseled_cap")
for(species in species_list) {
  for(predictor in spectral_predictors) {
    summer_mean <- mean_texture_data_df[mean_texture_data_df$model == "mean" & 
                                          mean_texture_data_df$species == species & 
                                          mean_texture_data_df$feature_set == predictor,"AUC"]
    summer_mean_stdDev <- mean_texture_data_df[mean_texture_data_df$model == "mean_texture" & 
                                                 mean_texture_data_df$species == species & 
                                                 mean_texture_data_df$feature_set == predictor,"AUC"]
    mean_texture_data_df[mean_texture_data_df$model == "mean_texture" & 
                           mean_texture_data_df$species == species & 
                           mean_texture_data_df$feature_set == predictor, "diff_from_mean"] <- 
      summer_mean_stdDev - summer_mean
  }
}


#########################################################
# 5a. Prepare data for statistical analysis - RAW BANDS 
#     Comparison
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "mean_allSeasons", "mean_stdDev", "mean_texture", "mean_stdDev_allSeasons",
            "mean_allSeasons_texture", "mean_stdDev_texture", "mean_stdDev_allSeasons_texture")
spectral_predictor <- "raw_bands"
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
#                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
#                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
#                  "Yellow-breasted Chat")

data_df <- all_metrics[all_metrics$model %in% models & all_metrics$feature_set==spectral_predictor &
                         all_metrics$species %in% species_list,]


#calculate percent difference in AUC from the raw bands AUC by species
for (species in species_list) {
  raw_band_val <- data_df[data_df$species==species & data_df$model=="mean" & data_df$feature_set=="raw_bands", "AUC"]
  for (model in models[-1]) {
    # calculate the difference of each index from the species raw bands summer mean 
    data_df$diff_from_raw_bands[data_df$species==species & data_df$model==model] <- 
      data_df[data_df$species==species & data_df$model==model, "AUC"] - raw_band_val
    # calculate the percent difference of each index from the species mean
    data_df$percent_diff_from_raw_bands[data_df$species==species & data_df$model==model] <- 
      (data_df[data_df$species==species & data_df$model==model, "AUC"] - raw_band_val) / raw_band_val
  }
}


#########################################################
# 5b. Visualize data - RAW BANDS COMPARISON
#########################################################
library(ggplot2)

# keep everything but the mean model
data_df <- data_df[!(data_df$model =="mean"),]

# box plots of percent difference from species means
qplot(model, diff_from_raw_bands, data = data_df, geom = "boxplot") + 
  labs(x = "Summary methods", y = "Increase in AUC from raw bands summer means", 
       title = "Increase in AUC from raw bands summer means averaged over species") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(guide=guide_axis(angle = 45),
                   limits = c("mean_allSeasons", "mean_stdDev", "mean_texture", "mean_stdDev_allSeasons", 
                              "mean_allSeasons_texture", "mean_stdDev_texture", "mean_stdDev_allSeasons_texture"),
                   labels = c("mean_allSeasons"="Seasons", "mean_stdDev"="Std. Dev.", "mean_texture"="Texture", 
                              "mean_stdDev_allSeasons"="Seasons + Std. Dev.", 
                              "mean_allSeasons_texture"="Seasons + Texture",
                              "mean_stdDev_texture"="Std. Dev. + Texture", 
                              "mean_stdDev_allSeasons_texture"="Seasons + Std. Dev. + Texture"))


#########################################################
# 6a. Prepare data for statistical analysis - NDVI 
#     Comparison
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "mean_allSeasons", "mean_stdDev", "mean_texture", "mean_stdDev_allSeasons",
            "mean_allSeasons_texture", "mean_stdDev_texture", "mean_stdDev_allSeasons_texture")
spectral_predictor <- "ndvi"
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
#                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
#                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
#                  "Yellow-breasted Chat")

data_df <- all_metrics[all_metrics$model %in% models & all_metrics$feature_set==spectral_predictor &
                         all_metrics$species %in% species_list,]


#calculate percent difference in AUC from the raw bands AUC by species
for (species in species_list) {
  ndvi_val <- data_df[data_df$species==species & data_df$model=="mean" & data_df$feature_set=="ndvi", "AUC"]
  for (model in models[-1]) {
    # calculate the difference of each index from the species raw bands summer mean 
    data_df$diff_from_raw_bands[data_df$species==species & data_df$model==model] <- 
      data_df[data_df$species==species & data_df$model==model, "AUC"] - ndvi_val
    # calculate the percent difference of each index from the species mean
    data_df$percent_diff_from_raw_bands[data_df$species==species & data_df$model==model] <- 
      (data_df[data_df$species==species & data_df$model==model, "AUC"] - ndvi_val) / ndvi_val
  }
}


#########################################################
# 6b. Visualize data - NDVI COMPARISON
#########################################################
library(ggplot2)

# keep everything but the mean model
data_df <- data_df[!(data_df$model =="mean"),]

# box plots of percent difference from species means
qplot(model, diff_from_raw_bands, data = data_df, geom = "boxplot") + 
  labs(x = "Summary methods", y = "Increase in AUC from NDVI summer mean", 
       title = "Increase in AUC from NDVI summer mean averaged over species") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(guide=guide_axis(angle = 45),
                   limits = c("mean_allSeasons", "mean_stdDev", "mean_texture", "mean_stdDev_allSeasons", 
                              "mean_allSeasons_texture", "mean_stdDev_texture", "mean_stdDev_allSeasons_texture"),
                   labels = c("mean_allSeasons"="Seasons", "mean_stdDev"="Std. Dev.", "mean_texture"="Texture", 
                              "mean_stdDev_allSeasons"="Seasons + Std. Dev.", 
                              "mean_allSeasons_texture"="Seasons + Texture",
                              "mean_stdDev_texture"="Std. Dev. + Texture", 
                              "mean_stdDev_allSeasons_texture"="Seasons + Std. Dev. + Texture"))


#########################################################
# 7a. Prepare data for statistical analysis - CLASSIFIED
#########################################################
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "modis", "nlcd", "nlcd_binary")
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
#                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
#                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
#                  "Yellow-breasted Chat")
classified_data_df <- all_metrics[all_metrics$model %in% models & all_metrics$species %in% species_list,]


#calculate absolute difference in AUC from the raw bands summer mean and the classified land cover predictors per species
classified_data_df$diff_from_raw_bands <- NA
for(species in species_list) {
  raw_band_mean <- classified_data_df[classified_data_df$model == "mean" & 
                                        classified_data_df$species == species & 
                                        classified_data_df$feature_set == "raw_bands","AUC"]
  for(model in models[-1]) {
    model_mean <- classified_data_df[classified_data_df$model == model & 
                                                 classified_data_df$species == species,"AUC"]
    classified_data_df[classified_data_df$model == model & classified_data_df$species == species, 
                       "diff_from_raw_bands"] <- model_mean - raw_band_mean
  }
}


#########################################################
# 7b. Visualize data - CLASSIFIED
#########################################################
library(ggplot2)

# box plots of AUCs excluding all indices expect for ndvi (highest performing index)
qplot(model, diff_from_raw_bands, data = classified_data_df[classified_data_df$model %in% models[-1],], 
      geom = "boxplot") + 
  labs(x = "Spectral predictors", y = "Difference in AUC", title = "Difference in AUC from raw bands summer means") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(labels = c("modis"="MODIS", "nlcd"="NLCD", "nlcd_binary"="discretized NLCD"))

# plot mean AUC for each predictor set
to_keep <- c("raw_bands", "modis", "nlcd", "nlcd_binary_all_classes")
to_plot <- classified_data_df[classified_data_df$feature_set %in% to_keep,]

qplot(feature_set, AUC, data = to_plot, geom = "boxplot") + 
  labs(x = "Spectral predictors", y = "Mean AUC", title = "Mean AUC across species") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18)) + 
  scale_x_discrete(labels = c("raw_bands"="Raw bands", "modis"="MODIS", "nlcd"="NLCD", "nlcd_binary"="discretized NLCD"))


#########################################################
# 7c. Summary statistics - CLASSIFIED
#########################################################
raw_bands_mean <- mean(classified_data_df[classified_data_df$model == "mean" & 
                                           classified_data_df$feature_set == "raw_bands" & 
                                          classified_data_df$species %in% species_list, 
                                           "AUC"]);raw_bands_mean

modis_diff <- mean(classified_data_df[classified_data_df$model == "modis" & 
                                            classified_data_df$species %in% species_list, 
                                      "diff_from_raw_bands"]);modis_diff

nlcd_diff <- mean(classified_data_df[classified_data_df$model == "nlcd" & 
                                        classified_data_df$species %in% species_list, 
                                      "diff_from_raw_bands"]);nlcd_diff

nlcd_binary_mean <- mean(classified_data_df[classified_data_df$model == "nlcd_binary" & 
                                        classified_data_df$species %in% species_list, 
                                      "diff_from_raw_bands"]);nlcd_binary_mean


#########################################################
# 7d. Friedma's test for non-parametric one-way repeated 
#    measures ANOVA - CLASSIFIED
#########################################################
# Need repeasted measures test since subjects (i.e. species) are repeatedly 
# measured across all groups (violates independence assumption)
#
# Need a non-parametric test since equal variance assumption is also 
# violated

# calculate percent difference from species mean to control for species
classified_data_df$diff_from_mean <- NA
classified_data_df$species_mean <- NA
for(species in species_list) {
  # get the mean AUC for the species
  species_mean <- mean(classified_data_df[classified_data_df$species==species, "AUC"])
  classified_data_df$species_mean[classified_data_df$species==species] <- species_mean
  # calculate the difference of each index from the species mean
  classified_data_df$diff_from_mean[classified_data_df$species==species] <- 
    classified_data_df[classified_data_df$species==species, "AUC"] - species_mean
  # calculate the percent difference of each index from the species mean
  classified_data_df$percent_diff_from_mean[classified_data_df$species==species] <- 
    (classified_data_df[classified_data_df$species==species, "AUC"] - species_mean) / species_mean
}

# select only relevant data
features_to_keep <- c("raw_bands", "modis", "nlcd", "nlcd_binary_all_classes")
classified_fried_df <- classified_data_df[classified_data_df$feature_set %in% features_to_keep,
                                          c("species", "model", "percent_diff_from_mean"),]
# reshape for friedman
reshpaed_classified_fried_df <- reshape(classified_fried_df, idvar = "species", timevar = "model", 
                                        direction ="wide")
# convert to matrix and drop species column
reshpaed_classified_fried_df <- as.matrix(reshpaed_classified_fried_df[c(-1)])
# run test
friedman.test(reshpaed_classified_fried_df)


#########################################################
# 7e. Nemenyi post-hoc analysis to determine which groups
#    are significantly different - CLASSIFIED
#########################################################
library(PMCMRplus)

#posthoc.friedman.nemenyi.test(reshpaed_classified_fried_df)
frdAllPairsNemenyiTest(reshpaed_classified_fried_df)


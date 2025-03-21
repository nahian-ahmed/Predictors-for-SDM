#################################################################################################
# Remote Sensing Indices Paper
# Summary: Compare habitat suitability models based on different satellite-derived indices
#
# Tyler Hallman
# 05/12/2020
# UPDATED 05/22/2020: Added for loop for running multiple models and necessary changing 
# name for model outputs.
# UPDATED 06/09/2020 by Laurel Hopkins: Adding code to parallelize species and models 
# (environmental variables)
#################################################################################################

# If running from RStudio, set get.inputs to FALSE and set the model and species on line 71 & 72.
# Otherwise, set get.inputs to TRUE and run the script with:
# R run_model.R --species "Ash-throated Flycatcher" --model "nlcd_binary_all_classes"
get.inputs <- FALSE  

model.buffers <- c(75, 600, 2400)
model.seasons <- c("summer", "spring", "fall")
model.segmentation <- "nbr"
model.summaries <- "mean_stdDev"  #"mean"

###### Set paths ###### 
Working.directory <-"C:\\Users\\Laurel\\Documents\\Oregon State\\Research\\ICB\\" # This is if I'm running it on my computer


#########################################################
# 1. Load packages
#########################################################

library(lme4)
library(MuMIn)
library(broom)
library(ggpubr)
library(reshape2)
library(DataCombine)
library(randomForest)
library(precrec)
library(data.table)
library(shinydashboard)
library(optparse)
library(hash)
library(blockCV)


#########################################################
# 2. Get user input
#########################################################

###### File Names ###### 
# Environmental files
environmental_file <- "Landsat_buffered_summaries"
modis_file <- "MODIS_MCD12Q1_proportions.csv"
nlcd_file <- "NLCD2016_proportions.csv"
nlcd_binary_file <- "NLCD2016_binary_indicators.csv"

# Species counts
species_file <- "Ash-throated_Flycatcher_species_counts.csv"

# blockCV spatial split
blockCV_file <- "sb1_Ash-throated_Flycatcher.RData"


###### Model parameters ######
if (model.params){
  option_list <- list(
    make_option("--model", action="store", default=NA, type="character", 
                metavar="model", help="Model to run (i.e. RS indices)"),
    make_option("--species", action="store", default=NA, type="character",
                metavar="species", help="Species to analyze")
  )
  args <- parse_args(OptionParser(option_list = option_list))
  
  model <- args$model
  species <- args$species

} else { # define  manually 
  
  # Model options: "raw_bands", "tasseled_cap", "ndvi", "ndmi", "nbr", "nbr_2", "evi", "savi", "msavi", "ndsi"
  model <- "raw_bands" 
  species <- "Ash-throated Flycatcher"
}

print(paste0("Model: ", model))
print(paste0("Species: ", species)) 


#########################################################
# 3. Load data 
#########################################################

###### Set the working directory ######
setwd(Working.directory)

print(blockCV_file)

###### Load species data ######

bird_data <- read.csv(file=species_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)

###### Load environmental data ######
if (model == "modis") {
  # use modis data for environmental features
  environ_data <- read.csv(file=modis_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)
} else if (model == "nlcd") {
  environ_data <- read.csv(file=nlcd_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)
} else if (model == "nlcd_binary") {
  environ_data <- read.csv(file=nlcd_binary_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)
} else {
  environ_data <- read.csv(file=environmental_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)
}

###### Load BlockCV data ######
load(blockCV_file)


#########################################################
# 4. Select environ data and merge with species counts
#    and blockCV folds
#########################################################
###### Select explanatory variables ###### 
## Buffers (6): 75, 600, 2400 m
## Seasons (4): Spring, Fall, Early Summer, Late Summer


## List of models to run
## List of corresponding variables
if (model == "raw_bands") {
  if (model.summaries == "mean") {
    model.indices <- c("B1_mean", "B2_mean", "B3_mean", "B4_mean", "B5_mean", "B7_mean")
  } else if (model.summaries == "mean_stdDev") {
    model.indices <- c("B1_mean", "B2_mean", "B3_mean", "B4_mean", "B5_mean", "B7_mean",
                     "B1_stdDev", "B2_stdDev", "B3_stdDev", "B4_stdDev", "B5_stdDev", "B7_stdDev")}
} else if (model == "tasseled_cap") {
  if (model.summaries == "mean") {
    model.indices <- c("TCA_mean", "TCB_mean", "TCG_mean", "TCW_mean")
  } else if (model.summaries == "mean_stdDev") {
    model.indices <- c("TCA_mean", "TCB_mean", "TCG_mean", "TCW_mean",
                     "TCA_stdDev", "TCB_stdDev", "TCG_stdDev", "TCW_stdDev")}
} else if (model == "ndvi") {
  if (model.summaries == "mean") { 
    model.indices <- c("NDVI_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("NDVI_mean", "NDVI_stdDev")}
} else if (model == "ndmi") {
  if (model.summaries == "mean") { 
    model.indices <- c("NDMI_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("NDMI_mean", "NDMI_stdDev")}
} else if (model == "ndsi") {
  if (model.summaries == "mean") { 
    model.indices <- c("NDSI_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("NDSI_mean", "NDSI_stdDev")}
} else if (model == "nbr") {
  if (model.summaries == "mean") { 
    model.indices <- c("NBR_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("NBR_mean", "NBR_stdDev")}
} else if (model == "nbr_2") {
  if (model.summaries == "mean") { 
    model.indices <- c("NBR_2_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("NBR_2_mean", "NBR_2_stdDev")}
} else if (model == "evi") {
  if (model.summaries == "mean") { 
    model.indices <- c("EVI_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("EVI_mean", "EVI_stdDev")}
} else if (model == "savi") {
  if (model.summaries == "mean") { 
    model.indices <- c("SAVI_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("SAVI_mean", "SAVI_stdDev")}
} else if (model == "msavi") {
  if (model.summaries == "mean") { 
    model.indices <- c("MSAVI_mean")
  } else if (model.summaries == "mean_stdDev") { 
    model.indices <- c("MSAVI_mean", "MSAVI_stdDev")}
} else if (model == "modis") {
  model.indices <- grep("pland", colnames(environ_data), value=TRUE)
} else if (model == "nlcd" | model == "nlcd_binary") {
  model.indices <- c("Barren_Land", "Cultivated_Crops", "Deciduous_Forest", "Developed_High_Intensity",
                    "Developed_Low_Intensity", "Developed_Medium_Intensity", "Developed_Open_Space",
                    "Emergent_Herbaceous_Wetlands", "Evergreen_Forest", "Grassland_Herbaceous",	
                    "Mixed_Forest", "Open_Water", "Pasture_Hay", "Perennial_Ice_Snow", "Shrub_Scrub",	
                    "Woody_Wetlands")
} else {
  print("Misidentified feature set")
  stop()
}

### Add texture means to model.indices
texture_means <- c("B4_Contrast_mean", "B4_Correlation_mean", "B4_Energy_mean", "B4_Entropy_mean",
                   "B4_Inertia_mean", "B4_Prominence_mean", "B4_Shade_mean", "B4_Variance_mean",
                   "B7_Contrast_mean", "B7_Correlation_mean", "B7_Energy_mean", "B7_Entropy_mean",
                   "B7_Inertia_mean", "B7_Prominence_mean", "B7_Shade_mean", "B7_Variance_mean")
texture_stdDevs <- c("B4_Contrast_stdDev", "B4_Correlation_stdDev", "B4_Energy_stdDev", "B4_Entropy_stdDev",
                     "B4_Inertia_stdDev", "B4_Prominence_stdDev", "B4_Shade_stdDev", "B4_Variance_stdDev",
                     "B7_Contrast_stdDev", "B7_Correlation_stdDev", "B7_Energy_stdDev", "B7_Entropy_stdDev",
                     "B7_Inertia_stdDev", "B7_Prominence_stdDev", "B7_Shade_stdDev", "B7_Variance_stdDev")
if (model.textures == TRUE) {
  model.indices <- append(model.indices, texture_means)
}
if (model.summaries == "mean_stdDev") { # also include standard deviations
  model.indices <- append(model.indices, texture_stdDevs)
}

# Select the indices and checklist ID
if (model == "modis") {
  model.variables <- model.indices
} else if (model == "nlcd" | model == "nlcd_binary") {
  model.variables <- c()
  for(i in 1:length(model.buffers)){
    model.variables <- append(model.variables, paste(model.indices, model.buffers[i], sep="_"))
  }
} else {
  model.variables <- c()
  for(i in 1:length(model.buffers)){
    for(j in 1:length(model.seasons)) {
      model.variables <- append(model.variables, paste(model.seasons[j], model.segmentation, model.indices, 
                                                       model.buffers[i], sep="_"))
    }
  }
}
print(model.variables)

###### Merge with counts and fold IDs ######
data <- merge(bird_data, environ_data[,c("Unique_Checklist_ID", model.variables)], by = "Unique_Checklist_ID")
data$Fold <- sb1$foldID 
data <- na.omit(data)

#########################################################
# 7. Run models
#########################################################
# On using classwt, strata and sampsize to compensate for extreme class imbalance

data$Occur <- as.factor(data$Occur) 

###### Prep to run models on each fold ###### 
# Extract fold IDs
folds <- data$Fold

# create a data.frame to store the predictions
testTable <- data[,c("Unique_Checklist_ID", "Fold", "Occur")]
testTable$Pred <- NA
num_predictors <- length(model.variables)

# Tune sampsize, RF defaults to using ~2/3 of data for each bootstrapped sample to build trees. 
# We'll use 2/3 of the presences and an equal sized sample of absences. We're doing this to compensate
# for the extreme class imbalance, rather than applying a downsampling scheme to address the imbalance. 
num_pres <- nrow(data[data$Occur==1,])
print("nrow(data[data$Occur==1,]):")
print(num_pres)
print(dim(data))
sampleSize = floor((2/3) * num_pres)

#full.variable.table <- data.frame(rn = NA, MeanDecreaseAccuracy = NA, k=NA)

###### Run models on each fold ###### 
predictors <- which(names(data) %in% model.variables)

for(k in seq_len(length(unique(folds)))){ 
  trainSet <- which(folds != k) # training set indices
  testSet <- which(folds == k) # testing set indices
  
  rf <- randomForest(data[trainSet, predictors], factor(data[trainSet, "Occur"]), ntree = 5000, 
                     sampsize=c(sampleSize, sampleSize))  #importance = TRUE
  testTable[testSet, "Pred"] <- predict(rf, data[testSet, predictors], type = "prob")[,2]  # predict the test set
  
  # OOB error
  #rf$err.rate[,1]
  #par(mfrow = c(2,1))
  #plot(rf$err.rate[,1], type = "l")
  #plot(rf)
  
  # Variable importance:
  #variable.table <- data.frame(importance(rf,type=1))
  #setDT(variable.table, keep.rownames = TRUE)[]
  #variable.table$k <- k
  
  #full.variable.table <- rbind(full.variable.table, variable.table)
  
  #jpeg(paste0("VariableImportance_TestFold", k, '.jpg'), width = 10, height = 6, units = "in", res = 300)
  #varImpPlot(rf)
  #dev.off()
}

# save predictions (testTable)
save(testTable, file = paste0(species, "_", model, "_", name, ".RData"))

# save model
#save(rf, file = paste0(species, "_", model, "_rf_model.RData"))


###### Assess Variable Importance ###### 
#full.variable.table = full.variable.table[-1,]
#names(full.variable.table)[1] <- "EnvironmentalVariable"
#env.variable.mean <- aggregate(MeanDecreaseAccuracy ~ EnvironmentalVariable, data = full.variable.table, FUN = mean)
#env.variable.sd <- aggregate(MeanDecreaseAccuracy ~ EnvironmentalVariable, data = full.variable.table, FUN = sd)
#variable.import.aggregated <- merge(env.variable.mean,env.variable.sd, by = "EnvironmentalVariable")
#names(variable.import.aggregated) <- c("EnvironmentalVariable", "MeanDecreaseAccuracy", "SDDecreaseAccuracy")
#write.csv(full.variable.table, paste0("VariableImportanceFullFolds_", model, ".csv"), row.names = FALSE)
#write.csv(variable.import.aggregated, paste0("VariableImportanceSummary_", model, ".csv"), row.names = FALSE)
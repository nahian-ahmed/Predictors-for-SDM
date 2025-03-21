#################################################################################################
# Given test tables from run_model.R, generate AUCs and CIs for models
#
# Laurel Hopkins 
# 11/2/2020
#################################################################################################

library(pROC)  
library(plotrix)
library(PRROC)
library(stringr)

Working.directory <-"C:\\Users\\Laurel\\Documents\\Oregon State\\Research\\ICB" 
setwd(Working.directory)


#########################################################
# 1. Calculate AUC and CIs given testTable
#########################################################
# list of species
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", "Hermit Warbler",  
#                  "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", "Swainson’s Thrush",  
#                  "Western Meadowlark", "Western Tanager", "Yellow Warbler", "Yellow-breasted Chat")
                  
models <- c("mean", "mean_stdDev", "mean_allSeasons", "mean_texture", "mean_stdDev_allSeasons",
            "mean_stdDev_texture", "mean_allSeasons_texture", "mean_stdDev_allSeasons_texture")
models <- c("modis", "nlcd", "nlcd_binary_all_classes")

for(model in models){
  get_metrics<- lapply(species_list, function(species) {
    print(species)
    files <- list.files(path=paste0("results/", model, "/", species), pattern="*.RData", full.names=TRUE, 
                        recursive=FALSE)
    
    eval <- lapply(files, function(x) {
      load(x)
      
      index <- unlist(strsplit(basename(x), ".RData"))[1]
      index <- unlist(strsplit(index, paste0(species, "_")))[2]
      
      #index <- substring(index, 6)
      #index <- str_remove(index, paste0("_", model))
      
      #index <- unlist(strsplit(index, "_"))[2]
      #index <- unlist(strsplit(index, paste0(species, "_")))[2]
      #index <- unlist(strsplit(index, paste0("_", model)))[1] # don't need this line if model is not appended to name
      print(index)
      
      # Calculate ROC
      #testTable <- testTable[testTable$Unique_Checklist_ID %in% test_ids,]  
      roc <- roc(response = as.factor(testTable$Occur), predictor = testTable$Pred, auc = TRUE, ci = TRUE)
      ci <- round(roc$ci, digit=5)
      data.frame (species=species, feature_set=index, AUC = ci[2], lowerCI = ci[1], upperCI = ci[3])  
    })
    do.call(rbind, eval)})
  
  metrics <- do.call(rbind, get_metrics)
  save(metrics, file = paste0("results/", model, "/", model, "_metrics.RData"))
}
# Alternatively, load metrics dataframe
#load(file = paste0("results/", model, "/", model, "_metrics.RData"))


#########################################################
# 2. Calculate AUC and CIs given testTable with a 
#    balanced number of presences and absences
#########################################################
# list of species
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Acorn Woodpecker", "Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
#                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
#                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "White-breasted Nuthatch", 
#                  "Yellow Warbler", "Yellow-breasted Chat")

models <- c("mean", "mean_stdDev", "mean_allSeasons", "mean_texture", "mean_stdDev_allSeasons",
            "mean_stdDev_texture", "mean_allSeasons_texture", "mean_stdDev_allSeasons_texture")
for(model in models){
  get_metrics<- lapply(species_list, function(species) {
    print(species)
    files <- list.files(path=paste0("results/", model, "/", species), pattern="*.RData", full.names=TRUE, 
                        recursive=FALSE)
    
    eval <- lapply(files, function(x) {
      load(x)
      
      index <- unlist(strsplit(basename(x), ".RData"))[1]
      index <- unlist(strsplit(index, paste0(species, "_")))[2]
      index <- unlist(strsplit(index, paste0("_", model)))[1] # don't need this line if model is not appeded to name
      print(index)
      
      # Randomly downsample absences to match number of presences
      testTable_pres <- testTable[testTable$Occur == 1,]
      testTable_abs <- testTable[testTable$Occur == 0,]
      num_pres <- nrow(testTable[testTable$Occur == 1,])
      abs_to_keep <- sample(1:nrow(testTable_abs), size=num_pres, replace=FALSE)
      testTable_abs <- testTable_abs[abs_to_keep,]
      testTable_balanced <- rbind(testTable_pres, testTable_abs)
      
      # Calculate ROC
      roc <- roc(response = as.factor(testTable_balanced$Occur), predictor = testTable_balanced$Pred, auc = TRUE, ci = TRUE)
      ci <- round(roc$ci, digit=5)
      data.frame (species=species, feature_set=index, AUC = ci[2], lowerCI = ci[1], upperCI = ci[3])  
    })
    do.call(rbind, eval)})
  
  metrics <- do.call(rbind, get_metrics)
  save(metrics, file = paste0("results/", model, "/", model, "_metrics_balanced.RData"))
}
# Alternatively, load metrics dataframe
#load(file = paste0("results/", model, "/", model, "_metrics.RData"))


#########################################################
# 3. Combine all metrics into a single data frame
#########################################################
models <- c("mean", "mean_stdDev", "mean_allSeasons", "mean_texture", "mean_stdDev_allSeasons",
            "mean_stdDev_texture", "mean_allSeasons_texture", "mean_stdDev_allSeasons_texture", 
            "modis", "nlcd", "nlcd_binary_all_classes")

all_metrics <- data.frame()
for(model in models){
  load(paste0("results/", model, "/", model, "_metrics.RData"))
  metrics$model <- model
  if(nrow(all_metrics)==0){
    all_metrics <- metrics
  } else {
    all_metrics <- rbind(all_metrics, metrics)
  }
}
save(all_metrics, file = paste0("results/all_metrics.RData"))


##### Do the same for the PR metrics
models <- c("mean", "modis", "nlcd", "nlcd_2400", "nlcd_binary", "nlcd_binary_top3")
all_pr_metrics <- data.frame()
for(model in models){
  load(paste0("results/", model, "/", model, "_pr_metrics.RData"))
  pr_metrics$model <- model
  if(nrow(all_pr_metrics)==0){
    all_pr_metrics <- pr_metrics
  } else {
    all_pr_metrics <- rbind(all_pr_metrics, pr_metrics)
  }
}
save(all_pr_metrics, file = paste0("results/all_pr_metrics.RData"))


#########################################################
# 4. Plot all index sets for a given species
#########################################################
metrics <- all_metrics
species <- "Ash-throated Flycatcher"

# extract parameters
AUC <- metrics[metrics$species==species, "AUC"]
uiw <- metrics[metrics$species==species, "upperCI"] - metrics[metrics$species==species, "AUC"]
liw <- metrics[metrics$species==species, "AUC"] - metrics[metrics$species==species, "lowerCI"]

features <- unlist(strsplit(as.character(metrics[metrics$species==species,"feature_set"]), "_600_vi"))

# make plot
par(mar=c(12.8,4.1,3,2.1)) # 10.2 for smaller plots, 12.8 for larger
plot(AUC, ylim = c(0.4, 1), xlab="", ylab="AUC", axes = FALSE, main = species,
     #plot(AUC, ylim = c(0.4, 1), xlab="", ylab="AUC", axes = FALSE, main =  paste(model, "T2V w/out SDM loss \n(Frag species)", sep=" - "),
     cex.lab=1.3, cex.main=1.5, cex.sub=1.3)
plotCI(AUC,y=NULL, uiw=uiw, liw=liw, err="y", pch=20, slty=3, scol = "black", add=TRUE, xlab="")
axis(side=2, cex.axis=1.3)
axis(side=1,at=1:length(features), label=features, las=2, cex.axis=1.3)  


#########################################################
# 5. Variable Importance
#########################################################
library(caret)
library(randomForest)

model <- "raw_bands_vi"
buffer <- "600"

species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", "Hermit Warbler", 
                  "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler",  
                  "Yellow-breasted Chat")

png(paste0("results/", model, "/Plots/", "Raw bands_VI.png"), 1500, 2000)
par(mfrow=c(4, 4)) #, xpd=NA)
for(species in species_list) {
  load(paste0("results/", model, "/", species, "/rf_fold10_", species, "_", model, "_", buffer, ".RData"))
  #load("results/tasseled_cap_vi/Western Tanager/rf_fold10_Western Tanager_tassled_cap_vi_600.RData")
  rows <- rownames(rf$importance)
  new_rows <- unlist(lapply(rows, function(x) unlist(strsplit(x, "_"))[[3]]))
  row.names(rf$importance) <- new_rows
  
  #varImp(rf)
  #importance(rf)
  
  #png(paste0("results/", model, "/Plots/", species, "_", model, "_", buffer, ".png"), 390, 600)
  if (species=="Ash-throated Flycatcher") {
    species <- "Ash-throated\nFlycatcher"
  }
  varImpPlot(rf, main = species, type = 1, cex=2.1)  # type 1: Accuracy, type 2: Gini
  #dev.off()
}
dev.off()

##### Partial Dependency Plots #####
species <- "Ash-throated Flycatcher"
model <- "tasseled_cap_vi"
model.indices <- "tasseled_cap"

model.buffers <- 2400 #c(75, 600, 2400)
model.seasons <- "summer"
model.segmentation <- "nbr"

# load RF model
load(paste0("results/", model, "/", species, "/rf_fold10_", species, "_", model, ".RData"))

# load environmental data 
environmental_file <- "reformatted_bird_buffer_samples_07_07_2020_20_06_30.csv"
environ_data <- read.csv(file=environmental_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)

# load species counts
species_file <- paste0("SpeciesCounts\\", species, "_Counts_062320.csv")
bird_data <- read.csv(file=species_file, header=TRUE, sep=",", dec=".", stringsAsFactors=FALSE, fill=TRUE)

if (model.indices == "raw_bands") {
  model.indices <- c("B1_mean", "B2_mean", "B3_mean", "B4_mean", "B5_mean", "B7_mean")
} else if (model.indices == "tasseled_cap") {
  model.indices <- c("TCA_mean", "TCB_mean", "TCG_mean", "TCW_mean")
}

# extract model variables
model.variables <- c()
for(i in 1:length(model.buffers)){
  for(j in 1:length(model.seasons)) {
    model.variables <- append(model.variables, paste(model.seasons[j], model.segmentation, model.indices, 
                                                     model.buffers[i], sep="_"))
  }
}

#merge envrionmental data and species counts
data <- merge(bird_data, environ_data[,c("Unique_Checklist_ID", model.variables)], by = "Unique_Checklist_ID")

#load(paste0("results/", model, "/", species, "/", species, "_", model, ".RData"))  # testTable
partialPlot(rf, data[,5:8], x.var=model.variables[4], which.class= "1", plot=TRUE, add=FALSE, 
            main=paste0("Partial dependence on ", model.variables[4]))
            #n.pt = min(length(unique(pred.data[, xname])), 51),
            #rug = TRUE, xlab=deparse(substitute(x.var)), ylab="",
            #main=paste("Partial Dependence on", deparse(substitute(x.var))),
            #...)
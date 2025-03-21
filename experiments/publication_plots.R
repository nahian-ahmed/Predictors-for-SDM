#################################################################################################
# Publication ready plots. Run this after generating all_metrics.RData with 
# evaluation.R
#
# Laurel Hopkins 
# 11/2/2020
#################################################################################################

Working.directory <-"C:\\Users\\Laurel\\Documents\\Oregon State\\Research\\ICB" 
setwd(Working.directory)

#########################################################
# 1. Prep data
#########################################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Use script evaluation.R to generate all_metrics.RData 
load(file = paste0("results/all_metrics.RData"))

models <- c("mean", "mean_allSeasons", "mean_stdDev", "mean_texture", "mean_stdDev_allSeasons",
            "mean_allSeasons_texture", "mean_stdDev_texture", "mean_stdDev_allSeasons_texture")

spectral_predictors <- c("evi", "msavi", "savi", "nbr", "nbr_2", "ndmi", "ndsi", "ndvi", "raw_bands", 
                        "tasseled_cap")
species_list <- c("Ash-throated Flycatcher")
#species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
#                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
#                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
#                  "Yellow-breasted Chat")

data_df <- all_metrics[all_metrics$model %in% models & all_metrics$feature_set %in% spectral_predictors &
                         all_metrics$species %in% species_list,]

# Rename models and feature sets for plotting
data_df$model[data_df$model == "mean"] <- "Summer means"
data_df$model[data_df$model == "mean_allSeasons"] <- "Sp/Su/Fa means"
data_df$model[data_df$model == "mean_stdDev"] <- "Summer means & SDs"
data_df$model[data_df$model == "mean_texture"] <- "Summer means & texture"
data_df$model[data_df$model == "mean_stdDev_allSeasons"] <- "Sp/Su/Fa means & SDs"
data_df$model[data_df$model == "mean_allSeasons_texture"] <- "Sp/Su/Fa means & texture"
data_df$model[data_df$model == "mean_stdDev_texture"] <- "Summer means, SDs & texture"
data_df$model[data_df$model == "mean_stdDev_allSeasons_texture"] <- "Sp/Su/Fa means, SDs & texture"
data_df$feature_set <- as.character(data_df$feature_set)
data_df$feature_set[data_df$feature_set == "evi"] <- "EVI"
data_df$feature_set[data_df$feature_set == "msavi"] <- "MSAVI"
data_df$feature_set[data_df$feature_set == "savi"] <- "SAVI"
data_df$feature_set[data_df$feature_set == "nbr"] <- "NBR"
data_df$feature_set[data_df$feature_set == "nbr_2"] <- "NBR2"
data_df$feature_set[data_df$feature_set == "ndmi"] <- "NDMI"
data_df$feature_set[data_df$feature_set == "ndsi"] <- "NDSI"
data_df$feature_set[data_df$feature_set == "ndvi"] <- "NDVI"
data_df$feature_set[data_df$feature_set == "raw_bands"] <- "Raw bands"
data_df$feature_set[data_df$feature_set == "tasseled_cap"] <- "Tasseled Cap"

# order of spectral predictors for plots
order_to_plot <- c("EVI", "MSAVI", "SAVI", "NBR", "NBR2", "NDMI", "NDSI", "NDVI", "Raw bands", "Tasseled Cap")
data_df$feature_set <- factor(data_df$feature_set, levels = order_to_plot)
levels(data_df$feature_set)

# reorder models for plotting 
model_neworder <- c("Summer means", "Sp/Su/Fa means", "Summer means & SDs", "Sp/Su/Fa means & SDs", 
                    "Summer means & texture", "Sp/Su/Fa means & texture", "Summer means, SDs & texture",
                    "Sp/Su/Fa means, SDs & texture")
data_df <- arrange(transform(data_df, model=factor(model, levels=model_neworder)), model)

# Species means for the summer means model
for (species in species_list) {
  species_mean <- mean(data_df[data_df$species == species, "AUC"])
  print(paste0(species, " mean: ", species_mean))
}


#########################################################
# 2. Figure 2: AUC and CIs for all spectral predictors  
#    for each species
#########################################################
## Plots
# make panel plot of all species for resubmission
species_to_plot <- species_list
to_plot <- data_df[data_df$model=="Summer means" & data_df$species %in% species_to_plot,]

### labeled
gp <- ggplot(to_plot, aes(x=feature_set, y=AUC)) + 
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.5, size=0.2) + 
  geom_line() + 
  geom_point(size=0.3) + 
  facet_wrap(~ species, ncol=4) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust=1), 
        axis.ticks.x=element_blank(), text=element_text(size=8)) 
gp

# save plot
ggsave("Fig2.tiff", gp, dpi=1200, width = 119, height = 119, units = "mm")
ggsave("Fig2.png", gp, dpi=1200, width = 119, height = 119, units = "mm")


# make two panel plots of species, 4 for paper, remaining in supplemental material
species_to_plot <- c("Ash-throated Flycatcher", "Hermit Thrush", "Sagebrush Sparrow", "Western Tanager")
to_plot <- data_df[data_df$model=="Summer means" & data_df$species %in% species_to_plot,]

### color coded
gp <- ggplot(to_plot, aes(x=feature_set, y=AUC, colour=feature_set, group=feature_set)) + 
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.8) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~ species, ncol=4) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        text = element_text(size = 22)) + 
  scale_color_brewer(name = "Spectral\nPredictors", palette="Paired")    
#scale_color_manual(name = "Spectral\nPredictors", values = clrs.hcl(10))
gp

### labeled
gp <- ggplot(to_plot, aes(x=feature_set, y=AUC)) + 
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.5, size=0.2) + 
  geom_line() + 
  geom_point(size=0.3) + 
  facet_wrap(~ species, ncol=4) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust=1), 
        axis.ticks.x=element_blank(), text=element_text(size=8)) 
gp

# save plot
ggsave("Fig2.tiff", gp, dpi=1200, width = 119, height = 43, units = "mm")
ggsave("Fig2.png", gp, dpi=1200, width = 119, height = 43, units = "mm")


# make panel plots for species to go in supplemental materials
other_species <- setdiff(species_list, species_to_plot)

to_plot <- data_df[data_df$model=="Summer means" & data_df$species %in% other_species,]

### color coded
gp_color <- ggplot(to_plot, aes(x=feature_set, y=AUC, colour=feature_set, group=feature_set)) + 
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.8) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~ species, ncol=3) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        text = element_text(size = 16.8)) + 
  scale_color_brewer(name = "Spectral\nPredictors", palette="Paired")    
#scale_color_manual(name = "Spectral\nPredictors", values = clrs.hcl(10))
gp_color

### labeled
gp <- ggplot(to_plot, aes(x=feature_set, y=AUC)) + 
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.5, size=0.2) + 
  geom_line() + 
  geom_point(size=0.3) + 
  facet_wrap(~ species, ncol=3) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust=1), 
        axis.ticks.x=element_blank(), text=element_text(size=8)) 
gp

# save plot
ggsave("FigS1.tiff", gp, dpi=1200, width = 119, height = 110, units = "mm")
ggsave("FigS1.png", gp, dpi=1200, width = 119, height = 110, units = "mm")


#########################################################
# 3. Figure 3. box plots of the mean AUCs of the   
#     predictors for the different summary methods 
#########################################################
library(plyr)
library(MASS)
library(colorspace)
library(Hmisc)
library(scales)
library(RColorBrewer)
# reorder models for plotting 
model_neworder <- c("Summer means", "Sp/Su/Fa means", "Summer means & SDs", "Sp/Su/Fa means & SDs", 
                    "Summer means & texture", "Sp/Su/Fa means & texture", "Summer means, SDs & texture",
                    "Sp/Su/Fa means, SDs & texture")
data_df <- arrange(transform(data_df, model=factor(model, levels=model_neworder)), model)

### Color plots ###
num_levels <- length(unique(data_df$feature_set))

# color palette
clrs.spec <- colorRampPalette(rev(brewer.pal(num_levels, "Spectral")))

clrs.hcl <- function(n) {
  hcl(h = seq(230, 0, length.out = n), 
      c = 60, l = seq(10, 90, length.out = n), 
      fixup = TRUE)
}

# function to plot a colour palette
pal <- function(col, border = "transparent", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

bw.ggplot <- ggplot(data_df, aes(x = feature_set, y = AUC)) +
  scale_x_discrete(limits = order_to_plot)   

g.bw <- bw.ggplot + 
  #scale_colour_manual(name="Spectral\nPredictors:", values=vars) +
  #theme(legend.direction="horizontal",legend.position="top", legend.box = "vertical") +
  facet_wrap(~ model, ncol=2) +
  llply(unique(data_df$feature_set), 
        function(i) geom_boxplot(#color="black",
          fill = clrs.hcl(num_levels)[i],
          data = subset(data_df, feature_set == i)))  +
  theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust=1), 
                     text = element_text(size = 12.5)) 
# show plot
print(g.bw)# save plot

ggsave("Fig3.tiff", g.bw, dpi=1200, width = 119, height = 195, units = "mm")
ggsave("Fig3.png", g.bw, dpi=1200, width = 119, height = 195, units = "mm")


#########################################################
# 3. Figure 4: Raw bands comparison and
#     NDVI
#########################################################
library(plyr)
library(MASS)
library(colorspace)
library(Hmisc)
library(scales)
library(RColorBrewer)

spectral_predictors_mini <- c("Raw bands", "NDVI")

data_df_mini <- data_df[data_df$feature_set %in% spectral_predictors_mini,]


### Plots
# color palette
num_levels <- length(unique(data_df_mini$model))
clrs.spec <- colorRampPalette(rev(brewer.pal(num_levels, "Spectral")))
clrs.hcl <- function(n) {
  hcl(h = seq(230, 0, length.out = n), 
      c = 60, l = seq(10, 90, length.out = n), 
      fixup = TRUE)
}

# function to plot a colour palette
pal <- function(col, border = "transparent", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}


# make boxplots 
bw.ggplot <- ggplot(data_df_mini, aes(x = model, y = AUC)) +
  scale_x_discrete(limits = model_neworder)   

g.bw <- bw.ggplot +
  facet_wrap(~ feature_set) +
  llply(unique(data_df_mini$model), 
        function(i) geom_boxplot(fill = clrs.hcl(num_levels)[i],
                                 data = subset(data_df_mini, model == i),
                                 width=0.5, size=0.2, outlier.size=0.3))  +
  theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust=1), #axis.ticks.x=element_blank(),
                     text = element_text(size = 8)) 
# show plot
g.bw
# save plot
ggsave("Fig4.tiff", g.bw, dpi=1200, width = 90, height = 60, units = "mm")
ggsave("Fig4.png", g.bw, dpi=1200, width = 90, height = 60, units = "mm")



#########################################################
# 5. Figure 4: Classified comparison
#########################################################
library(plyr)
library(MASS)
library(colorspace)
library(Hmisc)
library(scales)
library(RColorBrewer)

models <- c("mean", "modis", "nlcd", "nlcd_binary_all_classes")
species_list <- c("Ash-throated Flycatcher", "Gray Flycatcher", "Hermit Thrush", 
                  "Hermit Warbler", "Pacific Wren", "Sage Thrasher", "Sagebrush Sparrow", "Savannah Sparrow", 
                  "Swainson’s Thrush", "Western Meadowlark", "Western Tanager", "Yellow Warbler", 
                  "Yellow-breasted Chat")

data_df <- all_metrics[all_metrics$model %in% models & all_metrics$species %in% species_list,]
spectral_predictors <- c("raw_bands", "modis", "nlcd", "nlcd_binary_all_classes")
data_df <- data_df[data_df$feature_set %in% spectral_predictors,]

# remane models
data_df$feature_set <- as.character(data_df$feature_set)
data_df$feature_set[data_df$feature_set == "raw_bands"] <- "Raw bands"
data_df$feature_set[data_df$feature_set == "modis"] <- "MODIS"
data_df$feature_set[data_df$feature_set == "nlcd"] <- "NLCD"
data_df$feature_set[data_df$feature_set == "nlcd_binary_all_classes"] <- "NLCD\ndiscretized"


# order of spectral predictor box plots
order_to_plot <- c("Raw bands", "MODIS", "NLCD", "NLCD\ndiscretized")
data_df$feature_set <- factor(data_df$feature_set, levels = order_to_plot)
levels(data_df$feature_set)

### make plots
# color palette
clrs.spec <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
clrs.hcl <- function(n) {
  hcl(h = seq(230, 0, length.out = n), 
      c = 60, l = seq(10, 90, length.out = n), 
      fixup = TRUE)
}

# function to plot a colour palette
pal <- function(col, border = "transparent", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

vars <- c("Raw bands"=clrs.hcl(10)[1], "MODIS"=clrs.hcl(10)[5], "NLCD"=clrs.hcl(10)[2], 
          "NLCD discretized"=clrs.hcl(10)[9])

# make boxplots 
bw.ggplot <- ggplot(data_df, aes(x = feature_set, y = AUC)) +
  scale_x_discrete(limits = order_to_plot)   

g.bw <- bw.ggplot + 
  llply(unique(data_df$feature_set), 
        function(i) geom_boxplot(#color="black",
          fill = vars[i], #clrs.hcl(10)[i],
          data = subset(data_df, feature_set == i),
          width=0.5, size=0.2, outlier.size=0.3))  +
  theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust=1), #axis.ticks.x=element_blank(),
                     text = element_text(size = 6)) 
# show plot
print(g.bw)
# save plot
ggsave("Fig5.tiff", g.bw, dpi=1200, width = 45, height = 40, units = "mm")
ggsave("Fig5.png", g.bw, dpi=1200, width = 45, height = 40, units = "mm")


#########################################################
# 6. Plot all index sets for a given species
#########################################################
library(ggplot2)

model <- "tasseled_cap_indep"
load(file = paste0("results/", model, "/", model, "_metrics.RData")) 

metrics$feature_set <- as.character(metrics$feature_set)
if (model == "raw_bands_indep"){
  metrics$feature_set[metrics$feature_set == "B1_600_vi"] <- "B1"
  metrics$feature_set[metrics$feature_set == "B2_600_vi"] <- "B2"
  metrics$feature_set[metrics$feature_set == "B3_600_vi"] <- "B3"
  metrics$feature_set[metrics$feature_set == "B4_600_vi"] <- "B4"
  metrics$feature_set[metrics$feature_set == "B5_600_vi"] <- "B5"
  metrics$feature_set[metrics$feature_set == "B7_600_vi"] <- "B7"
} else if (model == "tasseled_cap_indep") {
  metrics$species <- as.character(metrics$species)
  metrics$species[metrics$species == "Ash-throated Flycatcher"] <- "Ash-throated\nFlycatcher"
  metrics$feature_set[metrics$feature_set == "TCA_600_vi"] <- "TCA"
  metrics$feature_set[metrics$feature_set == "TCB_600_vi"] <- "TCB"
  metrics$feature_set[metrics$feature_set == "TCG_600_vi"] <- "TCG"
  metrics$feature_set[metrics$feature_set == "TCW_600_vi"] <- "TCW"
}


### labeled
gp <- ggplot(metrics, aes(x=feature_set, y=AUC)) + 
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.8) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~ species, ncol=3) +
  theme_bw() + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust=1), 
        axis.ticks.x=element_blank(), text=element_text(size=24)) 
gp


# save plot
ggsave(".png", gp, dpi=320,
       width = 20, height = 12, units = "cm")

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
# 7. Variable Importance
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


##### Partial Dependency Plots - not included in paper##### 
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


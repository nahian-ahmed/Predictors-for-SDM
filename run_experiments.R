#################################################################################################
# Predictors for Habitat Suitability Models
# Summary: Compare habitat suitability models based on different satellite-derived indices
#          and embeddings from geospatial foundation models
#
# Nahian Ahmed
# 03/21/2025
#################################################################################################



# set working directory to current/source directory
library(here)
setwd(here::here())

f_path <- "results"
if (!dir.exists(f_path)) {
        dir.create(f_path)
}


run_experiments <- function(model_names, species_names, experiments=TRUE, summarize=TRUE, stat_test=TRUE, plot_results=TRUE){
    
    assign("model_names", model_names, envir = .GlobalEnv)
    assign("species_names", species_names, envir = .GlobalEnv)
    


    if (experiments){

        # Run experiments
        for (model_name in model_names){
            for (species_name in species_names){

                assign("model", model_name, envir = .GlobalEnv)
                assign("species", species_name, envir = .GlobalEnv)
                
                # Run experiments on specific model and species data
                source("experiments/run_model.R")

            }
        }
    }

    if (summarize){
        # Summarize results
        source("experiments/evaluation.R")
    }

    if (stat_test){
        # Run statistical tests
        source("experiments/statistical_analysis.R")
    }

    if (plot_results){
        # Plot results
        source("experiments/publication_plots.R")
    }
}


# Model options: "raw_bands", "tasseled_cap", "ndvi", "ndmi", "nbr", "nbr_2", "evi", "savi", "msavi", "ndsi"
model_names = c(
    "raw_bands",
    "tasseled_cap", 
    "ndvi", 
    "ndmi", 
    "nbr", 
    "nbr_2", 
    "evi", 
    "savi", 
    "msavi", 
    "ndsi",
    "embedding_fields"
)

species_names = c(
    "Acorn Woodpecker",
    "Ash-throated Flycatcher",
    "Gray Flycatcher",
    "Hermit Thrush",
    "Hermit Warbler",
    "Pacific Wren",
    "Sage Thrasher",
    "Sagebrush Sparrow",
    "Savannah Sparrow",
    "Swainson's Thrush",
    "Western Meadowlark",
    "Western Tanager",
    "White-breasted Nuthatch",
    "Yellow Warbler",
    "Yellow-breasted Chat"
)

run_experiments(
    model_names = model_names,
    species_names = species_names,
    experiments = TRUE,
    summarize = FALSE,
    stat_test = FALSE,
    plot_results = FALSE
)

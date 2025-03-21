# set working directory to current/source directory
library(here)
setwd(here::here())

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
    "ndsi"
)

species_names = c(
    "Ash-throated Flycatcher"
)

dirs_to_create = c(
    "results",
    "results/species_experiments",
    "results/species_experiments/raw",
    "results/species_experiments/raw/test_data_splits",
    "results/species_experiments/raw/clusterings",
    "results/species_experiments/raw/clustering_parameters",
    "results/species_experiments/raw/model_parameters",
    "results/species_experiments/raw/predictions",
    "results/species_experiments/raw/metrics",
    "results/species_experiments/summarized",
    "results/species_experiments/plotted",
    "results/species_experiments/plotted/plots",
    "results/species_experiments/plotted/maps"
)

for (f.path in dirs_to_create){
    if (!dir.exists(f.path)) {
        dir.create(f.path)
    }
}

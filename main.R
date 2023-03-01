"Instructions for use:
1. Install / load the libraries
2. Install the user-defined functions from: st_2_upload_data, st_4_regrf, and st_7_plotpredictions
3, Define parameters (taxalevel, Rarefaction level, environmental_variables, SampleSize)
4. Run the code via the wrapper functions.
  4.1 Run the Initial Data Analysis. Cleans and evaluates the data.
  4.2 Note that some uploading functions are to access the data from a remote drive
      and are not utilised here (data transferred to local drive from where it is uploaded)
  4.3 Run the training of the Random Forest algorithm and make plots.
"

#________________________________________________________________________________________________________________________________________#
## 1. Libraries ##
" Load the libraries required for analysis."

# Check if all required packages are installed, if not download them
packages = c("reshape2", "data.table", "splitstackshape", "strex", "vegan", "latticeExtra",
             "stringr", "R2wd", "gtools", "grid", "rstudioapi", "R2wd", "lattice",
             "randomForest", "tidyr", "ggplot2", "gridExtra", "dplyr",
             "ggpubr", "ecodist", "ggrepel", "RColorBrewer", "squash", "caret", "vegan", "plyr",
             "furrr")
new_packages = packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all the packages
lapply(packages, library, character.only = TRUE)

#________________________________________________________________________________________________________________________________________#
## 2. Load functions ##
" Load the functions written in the upload_data, initial_data_analysis, and random_forest files."
# Set working directory to source file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load the functions into R
source("upload_data.r")
source("initial_data_analysis.r")
source("random_forest.r")

#________________________________________________________________________________________________________________________________________#
## 3. Parameter definitions ##
" Define the parameters required. taxalevel is needed to isolate until what specificity we want to identify the data.
rarefaction_rate is needed to to set the random subsample size to give the expected species richness.
environmental_variables is needed to include extra variables to train the RF model"
## For the data set (data) creation ##
# Specify which sources to keep in the data set
# sources = c("Thorsten", "Shraveena", "MOWI", "MeioMetBar", "SSF", "VAW")
sources = c("Thorsten", "Shraveena", "MOWI", "MeioMetBar")

# Define the taxon specificity of the data set
# taxalevel can be "Kingdom", "Phylum", "Class", Order", or "Family"
taxalevel = "Family"

# Set sample size for the rarefaction_rate function
# 2500 excludes 3 Thorsten samples + all Shraveena, 10816 is the MeioM limit and excludes most of Thorsten's samples
# When set to NA no rarefaction is performed on the data set
"Check the number of read counts per sample for Rarefaction."
rarefaction_rate = 1170 # Determined based on the plot of the rarefaction function

## For the Random Forest algorithm ##
# Determine which environmental conditions should be used to train the RF model
# environmental_variables can be c() or any combination of c("MDTheta","OrganicCarbon","Month", "Distance")
environmental_variables = c()

#________________________________________________________________________________________________________________________________________#
## 4.1 Initial Data Analysis ##
" Check whether there are missing values and if there are any unrealistic values based on
the background knowledge on this topic."

"When data_upload is run we take the 16S sequence data of the pre-specified sources and
merge it with the known metadata from these sources. This creates one data set, MDS4, which
contains all the data obtained from all these sample sites."

## Rarefaction
'Generate the rarefaction curves for the unrarefied data.'
MDS4 = data_upload(sources, taxalevel, rarefaction_rate = NA)
rarefaction(MDS4)

## Check for negative values in the data set ##
"There are no negative values in the data set."
sum(MDS4 < 0)

## Count the number of observations per site ##
observations = data.frame(table(MDS4$Site))
colnames(observations) = c("Site", "Frequency")
observations

## Display hist, boxplots and barplots of IQI value per site coloured by source ##
# May be interesting to colour by Source
"In the boxplot, the bar represents the median of the dataset"
iqi_source(MDS4)

# Load the rarefied data
MDS4 = data_upload(sources, taxalevel, rarefaction_rate=rarefaction_rate)

## Alpha diversity ##
"Generate boxplots of Chao1 and Shannon alpha diveristy indices."
alpha_diversity(MDS4)

## Beta diversity ##
"Generate a Bray-Curtis dissimilarity matrix (that uses a rarefaction rate).
Generate a PCoa based on the Bray-Curtis matrix."
beta_diversity(MDS4)

#________________________________________________________________________________________________________________________________________#
## 4.3 Random Forest training ##
" Generate a Random Forest data set with the SampleIDs as row names. Uses the data set to run multiple RF models.

Per default, Cross-Validation with an 80% train and 20% test division based on the site is used.
This can be set equal to false in the Build functions. "

## Generate Random Forest model
# Run the randomForests in parallel
future::plan(multisession)

# Test the different rarefaction rates
# normalization = c(1170, 5244, 13493, 14987, 50384, NA)
normalization = c(1170, 5244, 13493, NA)

# Save sampleIDs that can reach 13493 read counts as "samples"
rarefy_13493 = data_upload(sources, taxalevel, rarefaction_rate = 13493)
samples = rarefy_13493$SampleID

# Ensure the training uses the same data set
rand_train = sample(1:nrow(rarefy_13493), size = nrow(rarefy_13493) * 0.8)

source_counts = future_map_dfr(normalization, function(normalization){
    for (i in 1:length(normalization)){rarefaction_rate = as.numeric(normalization[i])}

    # Upload data
    MDS4 = data_upload(sources, taxalevel, rarefaction_rate, samples = samples) # Set samples to FALSE if you want to rarefy over all samples

    # Count samples per source
    source_counts = data.frame(Source = unique(MDS4$Source),
                               normalization = normalization[i])
    source_counts = full_join(source_counts, count(MDS4, "Source"), by = "Source")

    # Generate the data for RF analysis
    RF_Data = prepare_rf_data(MDS4, environmental_variables)
 
    # Make the training and testing sets
    train = RF_Data[rand_train, ]
    test = MDS4[-rand_train, ]

    # Build the RF algorithmc
    RF = build_rf(train, sources, taxalevel, rarefaction_rate)

    # Make and plot the predictions on test data
    plot_iqi_prediction(test, RF$RF_final_reduced, sources, taxalevel, rarefaction_rate)

    return(source_counts)
  })

# Rename and save counts per source
source_counts = source_counts %>% pivot_wider(names_from = normalization,
                                              values_from = freq)
write.csv(source_counts, file = paste("figures/source_counts_sources_", toString(sources), ".csv"), row.names = FALSE)

#________________________________________________________________________________________________________________________________________#
## 4.3 Random Forest testing ##
" Load the top models and make predictions on the new data to check accuracy."

# Generate function to generate rarefied new test data and run plots
test_newdata = function(rarefaction_rate, model, samples = FALSE){
  sources = c("SSF", "VAW")
  data = data_upload(sources, taxalevel, rarefaction_rate, samples = samples)
  plot_iqi_prediction(data, model, sources, taxalevel, rarefaction_rate)
}

# Select the required samples
rarefy_13493_test = data_upload(sources, taxalevel, rarefaction_rate = 13493)
samples_test = rarefy_13493_test$SampleID

# Load all sources rarefaction 5244 model and make predictions
load("figures/randomForest_models/Rarefy_1170_Sources_Thorsten, MeioMetBar, Shraveena, MOWI_TaxaLevel_Family_Scaling_FALSE.Rdata")
test_newdata(NA, RF_final_reduced)

# Load all sources except Thorsten rarefaction 13493 model and make predictions
load("figures/randomForest_models/Rarefy_13493_Sources_MeioMetBar, Shraveena, MOWI_TaxaLevel_Family_Scaling_FALSE.Rdata")
test_newdata(13493, RF_final_reduced)

# Master Thesis 

The goal of this master thesis is to predict IQI values around salmon farms in Scotland. This is done by training the Random Forest algorithm on bacterial eDNA read counts using known measured IQI values as the response. To perform this analysis, one data set containing the sample names as rows and bacterial taxa with their occurence as columns was provided. Additionally, a data set with sample names and information about the sample such as location, IQI value, and environmental conditions was provided. Since the data is confidential it is not added to this repository.

**`main.R`** Contains the main script which loads the functions in upload_data, initial_data_anlysis, and random_forest. It then runs a script to perform the intial data analysis and train the Random Forest algorithm under different settings, where each setting is run over a different computer core to reduced computation time. In the end, the final models are loaded into the global environment and used to generate prediction plots on independent data sets.

**`upload_data`** Contains a main function that can be called to upload the bacterial reads data set and the metadata. It then combines these data sets based on the sample names. After which the mean is taken over replicate measurements and the data is rarefied, anonymised, and cleaned.

**`initial_data_analysis`** Contains functions to perform the intial data analysis such as plotting the mean IQI values over the different salmon farms, generating rarefaction plots, calculating alpha and beta diversity per site, and generating a PCoA based on the beta diversity. 

**`random_forest`** Contains functions to generate and save a tuned Random Forest model based on a 5-fold cross-validation excluding sites. Additionally, there is a function to generate IQI predictions based on the optimal model and plot the actual IQI values against the predicted IQI values.

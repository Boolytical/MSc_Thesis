#Function: prepare_rf_data
#Takes:MDS5 (SampleID by collated taxa)
#Returns:list of two dataframe, one taxa, one Env metadata.
#Calls: NA
#Job: Finds 'Site' (which identifies the start of the metadata) in MDS5,
#and then splits the data into taxa and env dataframes and rarefies
prepare_rf_data = function(MDS5, environmental_variables){
  message("Function: 'prepare_rf_data' start")
  # Reorder columns in the dataframe for RF analysis
  # Make a data frame of taxa and environmental columns for analysis
  RRFTaxa = MDS5[ , grepl("Bacteria", names(MDS5))] # Isolate taxa columns
  environmental_variables = ifelse(is.vector(environmental_variables), append(environmental_variables, "Site"), c("Site"))
  EnvCols = as.data.frame(MDS5[, environmental_variables])
  colnames(EnvCols) = environmental_variables

  # Fill in Unknown for NAs, partially complete environmental columns can be used now
  if(ncol(EnvCols)>0){ EnvCols[is.na(EnvCols), ] = "Unknown"}

  # Bind Taxa and Environment columns
  RF_Data = cbind(RRFTaxa, EnvCols)
  rownames(RF_Data) = MDS5$SampleID # Add SampleID as row names, allows easy metadata merge later on
  RF_Data$IQI = MDS5$IQI

  message("Function: 'prepare_rf_data' end")
  return(RF_Data)
}

#Function: custom_rf
#Takes: parameters for comparison, classificaion/regression type
#Returns: NA
#Calls: randomForest
#Job: Compares values of the mtry, maxnodes, and ntree sets
custom_rf = function(type = 'Regression'){
  # Make custom RF for comparing mtry and ntrees
  custom_rf = list(type = type, library = "randomForest", loop = NULL)

  custom_rf$parameters = data.frame(parameter = c("mtry", "maxnodes"), class = rep("numeric", 2, label = c("mtry", "maxnodes")))
  custom_rf$grid = function(x, y, len = NULL, search = "grid") {}
  custom_rf$fit = function(x, y, wts, param, lev, last, weights, classProbs, ...){randomForest(x, y, mtry = param$mtry, maxnodes = param$maxnodes, ...)}
  custom_rf$predict = function(modelFit, newdata, preProc = NULL, submodels = NULL){predict(modelFit, newdata)}
  custom_rf$prob = function(modelFit, newdata, preProc = NULL, submodels = NULL){predict(modelFit, newdata, type = "prob")}
  custom_rf$sort = function(x){x[order(x[, 1]),]}
  custom_rf$levels = function(x){x$classes}

  return(custom_rf)
}

#Function: build_rf
#Takes: rarefaction_rate
#Returns: NA
#Job: Generate a random forest model.
build_rf = function(RF_Data, sources, taxalevel, rarefaction_rate, CV = "Site"){
  message("Function: 'build_rf' start")
  set.seed(123) # Makes running the model reproducible

  if (CV == "Site"){
    # Set the 5-fold CV using Grid search, excluding entire sites
    trControl_grid = trainControl(method = "cv",
                                  index = groupKFold(RF_Data$Site, k = 5),
                                  number = 5,
                                  search = "grid")
  } else if (CV == "Sample"){
    # Set the 5-fold CV using Grid search
    trControl_grid = trainControl(method = "cv",
                                  number = 5,
                                  search = "grid")
  }
  RF_Data$Site = NULL

  # Generate the default RF model
  RF_default = train(IQI ~ .,
                     data = RF_Data,
                     method = "rf",
                     metric = "RMSE",
                     trControl = trControl_grid,
                     importance = TRUE)
  message(paste("RMSEs of default model = ", paste(RF_default$results$RMSE, collapse = ", ")))

  ## Tune the mtry value
  # Do a rough grid search on equally spaced mtry values
  tuneGrid = expand.grid(.mtry = round(seq(10, ncol(RF_Data), by = 30)))
  RF_grid = train(IQI ~ .,
                  data = RF_Data,
                  method = "rf",
                  metric = "RMSE",
                  tuneGrid = tuneGrid,
                  trControl = trControl_grid)

  # Select mtry that gave smallest RMSE value in default or grid search
  mtry = ifelse(RF_default$results[RF_default$results$mtry == RF_default$bestTune$mtry, "RMSE"] <
                  RF_grid$results[RF_grid$results$mtry == RF_grid$bestTune$mtry, "RMSE"], RF_default$bestTune$mtry, RF_grid$bestTune$mtry)
  message(paste("Best mtry value =", mtry))


  ## Tune the maxnodes value
  # randomForest automatically calculates maxnodes (# terminal nodes) to maxnodes (# all nodes)
  # using 2 * maxnodes - 1, so we have to convert back to this before optimizing
  maxnodes_default = round(dim(RF_default$finalModel$forest$nodestatus)[1] + 1) / 2
  store_maxnodes = list()
  tuneGrid = expand.grid(.mtry = mtry)
  for (maxnodes in seq(1, maxnodes_default, by = 10)){
    RF_maxnodes = train(IQI ~ .,
                        data = RF_Data,
                        method = "rf",
                        metric = "RMSE",
                        tuneGrid = tuneGrid,
                        trControl = trControl_grid,
                        maxnodes = maxnodes)
    iteration = toString(maxnodes)
    store_maxnodes[[iteration]] = RF_maxnodes
  }
  results_maxnodes = resamples(store_maxnodes)
  summary(results_maxnodes)
  maxnodes = as.numeric(names(which.min(summary(results_maxnodes)$statistics$RMSE[, 4])))
  message(paste("Best maxnodes value =", maxnodes))

  ## Tune the ntrees value
  store_maxtrees = list()
  for (ntree in c(500, 1000, 1500, 2000)){
    RF_maxtrees = train(IQI ~ .,
                        data = RF_Data,
                        method = "rf",
                        metric = "RMSE",
                        tuneGrid = tuneGrid,
                        trControl = trControl_grid,
                        maxnodes = maxnodes,
                        ntree = ntree)
    iteration = toString(ntree)
    store_maxtrees[[iteration]] = RF_maxtrees
  }
  results_tree = resamples(store_maxtrees)
  summary(results_tree)
  ntree = as.numeric(names(which.min(summary(results_tree)$statistics$RMSE[, 4])))
  message(paste("Best ntree value =", ntree))

  # Set the tuning grid
  message("Perform a grid search on mtry and maxnodes to find interaction between parameters")
  # Do a more accurate grid search around best mtry of the model; spans section between rough mtry of grid search
  tuneGrid = expand.grid(.mtry = seq(from = ifelse((mtry - 15) < 0, 0, mtry - 15),
                                     to = ifelse((mtry + 15) < ncol(RF_Data), mtry + 15, ncol(RF_Data)),
                                     by = 1),
                         .maxnodes = seq(from = ifelse((maxnodes - 10) < 0, 2, maxnodes - 10),
                                         to = ifelse((maxnodes + 10) > maxnodes_default, maxnodes_default, maxnodes + 10),
                                         by = 5))

  # Run the RF algorithm
  RF_compare = train(IQI ~ .,
                     data = RF_Data,
                     method = custom_rf(),
                     metric = "RMSE",
                     tuneGrid = tuneGrid,
                     trControl = trControl_grid,
                     ntree = ntree)
  summary(RF_compare)

  # Save the tune grid plot
  tune_grid = ggplot(RF_compare$results, aes(x = mtry, y = RMSE, color = as.factor(maxnodes))) +
    geom_line() +
    geom_point(shape = 1) +
    labs(color = "maxnodes") +
    theme(legend.position = "bottom", text = element_text(size = 15),
          legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
    xlab("# Prediction features (mtry)") +
    ylab("RMSE")
  tune_grid %>%
    ggexport(filename = paste("figures/Tuning_Rarefy_", rarefaction_rate,"_Sources_", toString(sources),
                              "_TaxaLevel_", taxalevel, ".pdf", sep = ''))

  # Fit the final RF Model with determined parameters
  message(paste('Fit final RF model with; mtry = ', mtry, 'maxnodes = ', maxnodes, ', ntree = ', ntree))
  tuneGrid = expand.grid(.mtry = RF_compare$bestTune$mtry)
  maxnodes = RF_compare$bestTune$maxnodes

  RF_final = train(IQI ~ .,
                   data = RF_Data,
                   method = "rf",
                   metric = "RMSE",
                   tuneGrid = tuneGrid,
                   trControl = trControl_grid,
                   importance = TRUE,
                   maxnodes = maxnodes,
                   ntree = ntree)
  message(paste("RMSE of final model = ", RF_final$results$RMSE))

  # Test a model of reduced features, starting with 2 features upwards, until RMSE becomes worse when adding a new feature
  results = data.frame(c(1, 1, 1), c(2, 2, 2)) # Add some random data to start the while loop
  colnames(results) = c(1, 2)
  rownames(results) = c("mtry", "RMSE", "Rsquared")
  features = 2

  # Add features until the RMSE doesn't change more than 5e-4 or we reached the max. column and the RMSE is still smaller than the full final model
  while(abs(results[2, ncol(results)] - results[2, ncol(results) - 1]) > 0.00005
        | features == ncol(RF_Data)
        & results[2, ncol(results)] < RF_final$results$RMSE){
    important_feat = order(RF_final$finalModel$importance[, 1], decreasing = T)[1:features]
    if(which(colnames(RF_Data) == "IQI") %in% important_feat){break} else{important_feat = append(important_feat, which(colnames(RF_Data) == "IQI"))}

    # Generate reduced final model
    RF_final_reduced = train(IQI ~ .,
                             data = RF_Data[, important_feat],
                             method = "rf",
                             metric = "RMSE",
                             trControl = trControl_grid,
                             importance = TRUE,
                             maxnodes = maxnodes,
                             ntree = ntree)

    # Add a new column with the results
    min = which.min(RF_final_reduced$results$RMSE)
    results[1, features] = RF_final_reduced$results$mtry[min]
    results[2, features] = RF_final_reduced$results$RMSE[min]
    results[3, features] = RF_final_reduced$results$Rsquared[min]
    colnames(results)[features] = features

    # Update the number of features used
    features = features + 1
  }
  results[, 1] = NULL

  # When done, select smallest RMSE and its mtry to run the final reduced model
  features = as.numeric(colnames(results)[which.min(results[2, ])])
  important_feat = order(RF_final$finalModel$importance[, 1], decreasing = TRUE)[1:features]
  if(which(colnames(RF_Data) == "IQI") %in% important_feat){break} else{important_feat = append(important_feat, which(colnames(RF_Data) == "IQI"))}
  tuneGrid = expand.grid(.mtry = results[1, which.min(results[2, ])])

  RF_final_reduced = train(IQI ~ .,
                           data = RF_Data[, important_feat],
                           method = "rf",
                           metric = "RMSE",
                           tuneGrid = tuneGrid,
                           trControl = trControl_grid,
                           importance = TRUE,
                           maxnodes = maxnodes,
                           ntree = ntree)

  message(paste("RMSE of reduced final model = ", RF_final_reduced$results$RMSE))
  write.csv(results, file = paste("figures/Importance_reduced_Rarefy_", rarefaction_rate,"_Sources_", toString(sources),
                                  "_TaxaLevel_", taxalevel, ".csv", sep = ''), row.names = TRUE)

  # Make a variance importance plot
  pdf(paste("figures/Importance_Rarefy_", rarefaction_rate,"_Sources_", toString(sources),
            "_TaxaLevel_", taxalevel, ".pdf", sep = ''), width = 1080/72, height = 12)
  randomForest :: varImpPlot(RF_final_reduced$finalModel, type = 1, main = NULL)
  dev.off()

  # Save the final model as randomForest object
  save(RF_final_reduced, file = paste("figures/randomForest_models/Rarefy_", rarefaction_rate,"_Sources_", toString(sources),
                                      "_TaxaLevel_", taxalevel, ".Rdata", sep = ''))

  message("Function: 'build_rf' end")
  return(list(RF_default = RF_default, RF_final = RF_final, RF_final_reduced = RF_final_reduced))
}

#Function: plot_iqi_prediction
#Takes: test data, rarefy, sources
#Returns: NA
#Job: Generate and save a predicted IQI over actual IQI plot.
plot_iqi_prediction = function(test, model, sources, taxalevel, rarefaction_rate){
  message("Function: 'plot_iqi_prediction' start")

  # Make predictions with the reduced model
  test$pred_IQI = predict(model, test)
  used_mtry = which.min(model$results$RMSE)
  RMSE = model$results$RMSE[used_mtry]
  r = model$results$Rsquared[used_mtry]
  Rsquared = 1 - (sum((test$IQI-test$pred_IQI)^2)/sum((test$IQI-mean(test$IQI))^2))

  # Plot predicted vs actual IQI
  pred_plot = ggplot(test, aes(x = pred_IQI, y = IQI)) +
    geom_point(aes(shape = Source, color = Source), size = 3.5) +
    theme(legend.position = "right", text = element_text(size = 18),
          legend.text = element_text(size = 18), legend.title = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    annotate(geom = "text", label = paste("RMSE = ", round(RMSE, 4),
                                          "\nr2 = ", round(r, 4),
                                          "\nR2 = ", round(Rsquared, 4)),
             x = 0, y = Inf, hjust = 0, vjust = 1, size = 6.5) +
    geom_abline(intercept = 0,
                slope = 1,
                color = "brown3") +
    scale_color_brewer(palette = "Paired", aesthetics = c("fill", "colour")) +
    geom_hline(yintercept = 0.64, col = "green") +
    xlab("Predicted IQI") + xlim(0, 1) +
    ylab("IQI") + ylim(0, 1)

  pred_plot %>%
    ggexport(filename = paste("figures/Prediction_Rarefy_", rarefaction_rate,"_Sources_", toString(sources),
                              "_TaxaLevel_", taxalevel, ".pdf", sep = ''),
             width = 10, height = 8)

  # Plot predicted vs actual IQI per distance
  test_long = pivot_longer(test, cols = ends_with("IQI"))
  dist_plot = ggplot(test_long, aes(x = Distance, y = value, col = name)) +
    geom_line() +
    geom_hline(yintercept = 0.64, col = "green") +
    scale_color_brewer(palette = 'Paired', aesthetics = c("colour", "fill"),
                       guide_legend(title="IQI"),
                       labels = c("Actual", "Predicted")) +
    ylab("IQI") + ylim(0, 1) +
    ggtitle("(a)") +
    theme(plot.title = element_text(size = 12)) +
    coord_cartesian(xlim = c(0, 250))

  # Plot predicted vs actual IQI per month
  month_plot = ggplot(test_long, aes(x = Month, y = value, fill = name)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.64, col = "green") +
    scale_color_brewer(palette = 'Paired', aesthetics = c("colour", "fill"),
                       guide_legend(title="IQI"),
                       labels = c("Actual", "Predicted")) +
    ylab("IQI") + ylim(0, 1) +
    ggtitle("(b)") +
    theme(plot.title = element_text(size = 12)) +
    scale_x_discrete(limits = month.abb)

  ggarrange(dist_plot, month_plot, nrow = 2, ncol = 1, common.legend = T, legend = "bottom") %>%
    ggexport(filename = paste("figures/Prediction_Distance_Month_Rarefy_", rarefaction_rate,"_Sources_", toString(sources),
                              "_TaxaLevel_", taxalevel, ".pdf", sep = ''), width = 6, height = 9)

  message("Function: 'plot_iqi_prediction' end")
}

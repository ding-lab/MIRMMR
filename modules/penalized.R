#Module performs penalized regression

#Set seed for repeatable results (cross-validation approach is stochastic)
if(opt$s != FALSE){
  set.seed(opt$s)
}

suppressPackageStartupMessage(library(glmnet))

#Data to be used
if( opt$train ){
  train_set <- sample(c(TRUE,FALSE), size=nrow(df), replace=TRUE, prob=c(opt$train_proportion, 1-opt$train_proportion))
} else{
  train_set <- rep(TRUE, nrow(df))
}

trainX <- as.matrix(df[train_set,fdc:ncol(df)])
trainY <- as.numeric(df[train_set,col])

if( opt$train ){
  testX <- as.matrix(df[!train_set,fdc:ncol(df)])
  testY <- as.matrix(df[!train_set,fdc:ncol(df)])
}

# If user specified to run consensus step
if( opt$consensus ){
  parameter_counts <- consensus_parameters(opt, trainX, trainY)
}

# Find best model using penalized regression
best_model <- best_lambda_model(opt, trainX, trainY)

# Use test data to evaluate training model
if( opt$train ){
  test_roc <- create_test_roc(best_model, testX, testY)
  test_auc <- auc( test_roc )
  parameter_names <- c("AUC","(Intercept)", names(df[,fdc:ncol(df)]]))
  beta_values <- c(test_auc, rep(0,length(parameter_names)))
  count <- 1
  for(param in parameter_names){
    count <- count + 1
    if( param %in% names(best_model$coefficients) ){
      beta_values[i] <- best_model$coefficients[[param]]
    }
  }
  output_df <- data.frame(value=beta_values)
  write.table(output_df, file=paste0(output_dir_prefix,".penalized_test.txt"), quote=FALSE, sep="\t") 
}

# Save best model
save(best_model, file=paste0(output_dir_prefix,".penalized_model.Robj"))

### PLOT THINGS ###

#Plot predicted model value vs. MSI status
plot_df <- data.frame( trainX, predicted=best_model$fitted.values )
plot_predicted( plot_df )

#Plot consensus model vs. best model 
if(opt$consensus){
  plot_df <- data.frame(parameter_counts, in_best_model=(parameter_counts$Parameter %in% names(best_model$coefficients)))
  plot_consensus( plot_df )
}

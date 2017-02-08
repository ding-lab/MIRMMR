#Module performs penalized regression

#Set seed for repeatable results (cross-validation approach is stochastic)
if(opt$set_seed != FALSE){
  set.seed(opt$set_seed)
}

suppressPackageStartupMessages(library(glmnet))

if( opt$parallel ){
  suppressPackageStartupMessages(library(doMC))
  registerDoMC(cores=opt$par_cores)
}

#Data to be used
if( opt$train ){
  train_set <- sample(c(TRUE,FALSE), size=nrow(df), replace=TRUE, prob=c(opt$train_proportion, 1-opt$train_proportion))
} else{
  train_set <- rep(TRUE, nrow(df))
}

trainX <- as.matrix(df[train_set,fdc:ncol(df)])
trainY <- as.matrix(df[train_set,col])

if( opt$train ){
  testX <- as.matrix(df[!train_set,fdc:ncol(df)])
  testY <- as.matrix(df[!train_set,col])
}

# If user specified to run consensus step
if( opt$consensus ){
  print('gets here')
  parameter_counts <- consensus_parameters(opt, trainX, trainY)
  print('gets here 2')
}

# Find best model using penalized regression
print('gets here 3')
best_model <- best_lambda_model(opt, trainX, trainY)
print('gets here 4')

# Use test data to evaluate training model
if( opt$train ){
  test_roc <- create_test_roc(best_model, as.data.frame(testX), testY)
  test_auc <- auc( test_roc )
  parameter_names <- c("AUC","(Intercept)", names(df[,fdc:ncol(df)]))
  beta_values <- rep(0,length(parameter_names))
  beta_values[1] <- test_auc
  count <- 1
  for(param in parameter_names[-1]){
    count <- count + 1
    if( param %in% c("(Intercept)",names(best_model$coefficients)) ){
      beta_values[count] <- best_model$coefficients[[param]]
    }
  }
  output_df <- data.frame(parameter=parameter_names,value=beta_values)
  write.table(output_df, file=paste0(output_dir_prefix,".penalized_test.txt"), quote=FALSE, sep="\t", row.names=FALSE) 
}

# Save best model
save(best_model, file=paste0(output_dir_prefix,".penalized_model.Robj"))

### PLOT THINGS ###
if( opt$plots ){
  suppressPackageStartupMessages(library(ggplot2))
  #Plot predicted model value vs. MSI status
  plot_df <- data.frame( status=df[train_set,col], predicted=best_model$fitted.values, group=df[train_set,opt$group] )
  plot_predicted( plot_df, xlab=opt$xlab, ylab=opt$ylab, color_indicates=opt$color_indicates, theme=opt$theme_bw)

  #Plot consensus model vs. best model 
  if(opt$consensus){
    plot_df <- data.frame(parameter_counts, in_best_model=(parameter_counts$Parameter %in% names(best_model$coefficients)))
    plot_consensus( plot_df, xlab="Number of models", ylab="Model variable", color_indicates="Included in\n'best' model", theme=opt$theme_bw )
  }
}

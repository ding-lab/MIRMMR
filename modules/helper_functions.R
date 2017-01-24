#Set of functions to help facilitate rest of scripts

#Function +++++++++++++++++++++++++++++
#Track how many times each variable appears in models
consensus_parameters <- function(opt, myX, myY){
  get_cvglmnet_variables <- function(opt, myX, myY){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$nfolds, family="binomial", type.measure=opt$type_measure, alpha=opt$alpha, parallel=opt$parallel)
    return(row.names(coef(cvfit, s=opt$lambda))[coef(cvfit, s=opt$lambda)[,1] != 0])
  }

  variables_in_model <- replicate(opt$repeats, get_cvglmnet_variables(opt, myX, myY))
  
  variables_in_model_df <- data.frame(table(unlist(variables_in_model)))
  names(variables_in_model_df) <- c("Parameter","Count")
  return(variables_in_model_df)
}
#--------------------------------------

#Find best lambda and return model with best lambda coefficients
#http://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results/173895#173895
best_lambda_model <- function(opt, myX, myY){
  get_lambdas_errors <- function(opt, myX, myY){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$nfolds, family="binomial", type.measure=opt$type_measure, alpha=opt$alpha, parallel=opt$parallel)
    return(data.frame(cvfit$lambda,cvfit$cvm))
  }
  
  lambdas <- replicate(opt$repeats, get_lambdas_errors(opt, myX, myY))
  
  # take mean cvm for each lambda
  lambdas <- aggregate(unlist(lambdas[2,]), list(unlist(lambdas[1,])), mean)

  # select the best lambda
  bestindex <- which(lambdas[,2]==min(lambdas[,2]))
  bestlambda <- min(lambdas[bestindex,1])

  # and now run glmnet once more with the best lambda
  glmnet.fit <- glmnet(x=myX, y=myY, lambda=bestlambda, family="binomial", alpha=opt$alpha)

  # Determine the coeffients in model derived from best lambda
  best_coef_indicator <- coef(glmnet.fit)[,1] != 0
  coefs_in_best_model <- rownames(coef(glmnet.fit))[best_coef_indicator]

  # Fit a glm model using coefficients from above
  model_df <- data.frame(myY, myX[, best_coef_indicator[-1]])
  names(model_df)[1] <- "outcome"
  model <- glm( outcome ~ ., data=model_df, family="binomial")
  return(model)
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Calculate points in a ROC curve
roc <- function( input ){
  #Input data frame has two columns: truth, score
  #Truth is what we assume to be true (like TRUE, FALSE)
  #Score is a numeric value on a scale (like a probability between 0,1)
  #Output is a data frame starting from 0,0 up to 1,1 of values in ROC curve
  output <- matrix(NA, nrow(input)+1, 2)
  output[1,] <- c(0,0)
  scores <- sort(input[,2], decreasing=TRUE)
  for(i in 1:nrow(input)){
    threshold <- scores[i]
    num_tp <- sum(input[,2] >= threshold & input[,1])
    num_con_true <- sum(input[,1])
    num_fp <- sum(input[,2] >= threshold & !input[,1])
    num_con_false <- sum(!input[,1])
    tp <- num_tp/num_con_true
    fp <- num_fp/num_con_false
    output[i+1,] <- c(tp, fp)
  }
  return(output)
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Calculate area under a ROC curve (AUC)
auc <- function( roc ){
  #Input ROC has two columns: tp, fp
  cum_auc <- 0
  for(i in 1:(nrow(roc)-1) ){
    rect <- roc[i,1]*(roc[i+1,2]-roc[i,2]) 
    tri <- (1/2)*(roc[i+1,2]-roc[i,2])*(roc[i+1,1]-roc[i,1])
    cum_auc <- cum_auc + rect + tri
  }
  return(cum_auc)
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Test accuracy of model using test data
create_test_roc <- function( model, testX, testY ){
  test_predictions <- predict( model, newdata=testX, type="response")
  test_roc <- data.frame(truth=testY, score=test_predictions)  
  return( roc( test_roc ) )
}
#--------------------------------------

### PLOTS ###
#Function +++++++++++++++++++++++++++++
#Plot a ROC curve
plot_roc <- function( plot_df, user_title, module){
  #Function plots a ROC curve given (possible) multiple curves
  #Input is a three column data frame (method, fpr, tpr)
  p <- ggplot(plot_df, aes(x=fpr, y=tpr, color=method))
  p <- p + geom_line(size=1)
  p <- p + geom_abline(intercept=0, slope=1, linetype=2, alpha=0.5)
  p <- p + labs(x="False positive rate", y="True positive rate", title=user_title, color="Method")
  pdf(paste0(output_dir_prefix,".",module,"_roc_curve.pdf"),10,10)
  print(p)
  dev.off()
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Plot predicted model value vs. MSI status
plot_predicted <- function( plot_df ){
  p <- ggplot(plot_df, aes(x=status, y=predicted))
  p <- p + geom_violin()
  p <- p + geom_jitter(aes(color=group), width=0.3, height=0)
  p <- p + labs(x="MSI status", y="Fitted probability of MSI status", title="Penalized regression probability of MSI status", color=opt$group)
  p <- p + scale_y_continuous(limits = c(0,1))
  pdf(paste0(output_dir_prefix,".penalized_predicted.pdf"),10,10)
  print(p)
  dev.off()
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Plot consensus model vs. best model
plot_consensus <- function( plot_df ){
  p <- ggplot(plot_df, aes(x=Count, y=reorder(Parameter,Count), color=in_best_model))
  p <- p + geom_point(size=3)
  p <- p + labs(x="Number of models", y="Model variable", title="Number of models in which each variable appears", color="Included in\n'best' model")
  pdf(paste0(output_dir_prefix,".penalized_consensus.pdf"),10,10)
  print(p)
  dev.off()
}
#--------------------------------------

### SANITY CHECKS ###

#Function +++++++++++++++++++++++++++++
#Check inputs for sanity and correctness, exit if any problem
sanity_checks <- function(opt){
  em <- function(old_message, new_message){
    return( paste0(old_message, new_message, sep="\n") )
  }
  error_message <- NULL
  
  #Primary checks of module and data_frame
  #module
  list_of_modules <- c("compare", "penalized", "predict", "stepwise", "univariate")
  if( !(opt$module %in% list_of_modules) ){
    error_message <- em(error_message, paste0("Error: module (-m) must be one of ", paste(list_of_modules, collapse=", "),"."))
  }
  if( !is.null(error_message) ){
    stop(error_message)
  }
  #data_frame
  if( is.null(opt$data_frame) ){
    error_message <- em(error_message, "Error: data frame (-f) is NULL and must be specified.")
  } else if( !(file.exists(opt$data_frame)) ){
    error_message <- em(error_message, "Error: data frame (-f) file does not exist.")
  } else if( file.access(opt$data_frame, mode=4)!=0 ){
    error_message <- em(error_message, "Error: data frame (-f) file is not readable.")
  } else if( TRUE ){
    default_warn <- options()$warn
    options(warn=2)
    df <- read.table(opt$data_frame, header=TRUE)
    options(warn=default_warn)
  }
  if( !is.null(error_message) ){
    stop(error_message)
  }

  #Secondary checks of required options (msi_status, first_data_column, output_prefix, output_directory, overwrite, plots, model)
  #msi_status
  if( !(opt$module %in% c("predict","compare") )){
    if( is.null(opt$msi_status) ){
      error_message <- em(error_message, "Error: MSI status column name (-i) must be specified.")
    } else if( !(opt$msi_status %in% names(df))  ){
      error_message <- em(error_message, "Error: MSI status column name (-i) not in data frame column names.")
    } else if( length(levels(droplevels(as.factor(df[,opt$msi_status]))))!=2 ){
      error_message <- em(error_message, "Error: MSI status column (-i) may only have two levels.")
    } 
    if( !is.null(error_message) ){
      stop(error_message)
    }  
  }
  #first_data_column
  if( is.null(opt$first_data_column) ){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be specified.")
  } else if( !is.numeric(opt$first_data_column) || opt$first_data_column!=round(opt$first_data_column) ){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be an integer.")
  } else if( opt$first_data_column >= ncol(df) ){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be less than the number of columns of the input data frame.")
  } else if( !(opt$module %in% c("predict","compare")) && opt$first_data_column <= which(names(df)==opt$msi_status )){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be greater than the column number of the MSI status column.")
  }
  #output_prefix
  if( is.null(opt$output_prefix) ){
    error_message <- em(error_message, "Error: Output file prefix (-o) must be specified.")
  }
  #output_directory
  if( is.null(opt$output_directory) ){
    error_message <- em(error_message, "Error: Output file directory (-d) must be specified.")
  }
  #overwrite
  if( !is.logical(opt$overwrite) ){
    error_message <- em(error_message, "Error: Overwrite (--overwrite) must be logical TRUE/FALSE, default is FALSE.")
  }
  #plots
  if( !is.logical(opt$plots) ){
    error_message <- em(error_message, "Error: Plots (--plots) must be logical TRUE/FALSE, default is FALSE.")
  }
  #model
  if( opt$module=="predict" ){
    if( is.null(opt$model) ){
      error_message <- em(error_message, "Error: model (--model) is NULL and must be specified.")
    } else if( !(file.exists(opt$model)) ){
      error_message <- em(error_message, "Error: model (--model) file does not exist.")
    } else if( file.access(opt$model, mode=4)!=0 ){
      error_message <- em(error_message, "Error: model (--model) file is not readable.")
    }
  }
  #Exit if error_message not NULL
  if( !is.null(error_message) ){
    stop(error_message)
  }

  #Tertiary chcks of penalized module options (alpha, consensus, lambda, nfolds, repeats, set_seed, train, train_proportion, type_measure)
  if( opt$module == "penalized" ){
    #plots
    if( !is.logical(opt$plots) ){
      error_message <- em(error_message, "Error: Plots (--plots) must be logical TRUE/FALSE, default is FALSE.")
    }
    #group
    if( !is.null(opt$group) && !(opt$group %in% names(df)) ){
      error_message <- em(error_message, "Error: Group (--group) is not found within the column names of the input data frame.")
    }
    #alpha
    if( is.na(opt$alpha) || opt$alpha < 0 | opt$alpha > 1 ){
      error_message <- em(error_message, "Error: Alpha (--alpha) must be a double beteween 0 and 1, default=0.9.")
    }
    #consensus
    if( is.na(opt$consensus) ){
      error_message <- em(error_message, "Error: Consensus (--consensus) must be logical TRUE/FALSE, default is FALSE.")
    }
    #lambda
    if( !(opt$lambda %in% c("lambda.1se","lambda.min")) ){
      error_message <- em(error_message, "Error: Lambda (--lambda) must be one of \"lambda.1se\" or \"lambda.min\", default is \"lambda.1se\".")
    }
    #nfolds
    if( is.na(opt$nfolds) || opt$nfolds!=round(opt$nfolds) | opt$nfolds < 1 | opt$nfolds > nrow(df) ){
      error_message <- em(error_message, "Error: Number of folds (--nfolds) must be an integer between 1 and the number of samples in the input data frame, default is 10.")
    } 
    #repeats
    if( is.na(opt$repeats) || opt$repeats!=round(opt$repeats) | opt$repeats < 1){
      error_message <- em(error_message, "Error: Number of repeats (--repeats) must be a positive integer, default is 1000.")
      error_message <- em(error_message, "Error: Parallel (--parallel) must be logical TRUE/FALSE, default is FALSE.")
    }
    #par_cores
    if( is.na(opt$par_cores) || opt$par_cores!=round(opt$par_cores) | opt$par_cores < 1 ){
      error_message <- em(error_message, "Error: Parallel cores (--par_cores) must be an integer greater than 0.")
    }
    #set_seed
    if( is.na(opt$set_seed) ){
      error_message <- em(error_message, "Error: Random seed (--seed) must be a number, default is 0 (not set).")
    }
    #train
    if( is.na(opt$train) ){
      error_message <- em(error_message, "Error: Train option (--train) must be logical TRUE or FALSE, default is FALSE.")
    }
    #train_proportion
    if( is.na(opt$train_proportion) || opt$train_proportion < 0 | opt$train_proportion > 1 ){
      error_message <- em(error_message, "Error: Train proportion (--train_proportion) must be a double between 0 and 1, default is 0.8.")
    }
    #type_measure
    if( !(opt$type_measure %in% c("mse", "deviance", "mae", "class", "auc")) ){
      error_message <- em(error_message, "Error: Type measure (--type_measure) must be one of \"mse\", \"deviance\", \"mae\", \"class\", \"auc\", default is \"class\".")
    }
  }
  #Exit if error_message not NULL
  if( !is.null(error_message) ){
    stop(error_message)
  }
}
#--------------------------------------
#Set of functions to help facilitate rest of scripts

#Function +++++++++++++++++++++++++++++

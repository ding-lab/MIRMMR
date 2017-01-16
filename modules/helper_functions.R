#Set of functions to help facilitate rest of scripts

#Function +++++++++++++++++++++++++++++
#Track how many times each variable appears in models
consensus_parameters <- function(opt, myX, myY){
  #TODO update opt parameter names
  variables_in_model <- vector("list", opt$n)
  for(i in 1:opt$n){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
    variables_in_model[[i]] <- row.names(coef(cvfit, s=opt$l))[coef(cvfit, s=opt$l)[,1] != 0]
  }
  variables_in_model_df <- data.frame(table(unlist(variables_in_model)))
  names(variables_in_model_df) <- c("Parameter","Count")
  return(variables_in_model_df)
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Find best lambda and return model with best lambda coefficients 
#Modified from stackexchange user 'Sideshow Bob', answer from 2016/03/31. Thanks!
#http://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results/173895#173895
best_lambda_model <- function(opt, myX, myY){
  #TODO update opt parameter names
  lambdas = NULL
  for (i in 1:opt$n){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
    #TODO: optimize this with prespecified data frame or list structure for lambdas
    errors <- data.frame(cvfit$lambda,cvfit$cvm)
    lambdas <- rbind(lambdas,errors)
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[,2], list(lambdas$cvfit.lambda), mean)

  # select the best lambda
  bestindex = which(lambdas[,2]==min(lambdas[,2]))
  bestlambda = lambdas[bestindex,1]

  # and now run glmnet once more with the best lambda
  glmnet.fit <- glmnet(x=myX, y=myY, lambda=bestlambda, family="binomial", alpha=opt$a)

  # Determine the coeffients in model derived from best lambda
  best_coef_indicator <- coef(glmnet.fit)[,1] != 0
  coefs_in_best_model <- rownames(coef(glmnet.fit))[best_coef_indicator]
  
  # Fit a glm model using coefficients from above
  model_df <- data.frame(myY, myX[, names(myX) %in% coefs_in_best_model])
  names(model_df)[1] <- "outcome"
  model <- glm( outcome ~ ., data=model_df, family=binomial)
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
  output <- matrix(NA, ncol(input)+1, 2)
  output[i,] <- c(0,0)
  scores <- sort(input[,2])
  for(i in 1:ncol(input)){
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
  for(i in 1:ncol(roc)-1 ){
    rect <- (roc[i+1,1]-roc[i,1])*(roc[i,2])*(roc[i,2]) 
    tri <- (1/2)*(roc[i+1,1]-roc[i,1])*(roc[i+1,2]-roc[i,2])
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

#Function +++++++++++++++++++++++++++++
#Plot a ROC curve
plot_roc <- function( plot_df ){
  #Function plots a ROC curve given (possible) multiple curves
  #Input is a three column data frame (method, fpr, tpr)
  p <- ggplot(plot_df, aes(x=fpr, y=tpr, color=method))
  p <- p + geom_line(size=1)
  p <- p + geom_abline(intercept=0, slope=1, linetype=2, alpha=0.5)
  p <- p + labs(x="False positive rate", y="True positive rate", title="ROC curve", color="Method")
  pdf()
  print(p)
  dev.off()
}
#--------------------------------------

### PLOTS ###

#Function +++++++++++++++++++++++++++++
#Plot predicted model value vs. MSI status
plot_predicted <- function( plot_df ){
  p <- ggplot(plot_df, aes())
  p <- p + geom_violin()
  p <- p + geom_jitter()
  p <- p + labs(x="", y="", title="", color="")
  pdf()
  print(p)
  dev.off()
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Plot consensus model vs. best model
plot_consensus <- function(){
  p <- ggplot(plot_df, aes())
  p <- p + geom_point()
  p <- p + labs(x="", y="", title="", color="")
  pdf()
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
  list_of_modules <- c("compare", "penalized", "stepwise", "univariate")
  if( !(opt$module %in% list_of_modules) ){
    error_message <- em(error_message, paste0("Error: module (-m) must be one of ", paste(list_of_modules, collapse=" "),"."))
  }
  if( !is.null(error_message) ){
    stop(error_message)
  }
  #data_frame
  if( is.null(opt$data_frame) ){
    error_message <- em(error_message, "Error: data frame (-d) is NULL and must be specified.")
  } else if( !(file.exists(opt$data_frame)) ){
    error_message <- em(error_message, "Error: data frame (-d) file does not exist.")
  } else if( !(file.access(opt$data_frame, mode=4)) ){
    errro_message <- em(error_message, "Error: data frame (-d) file is not readable.")
  } else if( TRUE ){
    default_warn <- options()$warn
    options(warn=2)
    df <- read.table(opt$data_frame, header=TRUE)
    options(warn=default_warn)
  }
  if( !is.null(error_message) ){
    stop(error_message)
  }

  #Secondary checks of required options (msi_status, first_data_column, output_prefix, output_directory, overwrite, plots)
  #msi_status
  if( is.null(opt$msi_status) ){
    error_message <- em(error_message, "Error: MSI status column name (-i) must be specified.")
  } else if( !(opt$msi_status %in% names(df))  ){
    error_message <- em(error_message, "Error: MSI status column name (-i) not in data frame column names.")
  } else if( length(levels(droplevels(as.factor(df[,opt$msi_status]))))!=2 ){
    error_message <- em(error_message, "Error: MSI status column (-i) may only have two levels.")
  } 
  #first_data_column
  if( is.null(opt$first_data_column) ){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be specified.")
  } else if( !is.double(opt$c) || opt$c!=round(opt$c) ){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be an integer.")
  } else if( opt$c >= ncol(df) | opt$c <= which(names(df)==opt$msi_status ){
    error_message <- em(error_message, "Error: Number of first data column (-c) must be greater than the column number of the MSI status column and must be less than number of columns of the input data frame.")
  }
  #output_prefix
  if( is.null(opt$output_prefix) ){
    error_message <- em(error_message, "Error: Output file prefix (-o) must be specified.")
  }
  #output_directory
  if( is.null(opt$output_directory) ){
    error_message <- em(error_message, "Error: Output file directory (-d) must be specified.")
  }
  #TODO check if file exists and will be overwritten if overwrite=FALSE
  #plots
  if( !is.logical(opt$plots) ){
    error_message <- em(error_message, "Error: Plots (--plots) must be logical TRUE/FALSE, default is FALSE.")
  }
  #Exit if error_message not NULL
  if( !is.null(error_message) ){
    stop(error_message)
  }

  #Tertiary chcks of penalized module options (alpha, consensus, lambda, nfolds, repeats, set_seed, train, train_proportion, type_measure)
  if( opt$module == "penalized" ){
    #alpha
    if( !is.double(opt$alpha) || opt$alpha < 0 | opt$alpha > 1 ){
      error_message <- em(error_message, "Error: Alpha (--alpha) must be a double beteween 0 and 1, default=0.9.")
    }
    #consensus
    if( !is.logical(opt$consensus) ){
      error_message <- em(error_message, "Error: Consensus (--consensus) must be logical TRUE/FALSE, default is FALSE.")
    }
    #lambda
    if( !(opt$lambda %in% c("lambda.1se","lambda.min")) ){
      error_message <- em(error_message, "Error: Lambda (--lambda) must be one of \"lambda.1se\" or \"lambda.min\", default is \"lambda.1se\".")
    }
    #nfolds
    if( !is.double(opt$nfolds) || opt$nfolds!=round(opt$nfolds) | opt$nfolds < 1 | opt$nfolds > nrow(df) ){
      error_message <- em(error_message, "Error: Number of folds (--nfolds) must be an integer between 1 and the number of samples in the input data frame, default is 10.")
    } 
    #repeats
    if( !is.double(opt$repeats) || opt$repeats!=round(opt$repeats) | opt$repeats < 1){
      error_message <- em(error_message, "Error: Number of repeats (--repeats) must be a positive integer, default is 1000.")
    }
    #set_seed
    if( !(is.logical(opt$set_seed) | is.double(opt$set_seed)) ){
      error_message <- em(error_message, "Error: Random seed (--seed) must be FALSE, TRUE, or any double, default is FALSE.")
    }
    #train
    if( !is.logical(opt$train) ){
      error_message <- em(error_message, "Error: Train option (--train) must be logical TRUE or FALSE, default is FALSE.")
    }
    #train_proportion
    if( !is.double(opt$train_proportion) || opt$train_proportion < 0 | opt$train_proportion > 1 ){
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

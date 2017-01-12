#Script checks command line inputs for sanity and returns an error message if things are not correct

#primary checks of command line arguments
error_message <- NULL
#check module
if( is.null(opt$module) | !(opt$module %in% c("compare", "penalized", "stepwise", "univariate")) ){
  error_message <- paste0(error_message, "Module (-m) must be one of compare, penalized, stepwise, univariate.\n")
} 
#check data_frame
if( is.null(opt$data_frame) | !file.exists(opt$data_frame) ){
  error_message <- paste0(error_message, "Data frame (-d) must be specified, or file does not exist.\n")
}

#proceed with importing data_frame
if( is.null(error_message) ){
  df <- read.table(opt$data_frame, header=TRUE)
  if( any(is.na(df)) ){
    error_message <- paste0(error_message, "Data frame (-d) contains missing data. User should make appropriate adjustments by removing samples or parameters with missing data or using imputation.")
    stop(error_message)
  } 
} else{
  stop(error_message)
}

#secondary checks of command line arguments
#check msi_score
if( is.null(opt$msi_score) ){
  error_message <- paste0(error_message, "MSI score (-i) column name or number must be specified.\n")
} else if( !(opt$msi_score %in% names(df)) | (is.numeric(opt$msi_score) & opt$msi_score > ncol(df)-1) ){
  error_message <- paste0(error_message, "MSI score (-i) column name does not exist in data frame or column number is out of range.\n")
} else if( length(levels(droplevels(as.factor(df[,opt$msi_score])))) < 2 ){
  error_message <- paste0(error_message, "MSI score (-i) column contains fewer than two levels (must contain two factors or have many numeric values).\n")
}

#check first_data_column
if( is.null(opt$c) | opt$c < 2 | opt$c > ncol(df) | opt$c!=round(opt$c)){
  error_message <- paste0(error_message, "First data column (-c) must be specified, or number is not an integer or out of bounds.\n")
}
#check verbose
if( !(opt$verbose==TRUE | opt$verbose==FALSE) ){
  error_message <- paste0(error_message, "Verbose (-v) must be either TRUE or FALSE, default is TRUE.\n") 
}

#proceed to next checks
if( !is.null(error_message) ){
  stop(error_message)
}

#tertiary check of command line arguments
if( opt$module == "penalized" ){
  #check alpha
  if( !is.numeric(opt$alpha) | opt$alpha < 0 | opt$alpha > 1 ){
    error_message <- paste0(error_message, "Penalized regression parameter alpha (-a) must be number between 0 and 1, default is 0.9.\n")
  }
  #check nfolds
  if( !is.numeric(opt$nfolds) | opt$nfolds!=round(opt$nfolds) | opt$nfolds < 0 ){
    error_message <- paste0(error_message, "Penalized regression parameter nfolds (-f) must be an integer greater than 0, default is 10.\n")
  }
  #check lambda
  if( !(opt$lambda %in% c("lambda.1se", "lambda.min")) ){
    error_message <- paste0(error_message, "Penalized regression parameter lambda (-l) must be one of lambda.1se or lambda.min, default is lambda.1se.\n")
  }
  #check number_repetitions
  if( !is.numeric(opt$n) | opt$n!=round(opt$n) | opt$n<0 ){
    error_messsage <- paste0(error_message, "Penalized regression parameter number_repetitions (-n) must be a positive integer.\n")
  }
  #check train_proportion
  if( !is.numeric(opt$p) | opt$p < 0.5 | opt$p > 1){
    error_message <- paste0(error_message, "Penalized regression parameter train_proportion (-p) must be a number between 0.5 and 1.\n")
  }
  #check set_seed
  if( !is.logical(opt$s) | !is.numeric(opt$s) ){
    error_message <- paste0(error_message, "Penalized regression parameter set_seed (-s) must either be TRUE, FALSE, or a number, default is FALSE.\n")
  }
  #check type_measure
  if( !(opt$type_measure %in% c("mse", "deviance", "mae", "class", "auc")) ){
    error_message <- paste0(error_message, "Penalized regression parameter type_measuge (-t) must be mse, deviance, mae, class, or auc, default is class.\n")
  }
  #TODO: IMPLEMENT MORE PARAMETER CHECKS
}

#First check if MSI score column is binary (logistic) or not (linear)
col <- opt$msi_score
fdc <- opt$first_data_column

if( length(levels(droplevels(as.factor(df[,col]))))==2 ){
  msi_binary <- TRUE
} else if( length(levels(droplevels(as.factor(df[,col]))))>2 ){
  msi_binary <- FALSE
}
else{
  error_message <- paste0(error_message, "MSI score column should have more than one level.\n")
}

#This script parses command line arguments and runs specified analysis module

suppressPackageStartupMessages(library("optparse"))

#set options
option_list <- list(
  make_option(c("-m", "--module"), default=NULL, type="character", help="What module should be run (compare, penalized, stepwise, univariate), must be specified"),
  make_option(c("-d", "--data_frame"), default=NULL, type="character", help="File path to data frame, must be specified"),
  make_option(c("-i", "--msi_score"), default=NULL, type="character", help="Column name referring to MSI score (stored as numeric or binary vector), must be specified"),
  make_option(c("-c", "--first_data_column"), default=NULL, type="integer", help="Column number referring to first parameter in data frame, assuming all higher columns are also parameters, must be specified"),
  make_option(c("-a", "--alpha"), default=0.9, type="double", help="Parameter alpha used in glmnet::glmnet, default=%default"),
  make_option(c("-f", "--nfolds"), default=10, type="integer", help="Parameter nfolds used in glmnet::cv.glmnet, default=%default"),
  make_option(c("-l", "--lambda"), default="lambda.1se", type="character", help="Parameter lambda used in glmnet::cv.glmnet, default=\"%default\", options: lambda.1se, lambda.min"),
  make_option(c("-n", "--number_repetitions"), default=1000, type="integer", help="Number of times to repeat testing to determine optimal lambda in penalized module, default=%default"),
  make_option(c("-p", "--train_proportion"), default=0.8, type="double", help="Proportion of samples to retain in training set in penalized module, default=%default"),
  make_option(c("-s", "--set_seed"), default=FALSE, type=NULL, help="Option to set seed before penalized module, seed can be set to any number, default=\"%default\""),
  make_option(c("-t", "--type_measure"), default="class", type="character", help="Parameter type.measure used in glmnet::cv.glmnet, default=\"%default\", options: mse, deviance, mae, class, auc"),
  make_option(c("-v", "--verbose"), default=TRUE, type="logical", help="Print extra output along the way, default=%default")
  make_option(c("--test"), default=FALSE, type="logical", help="Internal use: testing penalized module results."),
  make_option(c("--time"), default=FALSE, type="logical", help="Internal use: report time to complete analysis.")
)

#retrieve command line arguments
opt <- parse_args(OptionParser(usage="%prog -m MODULE -d DATA_FRAME -i MSI_SCORE -c FIRST_DATA_COLUMN [options]\n\nMSImodel assumes the data frame has meta information columns (e.g. sample name, cancer type, MSI score), followed by data columns (i.e. parameters that will be analyzed) and also assumes header=TRUE", option_list=option_list))

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
if( is.null(opt$msi_score) | !(opt$msi_score %in% names(df)) ){
  error_message <- paste0(error_message, "MSI score (-i) column name must be specified, or column name does not exist in data frame.\n")
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
}

#If there is no error message to report, run the specified module
if( !is.null(error_message) ){
  stop(error_message)
} else{
  #Run the specified module
  #Need this is refer to script directory MSImodel/modules
  #http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
  source(paste0("modules/", opt$module, ".R"))
}

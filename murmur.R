#This script parses command line arguments and runs specified analysis module

suppressPackageStartupMessages(library("optparse"))

#Set options
option_list <- list(
  make_option(c("-m", "--module"), default=NULL, type="character", help="What module should be run (compare, penalized, stepwise, univariate), must be specified"),
  make_option(c("-f", "--data_frame"), default=NULL, type="character", help="File path to data frame, must be specified"),
  make_option(c("-i", "--msi_status"), default=NULL, type="character", help="Column name referring to MSI status (stored as binary vector), must be specified"),
  make_option(c("-c", "--first_data_column"), default=NULL, type="double", help="Column number referring to first parameter in data frame, assuming all higher columns are also parameters, must be specified"),
  make_option(c("-o", "--output_prefix"), default=NULL, type="character", help="Output file name prefix, must be specified"),
  make_option(c("-d", "--output_directory"), default=NULL, type="character", help="Output directory, must be specified"),
  make_option(c("--overwrite"), default=FALSE, type="logical", help="Prevent overwriting existing files, turn off overwrite protection with TRUE, default=%default."),
  make_option(c("--plots"), default=TRUE, type="logical", help="Compare or Penalized module: Produce informative plots in the compare or penalized module, default=%default"),
  make_option(c("--group"), default=NULL, type="character", help="Penalized module: Column name referring to a group identifier (e.g. cancer type) used in plotting, default=%default"),
  make_option(c("--alpha"), default=0.9, type="double", help="Penalized module: Parameter alpha used in glmnet::glmnet, default=%default"),
  make_option(c("--consensus"), default=FALSE, type="logical", help="Penalized module: Use consensus method in additional to best lambda method to find set of coefficients that appear in most models, only useful when used with --plots, default=%default"),
  make_option(c("--lambda"), default="lambda.min", type="character", help="Penalized module: Parameter lambda used in glmnet::cv.glmnet, default=\"%default\", options: lambda.1se, lambda.min"),
  make_option(c("--nfolds"), default=10, type="double", help="Penalized module: Parameter nfolds used in glmnet::cv.glmnet, default=%default"),
  make_option(c("--parallel"), default=FALSE, type="logical", help="Penalized module: Option to use multiple cores when fitting different folds in glmnet::cv.glmnet, default=%default"),
  make_option(c("--par_cores"), default=1, type="double", help="Penalized module: Number of cores to employ when using parallel (--parallel) option, should not exceed the number of folds (--nfolds), default=%default"),
  make_option(c("--repeats"), default=1000, type="double", help="Penalized module: Number of times to repeat testing to determine optimal lambda in penalized module, default=%default"),
  make_option(c("--set_seed"), default=0, type="double", help="Penalized module: Option to set seed before penalized module, seed can be set to any number except 0, default=%default"),
  make_option(c("--train"), default=FALSE, type="logical", help="Penalized module: Option to use a training set/test set approach to measure accuracy of best lambda approach, default=%default"),
  make_option(c("--train_proportion"), default=0.8, type="double", help="Penalized module: With --train=TRUE, the proportion of samples to keep in the training set, default=%default"),
  make_option(c("--type_measure"), default="class", type="character", help="Penalized module: Parameter type.measure used in glmnet::cv.glmnet (options: mse, deviance, mae, class, auc), default=%default"),
  make_option(c("--model"), default=NULL, type="character", help="Predict module: Path to model object (e.g. .Robj) to be used to predict status of new data, must be specified if using predict module.")
)

#Retrieve command line arguments
opt <- parse_args(OptionParser(usage="%prog -m MODULE -f DATA_FRAME -i MSI_STATUS -c FIRST_DATA_COLUMN -o OUTPUT_PREFIX -d OUTPUT_DIRECTORY [options]\n\nAssumptions:\n1. The data frame has meta information columns (e.g. sample name, cancer type, MSI score)\n2. Followed by data columns (i.e. predictors in regression models)\n3. And that the input data frame has column headers.", option_list=option_list))

#Get path of this script to run sources
#StackOverflow: http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script/1815743#1815743 
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#Source these helper functions
source(paste0(script.basename,"/modules/helper_functions.R"))

#Sanity check on the parameter inputs
sanity_checks(opt)

#df is the name of the input data frame
df <- read.table(opt$data_frame, header=TRUE)
#col is the column name referring to MSI status
col <- opt$msi_status
#fdc is the column number of the first data column
fdc <- opt$first_data_column
#output directory and file prefix
output_directory <- paste(strsplit(gsub("/+","/",opt$output_directory),"/")[[1]],collapse="/")
file_prefix <- opt$output_prefix
output_dir_prefix <- paste0(output_directory,"/",file_prefix)

#For the specified model, see if the output files already exist and exit if --overwrite=FALSE
if( !opt$overwrite ){
  overwrite_message <- NULL
  if( opt$module == "univariate" ){
    if( file.exists(paste0(output_dir_prefix,".univariate_summary.txt")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".univariate_summary.txt"," file already exists, set --overwrite=TRUE to overwrite.")
    }
    if( file.exists(paste0(output_dir_prefix,".univariate_models.Robj")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".univariate_models.Robj"," file already exists, set --overwrite=TRUE to overwrite.")
    }
  } else if( opt$module == "stepwise" ){
    if( file.exists(paste0(output_dir_prefix,".stepwise_model_summary.txt")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".stepwise_model_summary.txt"," file already exists, set --overwrite=TRUE to overwrite.")      
    }
    if( file.exists(paste0(output_dir_prefix,".stepwise_model.Robj")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".stepwise_model.Robj"," file already exists, set --overwrite=TRUE to overwrite.")
    }
  } else if( opt$module == "penalized" ){
    if( opt$train & file.exists(paste0(output_dir_prefix,".penalized_test.txt")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".penalized_test.txt"," file already exists, set --overwrite=TRUE to overwrite.")
    }
    if( opt$consensus & file.exists(paste0(output_dir_prefix,".penalized_consensus.pdf")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".penalized_consensus.pdf"," file already exists, set --overwrite=TRUE to overwrite.")
    }
    if( file.exists(paste0(output_dir_prefix,".penalized_predicted.pdf")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".penalized_predicted.pdf"," file already exists, set --overwrite=TRUE to overwrite.")
    }
    if( file.exists(paste0(output_dir_prefix,".penalized_model.Robj")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".penalized_model.Robj"," file already exists, set --overwrite=TRUE to overwrite.")
    }
  } else if( opt$module == "compare" ){
    if( file.exists(paste0(output_dir_prefix,".compare_models_roc.pdf")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".compare_models_roc.pdf"," file already exists, set --overwrite=TRUE to overwrite.")
    }
    if( file.exists(paste0(output_dir_prefix,".compare_models_discordant.pdf")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".compare_models_discordant.pdf"," file already exists, set --overwrite=TRUE to overwrite.")
    }
  } else if( opt$module == "predict" ){
    if( file.exists(paste0(output_dir_prefix,".predict.txt")) ){
      overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".predict.txt"," file already exists, set --overwrite=TRUE to overwrite.")
    }
  }
  if( !is.null(overwrite_message) ){
    stop(overwrite_message)
  } 
}

#Now run specified model
source(paste0(script.basename,"/modules/", opt$module, ".R"))


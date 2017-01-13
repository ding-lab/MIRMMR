#This script parses command line arguments and runs specified analysis module

suppressPackageStartupMessages(library("optparse"))

#Set options
option_list <- list(
  make_option(c("-m", "--module"), default=NULL, type="character", help="What module should be run (compare, penalized, stepwise, univariate), must be specified"),
  make_option(c("-d", "--data_frame"), default=NULL, type="character", help="File path to data frame, must be specified"),
  make_option(c("-i", "--msi_status"), default=NULL, type="character", help="Column name referring to MSI status (stored as binary vector), must be specified"),
  make_option(c("-c", "--first_data_column"), default=NULL, type="integer", help="Column number referring to first parameter in data frame, assuming all higher columns are also parameters, must be specified"),
  make_option(c("-o", "--output_prefix"), default=NULL, type="character", help="Output file name prefix, must be specified"),
  make_option(c("-d", "--output_directory"), default="", type="character", help="Output directory, must be specified"),
  make_option(c("--plots"), default=FALSE, type="logical", help="Produce informative plots in the penalized module, default=%default"),
  make_option(c("--alpha"), default=0.9, type="double", help="Penalized module: Parameter alpha used in glmnet::glmnet, default=%default"),
  make_option(c("--consensus"), default=FALSE, type="logical", help="Penalized module: Use consensus method in additional to best lambda method to find set of coefficients that appear in most models, only useful when used with --plots, default=%default"),
  make_option(c("--lambda"), default="lambda.1se", type="character", help="Penalized module: Parameter lambda used in glmnet::cv.glmnet, default=\"%default\", options: lambda.1se, lambda.min"),
  make_option(c("--nfolds"), default=10, type="integer", help="Penalized module: Parameter nfolds used in glmnet::cv.glmnet, default=%default"),
  make_option(c("--repeats"), default=1000, type="integer", help="Penalized module: Number of times to repeat testing to determine optimal lambda in penalized module, default=%default"),
  make_option(c("--set_seed"), default=FALSE, type=NULL, help="Penalized module: Option to set seed before penalized module, seed can be set to any number or TRUE (1), default=%default"),
  make_option(c("--train"), default=FALSE, type=NULL, help="Penalized module: Option to use a training set/test set approach to measure accuracy of best lambda approach, default=%default"),
  make_option(c("--train_proportion"), default=0.8, type="double", help="Penalized module: With --train=TRUE, the proportion of samples to keep in the training set, default=%default"),
  make_option(c("--type_measure"), default="class", type="character", help="Penalized module: Parameter type.measure used in glmnet::cv.glmnet (options: mse, deviance, mae, class, auc), default=%default")
)

#Retrieve command line arguments
opt <- parse_args(OptionParser(usage="%prog -m MODULE -d DATA_FRAME -i MSI_STATUS -c FIRST_DATA_COLUMN -o OUTPUT_PREFIX -d OUTPUT_DIRECTORY [options]\n\nAssumptions:\n1. The data frame has meta information columns (e.g. sample name, cancer type, MSI score)\n2. Followed by data columns (i.e. predictors in regression models)\n3. And that the input data frame has column headers.", option_list=option_list))

#Need these to refer to script directory MSImodel/modules
#http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script

#Source these helper functions
source("modules/helper_functions.R")

#Sanity check on the parameter inputs
sanity_checks(opt)

#Now run specified model
source(paste0("modules/", opt$module, ".R"))


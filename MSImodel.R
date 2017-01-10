#This script parses command line arguments and runs specified analysis module

suppressPackageStartupMessages(library("optparse"))

#set options
option_list <- list(
  make_option(c("-m", "--module"), default=NULL, type="character", help="What module should be run (compare, penalized, stepwise, univariate), must be specified"),
  make_option(c("-d", "--data_frame"), default=NULL, type="character", help="File path to data frame, must be specified"),
  make_option(c("-i", "--msi_score"), default=NULL, type=NULL, help="Column name or number referring to MSI score (numeric, or binary vector), must be specified"),
  make_option(c("-c", "--first_data_column"), default=NULL, type="integer", help="Column number referring to first parameter in data frame, assuming all higher columns are also parameters, must be specified"),
  make_option(c("-a", "--alpha"), default=0.9, type="double", help="Parameter alpha used in glmnet::glmnet, default=%default"),
  make_option(c("-f", "--nfolds"), default=10, type="integer", help="Parameter nfolds used in glmnet::cv.glmnet, default=%default"),
  make_option(c("-l", "--lambda"), default="lambda.1se", type="character", help="Parameter lambda used in glmnet::cv.glmnet, default=\"%default\", options: lambda.1se, lambda.min"),
  make_option(c("-n", "--number_repetitions"), default=1000, type="integer", help="Number of times to repeat testing to determine optimal lambda in penalized module, default=%default"),
  make_option(c("-p", "--train_proportion"), default=0.8, type="double", help="Proportion of samples to retain in training set in penalized module, default=%default"),
  make_option(c("-s", "--set_seed"), default=FALSE, type=NULL, help="Option to set seed before penalized module, seed can be set to any number, default=\"%default\""),
  make_option(c("-t", "--type_measure"), default="class", type="character", help="Parameter type.measure used in glmnet::cv.glmnet, default=\"%default\", options: mse, deviance, mae, class, auc"),
  make_option(c("-v", "--verbose"), default=TRUE, type="logical", help="Print extra output along the way, default=%default")
)

#retrieve command line arguments
opt <- parse_args(OptionParser(usage="%prog -m MODULE -d DATA_FRAME -i MSI_SCORE -c FIRST_DATA_COLUMN [options]\n\nMSImodel assumes the data frame has meta information columns (e.g. sample name, cancer type, MSI score), followed by data columns (i.e. parameters that will be analyzed)", option_list=option_list))

#How to print verbose output to stderr
if( opt$verbose ){
  write("Write some message here...\n", stderr())
}

#Run the specified module
#Need this is refer to script directory MSImodel/modules
#http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source(paste0("modules/", opt$module, ".R"))

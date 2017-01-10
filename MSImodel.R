#This script parses command line arguments and runs specified analysis module

suppressPackageStartupMessages(require("argparse"))

#create parser object
parser <- ArgumentParser()

#set options
#Related to all modules
parser$add_argument("-c", "--first_data_column", numeric) #column number of first model parameter (assumes all columns above this one are also model parameters)
parser$add_argument("-d", "--data_frame", ) #location of the data frame
parser$add_argument("-i", "--msi_score", ) #column name or number of the data frame refering to true MSI measurement
parser$add_argument("-m", "--module", action) #what module to run
#Related to penalized.R module
parser$add_argument("-a", "--alpha", default 0.9, alpha for cv.glmnet())
parser$add_argument("-f", "--folds", default 10, number of folder for cv.glmnet())
parser$add_argument("-l", "--lambda", default NONE, use this lambda instead of cv calculation)
parser$add_argument("-n", "--number_repititions", default 1000, number of times to repeat penalized testing)
parser$add_argument("-p", "--train_proportion", default 0.8, proportion of data to use as training set)
parser$add_argument("-s", "--set_seed", defaults to False, unless given numeric or TRUE alternative)
parser$add_argument("-t", "--type_measure", default class, which accuracy measure to use for cv.glmnet())


#retrieve command line arguments
args <- parser$parse_args()

#How to print verbose output to stderr
if( args$verbose ){
  write("Write some message here...\n", stderr())
}

#Run the specified module
#Need this is refer to script directory MSImodel/modules
#http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
run set seed if TRUE or numeric 
source(paste0("modules/", args$module, ".R"))

#This script parses command line arguments and runs specified analysis module

suppressPackageStartupMessages(require("argparse"))

#create parser object
parser <- ArgumentParser()

#set options
parser$add_argument("-m", "--module", action) #what module to run
parser$add_argument("-d", "--data_frame", ) #location of the data frame
parser$add_argument("-i", "--msi_score", ) #column name or number of the data frame refering to MSI measurement
parser$add_argument("-f", "--first_data_column", numeric) #column number of first model parameter (assumes all columns above this one are also model parameters)

#retrieve command line arguments
args <- parser$parse_args()

#How to print verbose output to stderr
if( args$verbose ){
  write("Write some message here...\n", stderr())
}

#Run the specified module
#Need this is refer to script directory MSImodel/modules
#http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
source(paste0("modules/", args$module, ".R"))

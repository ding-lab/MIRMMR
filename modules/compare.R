#Compare results from various MSI status prediction methods
suppressPackageStartupMessages(library(ggplot2))

binary_columns <- which(apply(df[,fdc:ncol(df)],2,function(x) length(levels(as.factor(x))))==2)
numeric_columns <- which(apply(df[,fdc:ncol(df)],2,function(x) length(levels(as.factor(x))))!=2)

plot_df <- data.frame(tpr=NULL,fpr=NULL,method=NULL)
for(i in numeric_columns){
  temp_df <- data.frame(truth=df[,col], score=df[,i+fdc-1])
  method <- names(df)[i+fdc-1]
  roc_df <- roc( temp_df )
  auc_measure <- auc( roc_df )
  plot_df <- rbind(plot_df, data.frame(tpr=roc_df[,1], fpr=roc_df[,2], method=rep(paste0(method, "\nAUC=", format(auc_measure,digits=4)),nrow(roc_df))))
}

plot_roc( plot_df, "ROC comparison of different methods" )

#discordant_df <- subset(df, !apply(data.frame(df[,col], df[,binary_columns+fdc-1), 1, all), )
#write.table(discordant_df, paste0(output_dir_prefix,".compare_models_discordant.txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

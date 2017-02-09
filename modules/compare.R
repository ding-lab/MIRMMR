#Compare results from various MSI status prediction methods
suppressPackageStartupMessages(library(ggplot2))

binary_columns <- which(apply(df[,fdc:ncol(df)],2,function(x) length(levels(as.factor(x))))==2)
numeric_columns <- which(apply(df[,fdc:ncol(df)],2,function(x) length(levels(as.factor(x))))!=2)

numeric_df <- df[,fdc+numeric_columns-1]
n_numeric <- ncol(numeric_df)

plot_df <- data.frame(tpr=NULL,fpr=NULL,method=NULL)
cutoffs = vector()
for(i in numeric_columns){
  temp_df <- data.frame(truth=df[,col], score=numeric_df[,i])
  method <- names(numeric_df)[i]
  roc_df <- roc( temp_df )
  cutoffs[i] <- roc_df[which(max(roc_df[,1]+1-roc_df[,2])==roc_df[,1]+1-roc_df[,2]),3]
  auc_measure <- auc( roc_df )
  plot_df <- rbind(plot_df, data.frame(tpr=roc_df[,1], fpr=roc_df[,2], method=rep(paste0(method, "\nAUC=", format(auc_measure,digits=4), "\nCutoff=", format(cutoffs[i],digits=4)),nrow(roc_df))))
}

plot_roc( plot_df, xlab="False positive rate", ylab="True positive rate", color_indicates="Method", theme=opt$theme_bw )

if(n_numeric>1){
  for(i in 1:(n_numeric-1)){
    for(j in (i+1):n_numeric){
      #if(i==j){
      #  next
      #}
      #compare column i and j
      x <- numeric_df[,i]
      y <- numeric_df[,j]
      xcall <- (x >= cutoffs[i])
      ycall <- (y >= cutoffs[j])
      dis12 <- as.numeric(xcall!=ycall)+1
      plot_df <- data.frame(truth=df[,col], x, y, shape=c("No","Yes")[dis12])
      names(plot_df)[2] <- names(numeric_df)[i]
      names(plot_df)[3] <- names(numeric_df)[j]
      plot_compare( plot_df, xlab=names(plot_df)[2], ylab=names(plot_df)[3], color_indicates=opt$color_indicates, theme=opt$theme_bw, xcutoff=cutoffs[i], ycutoff=cutoffs[j])
    }
  }
}


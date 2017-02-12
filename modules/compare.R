#Compare results from various MSI status prediction methods
suppressPackageStartupMessages(library(ggplot2))

binary_columns <- which(apply(df[,fdc:ncol(df)],2,function(x) length(levels(as.factor(x))))==2)
numeric_columns <- which(apply(df[,fdc:ncol(df)],2,function(x) length(levels(as.factor(x))))!=2)
binary_df <- as.data.frame(apply(as.data.frame(df[,fdc+binary_columns-1]), 2, function(x) as.logical(as.numeric(as.factor(x))-1)))
names(binary_df) <- names(binary_columns)
numeric_df <- as.data.frame(df[,fdc+numeric_columns-1])
n_numeric <- ncol(numeric_df)

plot_df <- data.frame(tpr=NULL,fpr=NULL,method=NULL)
cutoffs = vector()
for(i in 1:n_numeric){
  temp_df <- data.frame(truth=df[,col], score=numeric_df[,i])
  temp_df <- temp_df[complete.cases(temp_df),]
  method <- names(numeric_df)[i]
  roc_df <- roc( temp_df )
  cutoffs[i] <- roc_df[which(max(roc_df[,1]+1-roc_df[,2])==roc_df[,1]+1-roc_df[,2]),3]
  auc_measure <- auc( roc_df )
  plot_df <- rbind(plot_df, data.frame(tpr=roc_df[,1], fpr=roc_df[,2], method=rep(paste0("\n",method, "\nAUC=", format(auc_measure,digits=4), "\nCutoff=", format(cutoffs[i],digits=4),"\n"),nrow(roc_df))))
}

plot_roc( plot_df, xlab="False positive rate", ylab="True positive rate", color_indicates="Method", theme=opt$theme_bw )

numeric_as_binary_df <- rep(TRUE, nrow(numeric_df))
for(i in 1:n_numeric){
  numeric_as_binary_df <- data.frame(numeric_as_binary_df, numeric_df[,i] >= cutoffs[i])
  names(numeric_as_binary_df)[i] <- names(numeric_df)[i]
}
numeric_as_binary_df <- numeric_as_binary_df[,-1]

truth_table <- df[,col]
if( nrow(binary_df) > 0 ){
  truth_table <- data.frame(truth_table, binary_df)
}
if( ncol(numeric_as_binary_df) > 0 ){
  truth_table <- data.frame(truth_table, numeric_as_binary_df)
}
discordant_vector <- apply(truth_table,1,function(x) !(all(x,na.rm=T) | all(!x,na.rm=T)))
write.table(df[discordant_vector,], paste0(output_dir_prefix,".compare_discordant_samples.txt"),quote=F,row.names=F,col.names=T,sep="\t")


### Plot things

if(n_numeric>1){
  for(i in 1:(n_numeric-1)){
    for(j in (i+1):n_numeric){
      x <- numeric_df[,i]
      y <- numeric_df[,j]
      xcall <- (x >= cutoffs[i])
      ycall <- (y >= cutoffs[j])
      dis12 <- as.numeric(xcall!=ycall)+1
      nrow_df <- nrow(df)
      plot_df <- data.frame(truth=df[,col], x, y, discordant=c("No","Yes")[dis12])
      label_df <- data.frame(xpos=rep(NA, nrow_df), ypos=rep(NA, nrow_df), xtext=rep(NA, nrow_df), ytext=rep(NA, nrow_df))
      label_df <- label_df[complete.cases(plot_df),]
      if(nrow(label_df)>=2){
        label_df$xpos[1:2] <- c(cutoffs[i], max(x,na.rm=TRUE))
        label_df$ypos[1:2] <- c(min(y, na.rm=TRUE), cutoffs[j])
        label_df$xtext[1] <- format(cutoffs[i], digits=4)
        label_df$ytext[2] <- format(cutoffs[j], digits=4)
      }
      plot_df <- plot_df[complete.cases(plot_df),]
      plot_df <- data.frame(plot_df, label_df)
      names(plot_df)[2] <- names(numeric_df)[i]
      names(plot_df)[3] <- names(numeric_df)[j]
      method1=names(plot_df)[2]
      method2=names(plot_df)[3]
      if( !opt$overwrite & file.exists(paste0(output_dir_prefix,".compare_models_discordant.", method1, "-", method2, ".pdf")) ){
        overwrite_message <- paste0(overwrite_message,"\n",output_dir_prefix,".compare_models_discordant.", method1, "-", method2, ".pdf"," file already exists, set --overwrite=TRUE to overwrite.")
        stop(overwrite_message)
      }
      plot_compare( plot_df, xlab=method1, ylab=method2, color_indicates=opt$color_indicates, theme=opt$theme_bw, xcutoff=cutoffs[i], ycutoff=cutoffs[j])
    }
  }
}


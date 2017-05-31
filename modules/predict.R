#Module returns model predictions of new data

model <- get(load(opt$model))
if( all( names(model$coefficients)[-1] %in% names(df)[fdc:ncol(df)] )){
  output <- data.frame(df[,1:(fdc-1)], MIRMMR=predict(model, newdata=df[,fdc:ncol(df)], type="response"), df[, names(df) %in% names(model$coefficients)[-1] ])
  write.table(output, file=paste0(output_dir_prefix,".predict.txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

} else{
  stop("Error: There are terms in the model that do not exist in the new data frame.")
}


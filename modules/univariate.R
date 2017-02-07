#Module returns univariate model p-values for each parameter in data frame

#Vector of parameters
params <- names(df)[fdc:ncol(df)]
#store results here
results <- vector("list", length(params))
names(results) <- params

for( p in params ){
  model <- NULL
  #Run univariate logistic regression
  #col is the MSI status column name
  #p is the parameter of interest column name
  results[[p]] <- glm( df[,col] ~ df[,p], data=df, family=binomial)
}

#Write summary table
output <- matrix(NA,length(params)+1,9)
output[1,] <- c("Parameter",paste0("Intercept",c("_Estimate","_Std._Error", "_z_value", "_P(>|z|)")),paste0("Parameter",c("_Estimate","_Std._Error", "_z_value", "_P(>|z|)")))
count <- 1
for( p in params ){
  count <- count + 1
  output[count,] <- c(p,summary(results[[p]])$coefficients[1,],summary(results[[p]])$coefficients[2,])
}
write.table(output,file=paste0(output_dir_prefix,".univariate_summary.txt"),quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#Write models to output 
save(results, file=paste0(output_dir_prefix,".univariate_models.Robj"))

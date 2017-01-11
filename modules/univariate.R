#Module returns univariate model p-values for each parameter in data frame
col <- opt$msi_score
fdc <- opt$first_data_column

#Check if MSI score column is binary (use logistic regression) or not (use linear regression)
if( length(levels(droplevels(as.factor(df[,col]))))==2 ){
  msi_binary <- TRUE
} else if( length(levels(droplevels(as.factor(df[,col]))))>2 ){
  msi_binary <- FALSE
}
else{
  stop("Whoa, MSI score column should have more than one level.\n") 
}

#Vector of parameters
params <- names(df)[fdc:ncol(df)]
#store results here
results <- vector("list", length(params))
names(results) <- params

for( p in params ){
  model <- NULL
  #Run univariate logistic or linear regression
  if( msi_binary ){ #logistic
    model <- glm( df[,col] ~ df[,p], data=df, family="logistic")
  } else if( !msi_binary ){ #linear
    model <- lm( df[,col] ~ df[,p], data=df)
  }
  results[[p]] <- summary(model)
}

return(results)

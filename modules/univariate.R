#Module returns univariate model p-values for each parameter in data frame
#col <- opt$msi_score
#fdc <- opt$first_data_column
#msi_binary is TRUE if only two factors in MSI score column, FALSE otherwise

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

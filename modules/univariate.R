#Module returns univariate model p-values for each parameter in data frame

#Vector of parameters
params <- names(df)[fdc:ncol(df)]
#store results here
results <- vector("list", length(params))
names(results) <- params

for( p in params ){
  model <- NULL
  #Run univariate logistic or linear regression
  model <- glm( df[,col] ~ df[,p], data=df, family="logistic")
  results[[p]] <- summary(model)
}

#Write results to output 
return(results)

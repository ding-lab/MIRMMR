#Module returns univariate model p-values for each parameter in data frame
msi_column
first_parameter

#Pseudo code
results_df <- EMPTY
parameter_names <- df[,first_paramter:end]
if( msi_column type is 0/1 or TRUE/FALSE or two factors ){
  for(name in parameter_names){
    check name is numeric
    results_df <- summary(glm(msi_column ~ name, data=df, family="binomial"))$p.value
  }
}
else{
  for(name in parameter_names){
    check name is numeric
    results_df <- summary(lm(msi_column ~ name, data=df))$p.value
  }  
}
print or write(results_df)

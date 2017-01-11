#Module returns best model using stepwise regression (MASS::stepAIC)
#col <- opt$msi_score
#fdc <- opt$first_data_column
#msi_binary is TRUE if only two factors in MSI score column, FALSE otherwise

suppressPackageStartupMessage(library(MASS))

#Model formula
f <- paste( names(df)[col], paste(names(df)[fdc:ncol(df)], collapse=" + "), sep=" ~ ")

#Fit full model
if( msi_binary ){ #logistic
  model <- glm(f, data=df, family="logistic")
} else if( !msi_binary ){ #linear
  model <-  lm(f, data=df)
}

#Find best model by AIC
stepwise_model <- stepAIC(model, direction="both")

return(stepwise_model)


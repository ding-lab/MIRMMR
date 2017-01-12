#Module returns best model using stepwise regression (MASS::stepAIC)

suppressPackageStartupMessage(library(MASS))

#Model formula
f <- paste( names(df)[col], paste(names(df)[fdc:ncol(df)], collapse=" + "), sep=" ~ ")

#Fit full model
model <- glm(f, data=df, family="logistic")

#Find best model by AIC
stepwise_model <- stepAIC(model, direction="both")

#TODO: Predict MSI-H probabilities based on model
#Calculate AUC, plot ROC, etc. 
#Write model summary to output
#Save model in .Robj file


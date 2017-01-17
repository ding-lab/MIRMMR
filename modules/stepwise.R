#Module returns best model using stepwise regression (MASS::stepAIC)

#Load MASS package
suppressPackageStartupMessage(library(MASS))

#Model formula
f <- paste( names(df)[col], paste(names(df)[fdc:ncol(df)], collapse=" + "), sep=" ~ ")

#Fit full model
model <- glm(f, data=df, family=binomial)

#Find best model by AIC
stepwise_model <- stepAIC(model, direction="both")

#Predict MSI-H probabilities based on model
probs <- predict(stepwise_model, newdata=df, type="response")

#ROC
roc_df <- roc( data.frame(df[,col], probs) )

#AUC
auc <- auc( roc_df )

#Plot ROC
p <- plot_roc( roc_df, user_title="ROC curve of stepwise regression model")

#Write model summary to output
write(summary(stepwise_model), file=paste0(output_dir_prefix,".stepwise_model_summary.txt"))

#Save model in .Robj file
save(stepwise_model, file=paste0(output_dir_prefix,".stepwise_model.Robj"))

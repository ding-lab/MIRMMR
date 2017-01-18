#Module returns best model using stepwise regression (MASS::stepAIC)

#Load MASS package
suppressPackageStartupMessages(library(MASS))

#Model formula
f <- as.formula( paste( col, paste(names(df)[fdc:ncol(df)], collapse=" + "), sep=" ~ "))

#Fit full model
model <- glm(f, data=df, family=binomial)

#Find best model by AIC
all_models <- stepAIC(model, direction="both")
best_model <- glm(all_models$formula, data=df, family=binomial)

#Predict MSI-H probabilities based on model
#probs <- predict(best_model, newdata=df, type="response")

#ROC
#roc_df <- roc( data.frame(df[,col], probs) )

#AUC
#auc <- auc( roc_df )

#Write model summary to output
sink(paste0(output_dir_prefix,".stepwise_model_summary.txt"))
print(summary(best_model))
sink()

#Save model in .Robj file
save(best_model, file=paste0(output_dir_prefix,".stepwise_model.Robj"))

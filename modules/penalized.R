#Module performs penalized regression

#Set seed for repeatable results (cross-validation approach is stochastic)
if(opt$s != FALSE){
  set.seed(opt$s)
}

suppressPackageStartupMessage(library(glmnet))

#Data to be used
if( opt$train ){
  train_set <- sample(c(TRUE,FALSE), size=nrow(df), replace=TRUE, prob=c(opt$train_proportion, 1-opt$train_proportion))
} else{
  train_set <- rep(TRUE, nrow(df))
}

trainX <- as.matrix(df[train,fdc:ncol(df)])
trainY <- as.numeric(df[train,col])

if( opt$train ){
  testX <- as.matrix(df[!train,fdc:ncol(df)])
  testY <- as.matrix(df[!train,fdc:ncol(df)])
}

# If user specified to run consensus step
if( opt$consensus ){
  parameter_counts <- consensus_parameters(opt, trainX, trainY)
}

# Find best model using penalized regression
best_model <- best_lambda_model(opt, trainX, trainY)

# Use test data to evaluate training model
if( opt$train ){
  test_roc <- test_accuracy(best_model, testX, testY)
  test_auc <- auc( test_roc )
  #TODO return something
}

############PLOT THINGS

#Plot to visualize how many times each variable appeared
vars <- as.matrix(unique(sort(variables_in_model)))
num <- apply(vars,1,function(y) sum(variables_in_model==y))
plot_df <- data.frame(vars, num, best_coefs=(vars %in% coefs_in_best_model))
p <- ggplot(plot_df, aes(x=num, y=reorder(vars,num), color=best_coefs)) + geom_point(size=3) + labs(x="Number of models", y="Model variable", title="Number of models each variable appeared in", color="Included in\n'best' model")
pdf("plots/variables_in_model.pdf",10,10)
print(p)
dev.off()

model_df <- xcc_tcga_msi[,names(xcc_tcga_msi) %in% c("tcga_msih", coefs_in_best_model)]
glm.bestmodel.fit <- glm( tcga_msih ~ ., data=model_df, family=binomial)
summary(glm.bestmodel.fit)

plot_df <- data.frame(predicted=glm.bestmodel.fit$fitted.values, actual=xcc_tcga_msi$tcga_msih, cancer_type=xcc_tcga_msi$cancer_type)
p <- ggplot(plot_df) + aes(x=actual, y=predicted, color=cancer_type) + geom_jitter(width=0.5, height=0) + labs(x="TCGA MSI-H status", y="Fitted probabilty of being MSI-H", title="Probabilty of MSI-H", color="Cancer type")
pdf("plots/probability_of_msih_dots.pdf",10,10)
print(p)
dev.off()
p <- ggplot(plot_df) + aes(x=actual, y=predicted) + geom_violin() + labs(x="TCGA MSI-H status", y="Fitted probabilty of being MSI-H", title="Probabilty of MSI-H")
pdf("plots/probability_of_msih_violin.pdf",10,10)
print(p)
dev.off()

#Calculate ROC curves for model and MSIsensor
num_obs = length(msisensor_range)
model_range <- sort(glm.bestmodel.fit$fitted.values)
msisensor_range <- sort(xcc_tcga_msi$msisensor)
model_roc <- matrix(NA,num_obs,3)
msisensor_roc <- matrix(NA,num_obs,3)
for(i in 1:num_obs){
  TT = sum(xcc_tcga_msi$tcga_msih & glm.bestmodel.fit$fitted.values >= model_range[i])
  TF = sum(xcc_tcga_msi$tcga_msih & !(glm.bestmodel.fit$fitted.values >= model_range[i]))
  FT = sum(!xcc_tcga_msi$tcga_msih & glm.bestmodel.fit$fitted.values >= model_range[i])
  FF = sum(!xcc_tcga_msi$tcga_msih & !(glm.bestmodel.fit$fitted.values >= model_range[i]))
  model_tab <- matrix(c(TT, FT, TF, FF),2,2)
  TT = sum(xcc_tcga_msi$tcga_msih & xcc_tcga_msi$msisensor > msisensor_range[i])
  TF = sum(xcc_tcga_msi$tcga_msih & !(xcc_tcga_msi$msisensor > msisensor_range[i]))
  FT = sum(!xcc_tcga_msi$tcga_msih & xcc_tcga_msi$msisensor > msisensor_range[i])
  FF = sum(!xcc_tcga_msi$tcga_msih & !(xcc_tcga_msi$msisensor > msisensor_range[i]))
  msisensor_tab <- matrix(c(TT, FT, TF, FF),2,2)
  #TPR, FPR for each type
  model_roc[i,] <- c(model_range[i], model_tab[1,1]/sum(model_tab[1,]), model_tab[2,1]/sum(model_tab[2,]))
  msisensor_roc[i,] <- c(msisensor_range[i], msisensor_tab[1,1]/sum(msisensor_tab[1,]), msisensor_tab[2,1]/sum(msisensor_tab[2,]))
}
plot_df <- rbind(data.frame(value = model_roc[,1], tpr = model_roc[,2], fpr = model_roc[,3], method=rep("Model",num_obs)), data.frame(value = msisensor_roc[,1], tpr = msisensor_roc[,2], fpr = msisensor_roc[,3], method=rep("MSIsensor",num_obs)))
p <- ggplot(plot_df) + aes(x=fpr, y=tpr, color=method) + geom_line(size=1) + geom_abline(intercept=0, slope=1, linetype=2, alpha=0.5) + labs(x="False positive rate", y="True positive rate", title="ROC curve of regression model vs. MSIsensor", color="Method")
pdf("plots/roc_curve_model_msisensor.pdf",10,10)
print(p)
dev.off()

#Module performs penalized regression
#col <- opt$msi_score
#fdc <- opt$first_data_column
#msi_binary is TRUE if only two factors in MSI score column, FALSE otherwise

#Set seed for repeatable results (cross-validation approach is stochastic)
if(opt$s != FALSE){
set.seed(opt$s)
}

suppressPackageStartupMessage(library(glmnet))

#Data to be used
myX <- as.matrix(df[,fdc:ncol(df)])
myY <- as.numeric(df[,col])

#Function +++++++++++++++++++++++++++++
#Track how many times each variable appears in models
consensus_parameters <- function(opt, msi_binary, myX, myY){
  variables_in_model <- vector("list", opt$n)
  for(i in 1:opt$n){
    if( msi_binary ){ 
      cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
    } else{
      #FIT linear model
    }
    variables_in_model[[i]] <- row.names(coef(cvfit, s=opt$l))[coef(cvfit, s=opt$l)[,1] != 0]
  }
  variables_in_model_df <- data.frame(table(unlist(variables_in_model)))
  names(variables_in_model_df) <- c("Parameter","Count")
  return(variables_in_model_df)
}

#Function +++++++++++++++++++++++++++++
#Find best lambda and return model with best lambda coefficients 
#Modified from stackexchange user 'Sideshow Bob', answer from 2016/03/31. Thanks!
#http://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results/173895#173895
best_lambda_model <- function(opt, msi_binary, myX, myY){
  lambdas = NULL
  for (i in 1:opt$n){
    if( msi_binary ){
      cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
    } else{
      #FIT linear model
    }
    #TODO: optimize this with prespecified data frame or list structure for lambdas
    errors <- data.frame(cvfit$lambda,cvfit$cvm)
    lambdas <- rbind(lambdas,errors)
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[,2], list(lambdas$cvfit.lambda), mean)

  # select the best one
  bestindex = which(lambdas[,2]==min(lambdas[,2]))
  bestlambda = lambdas[bestindex,1]

  # and now run glmnet once more with it
  if( msi_binary ){
    glmnet.fit <- glmnet(x=myX, y=myY, lambda=bestlambda, family="binomial", alpha=opt$a)
  } else{
    #FIT linear model
  }

  # Determine the coeffients in model derived from best lambda
  best_coef_indicator <- coef(glmnet.fit)[,1] != 0
  coefs_in_best_model <- rownames(coef(glmnet.fit))[best_coef_indicator]
  
  # Fit a model lm() or glm() using coefficients from above
  model_df <- data.frame(myY, myX[, names(myX) %in% coefs_in_best_model])
  names(model_df)[1] <- "outcome"
  if( msi_binary ){
    model <- glm( outcome ~ ., data=model_df, family=binomial)
  } else{
    model <- lm( outcome ~ ., data=model_df)
  }
  return(model)
}

#Function +++++++++++++++++++++++++++++
#Return many models fit using best lambda approach but fit with training data
train_lambda_models <- function(opt, msi_binary, myX, myY){
  train_models <- vector("list", opt$repeat_tests)
  for(i in 1:opt$repeat_tests){
    train <- sample(c(TRUE,FALSE), size=nrow(myX), replace=TRUE, prob=c(opt$p,1-opt$p))
    trainX <- myX[train,]
    trainY <- myY[train]
    test_df <- data.frame(myY, myX)[!train]
    names(test_df)[1] <- "outcome"
    train_model <- best_lambda_model(opt, msi_binary, trainX, trainY)
    test_predictions <- predict(train_model, newdata=test_df, type="response")
    test_auc <- calculate_auc(myY[!train],test_predictions)
    train_models[[i]] <- list(train_model, test_auc)
  }
}

#Function +++++++++++++++++++++++++++++
#Make comparison between best model using all data and models derived from training data

# End of functions ++++++++++++++++++++
# If user specified to run consensus step
if( opt$consensus ){
  parameter_counts <- consensus_parameters(opt, msi_binary, myX, myY)
}

# Find best model using penalized regression
best_model <- best_lambda_model(opt, msi_binary, myX, myY)

# Run testing section if specified
train_models <- test_lambda_models(opt, msi_binary, myX, myY)





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

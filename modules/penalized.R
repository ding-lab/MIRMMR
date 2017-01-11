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

#Keep track of which variables in each model
variables_in_model <- vector("list", opt$n)
for(i in 1:opt$n){
  print(i)
  if( msi_binary){ 
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
  } else{
    #FIT linear model
  }
  variables_in_model[[i]] <- row.names(coef(cvfit, s=opt$l))[coef(cvfit, s=opt$l)[,1] != 0]
}

#Find mean lambda
#Modified from stackexchange user 'Sideshow Bob', answer from 2016/03/31.
#http://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results/173895#173895
lambdas = NULL
for (i in 1:opt$n){
  print(i)
  if( msi_binary ){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
  } else{
    #FIT linear model
  }
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
return(glmnet.fit)
best_coef_indicator <- coef(glmnet.bestLambda.fit)[,1] != 0
coefs_in_best_model <- rownames(coef(glmnet.bestLambda.fit))[best_coef_indicator]

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

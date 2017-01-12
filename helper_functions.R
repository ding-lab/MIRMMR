#Set of functions to help facilitate rest of scripts

#Function +++++++++++++++++++++++++++++
#Check inputs for sanity and correctness, exit if any problem
sanity_checks <- function(opt){
  #TODO check all input parameters
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Track how many times each variable appears in models
consensus_parameters <- function(opt, myX, myY){
  #TODO update opt parameter names
  variables_in_model <- vector("list", opt$n)
  for(i in 1:opt$n){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
    variables_in_model[[i]] <- row.names(coef(cvfit, s=opt$l))[coef(cvfit, s=opt$l)[,1] != 0]
  }
  variables_in_model_df <- data.frame(table(unlist(variables_in_model)))
  names(variables_in_model_df) <- c("Parameter","Count")
  return(variables_in_model_df)
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Find best lambda and return model with best lambda coefficients 
#Modified from stackexchange user 'Sideshow Bob', answer from 2016/03/31. Thanks!
#http://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results/173895#173895
best_lambda_model <- function(opt, myX, myY){
  #TODO update opt parameter names
  lambdas = NULL
  for (i in 1:opt$n){
    cvfit <- cv.glmnet(x=myX, y=myY, nfolds=opt$f, family="binomial", type.measure=opt$t, alpha=opt$a)
    #TODO: optimize this with prespecified data frame or list structure for lambdas
    errors <- data.frame(cvfit$lambda,cvfit$cvm)
    lambdas <- rbind(lambdas,errors)
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[,2], list(lambdas$cvfit.lambda), mean)

  # select the best lambda
  bestindex = which(lambdas[,2]==min(lambdas[,2]))
  bestlambda = lambdas[bestindex,1]

  # and now run glmnet once more with the best lambda
  glmnet.fit <- glmnet(x=myX, y=myY, lambda=bestlambda, family="binomial", alpha=opt$a)

  # Determine the coeffients in model derived from best lambda
  best_coef_indicator <- coef(glmnet.fit)[,1] != 0
  coefs_in_best_model <- rownames(coef(glmnet.fit))[best_coef_indicator]
  
  # Fit a glm model using coefficients from above
  model_df <- data.frame(myY, myX[, names(myX) %in% coefs_in_best_model])
  names(model_df)[1] <- "outcome"
  model <- glm( outcome ~ ., data=model_df, family=binomial)
  return(model)
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Calculate points in a ROC curve
roc <- function(){
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Calculate area under a ROC curve (AUC)
auc <- function(){
}
#--------------------------------------

#Function +++++++++++++++++++++++++++++
#Plot a ROC curve
plot_roc <- function(){
}
#--------------------------------------

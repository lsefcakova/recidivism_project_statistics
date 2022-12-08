#########################################################################################
##
## R ROUTINES FOR SEMNINAR 1 - STATISTICAL MODELLING AND INFERENCE.
##
## AUTHOR: DAVID ROSSELL, PAUL ROGNON, UNIVERSITAT POMPEU FABRA
## REDONE FOR BINOMIAL FAMILY
#########################################################################################

# INDEX
# 1. SETTING PENALIZATION PARAMETER VIA BIC
# 2. CROSS-VALIDATION
# 3. LASSO POST-SELECTION INFERENCE
# 4. ROUTINES TO SIMULATE DATA

#########################################################################################
## 1. SETTING PENALIZATION PARAMETER VIA BIC
#########################################################################################

lasso.bic <- function(y,x,extended=FALSE) {
  #Select model in LASSO path with best BIC (using LASSO regression estimates)
  #Input
  # - y: vector with response variable
  # - x: design matrix
  #
  #Output: list with the following elements
  # - coef: LASSO-estimated regression coefficient with lambda set via BIC
  # - ypred: predicted y
  # - lambda.opt: optimal value of lambda
  # - lambda: data.frame with bic and number of selected variables for each value of lambda
  require(glmnet)
  fit <- glmnet(x=x,y=y,family='binomial',alpha=1)
  pred <- cbind(1,x) %*% rbind(fit$a0,fit$beta)
  n <- length(y)
  p <- colSums(fit$beta!=0) + 1
  if (!extended){
    bic <- deviance(fit) + log(n)*p 
  } else {
    bic <- deviance(fit) + log(n)*p + 2*log(choose(ncol(x),p))
  }
  
  sel <- which.min(bic)
  beta <- c(fit$a0[sel],fit$beta[,sel]); names(beta)[1]= 'Intercept'
  ypred <- pred[,sel]
  ans <- list(coef=beta,ypred=ypred,lambda.opt=fit$lambda[sel],lambda=data.frame(lambda=fit$lambda,bic=bic,nvars=p))
  return(ans)
}





#########################################################################################
## 2. CROSS-VALIDATION
#########################################################################################


kfoldCV.lasso <- function(y,x,K=10,seed,criterion='cv') {
  ## Perform K-fold cross-validation for LASSO regression estimate (lambda set either via cross-val or BIC or EBIC)
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - seed: random number generator seed (optional)
  ## - criterion: the criterion to select the penalization parameter, either cross-val or BIC or EBIC
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  require(glmnet)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  cat("Starting cross-validation")
  if (ncol(x)>0) {  #if there are some covariates
    for (k in 1:K) {
      sel <- subset==k
      if (criterion=='cv') {
        fit <- cv.glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], alpha = 1, nfolds=10,family = 'binomial',type_measure='auc')
        b= as.vector(coef(fit,s='lambda.min'))
        pred[sel] <- b[1] + x[sel,,drop=FALSE] %*% as.matrix(b[-1])
      } else if (criterion=='bic'){
        fit <- lasso.bic(y=y[!sel],x=x[!sel,,drop=FALSE])
        pred[sel] <- fit$coef[1] + x[sel,,drop=FALSE] %*% matrix(fit$coef[-1],ncol=1)
      } else if (criterion=='ebic'){
        fit <- lasso.bic(y=y[!sel],x=x[!sel,,drop=FALSE],extended = TRUE)
        pred[sel] <- fit$coef[1] + x[sel,,drop=FALSE] %*% matrix(fit$coef[-1],ncol=1) 
      } else { stop("method.lambda not implemented") }
      cat(".")
    }
  } else { #if there are no covariates, just use the intercept
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  cat("\n")
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}
kfoldCV.mle <- function(y,x,K=10,seed) {
  ## Perform K-fold cross-validation for least-squares regression estimate
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - seed: random number generator seed (optional)
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  if (ncol(x)>0) {
    for (k in 1:K) {
      sel <- subset==k
      fit <- glmnet(x=x[!sel,,drop=FALSE], y=y[!sel],lambda=0,family = 'binomial')
      pred[sel] <- predict(fit, newx=x[sel,,drop=FALSE])
    }
  } else {
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  return(list(pred=pred))
}



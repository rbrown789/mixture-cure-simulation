
###############################################################################
# This code defines all the functions that will be used in Simulation 1.  These include
# functions generating data from the cure rate mixture model, as well as the functions
# to run the models, and the outermost functions for running the simulation study itself
#
###############################################################################



### Some functions ###

logit <- function(x) { return(log(x/(1-x))) }  
expit <- function(x) {  return( exp(x)/ ( 1+exp(x) ) ) } # the expit function (inverse of logit)



findint <- function(x,y,value)
{
  ###################################################################
  # x and y are vectors of equal length (y must be numeric), value must be a number
  # returns x entry corresponding to y entry closest to "value"
  ###################################################################
  
  if( length(y) != length(x)) {stop("x and y must be of equal length")}
  if( !is.numeric(y) ) { stop("y must be numeric") }
  yy <- abs(y-value)
  return(x[which.min(yy)])
}



rlnorm.trunc <- function(n,meanlog,sdlog,trunc.pt=NULL)
{
  #########################################################################################
  #  Generate samples from right-truncated log normal distirbution
  #    - trunc.pt is the truncation point, if trunc.pt=NULL, then generates from regular lognormal
  #    - uses rejection sampling
  #    - is partially vectorized over meanlog and sdlog ( can input a vector of length n or a single value)
  #########################################################################################
  
  if( !is.numeric(trunc.pt) & !is.null(trunc.pt) ) {  stop("trunc.pt must be numeric or NULL")}
  out <- rlnorm(n,meanlog,sdlog)
  
  if(is.numeric(trunc.pt))
  {
    if(trunc.pt <=0) { stop("trunc.pt must be greater than 0") }
    
    req <- (out > trunc.pt)
    
    repeat
    {
      if(length(meanlog)==1) { meanlog.new <- meanlog }  else { meanlog.new <- meanlog[req] }
      if(length(sdlog)==1) { sdlog.new <- sdlog }  else { sdlog.new <- sdlog[req] }
      
      out[req] <- rlnorm( length(out[req]),meanlog.new,sdlog.new )
      req <- (out > trunc.pt)
      if( sum(req) == 0 ) {break}
    }
  }
  
  return(out)
}


rweibull.trunc <- function(n,shape,scale,trunc.pt=NULL)
{
  #########################################################################################
  #  Generate samples from right-truncated weibull distirbution
  #    - trunc.pt is the truncation point, if trunc.pt=NULL, then generates from regular weibull
  #    - uses rejection sampling
  #    - is partially vectorized over shape and scale ( can input a vector of length n or a single value)
  #########################################################################################
  
  if( !is.numeric(trunc.pt) & !is.null(trunc.pt) ) {  stop("trunc.pt must be numeric or NULL") }
  
  out <- rweibull(n,shape,scale)
  
  
  if(is.numeric(trunc.pt))
  {
    if(trunc.pt <=0) { stop("trunc.pt must be greater than 0") }
    
    req <- (out > trunc.pt)
    
    repeat
    {
      if(length(shape)==1) { shape.new <- shape }  else { shape.new <- shape[req] }
      if(length(scale)==1) { scale.new <- scale }  else { scale.new <- scale[req] }
      
      out[req] <- rweibull( length(out[req]),shape.new,scale.new )
      req <- (out > trunc.pt)
      if( sum(req) == 0 ) {break}
    }
  }
  
  return(out)
}


rllogis.trunc <- function(n,location,scale,trunc.pt=NULL)
{
  #########################################################################################
  #  Generate samples from right-truncated log logistic distirbution
  #    - trunc.pt is the truncation point, if trunc.pt=NULL, then generates from regular loglogistic
  #    - uses rejection sampling
  #    - is partially vectorized over location and scale ( can input a vector of length n or a single value)
  #########################################################################################
  
  if( !is.numeric(trunc.pt) & !is.null(trunc.pt) ) {  stop("trunc.pt must be numeric or NULL") }
  out <- exp( rlogis(n,location,scale) )
  
  req <- (out > trunc.pt)
  
  if(is.numeric(trunc.pt))
  {
    if(trunc.pt <=0) { stop("trunc.pt must be greater than 0") }
    
    repeat
    {
      if(length(location)==1) { location.new <- location }  else { location.new <- location[req] }
      if(length(scale)==1) { scale.new <- scale }  else { scale.new <- scale[req] }
      
      out[req] <- exp( rlogis( length(out[req]),location.new,scale.new ) )
      req <- (out > trunc.pt)
      if( sum(req) == 0 ) {break}
    }
  }
  
  return(out)
}




mixture.sim <- function( design, beta.ber, binary.link="logit",              
                         beta.surv,surv.dist="lognormal",scale,trunc.pt=NULL, 
                         fixed.cens=120  )
{
  ##########################################################################################
  # Function to simulate data from a cure-survival time mixture model with covariates and fixed censoring
  #
  #  - design: the design (X) matrix, with an intercept column of 1s (note that this determines sample size)
  #
  #  - beta.ber: the vector of coefficients associated with non-detection rate
  #
  #  - binary.link is link function used to map the linear predictor to the detection probability
  #
  #  - beta.surv: the vector of coefficients associated with survival time, conditional upon detection
  #
  #  - surv.dist: the survival time distribution type, must be one of c('lognormal','loglogistic','weibull')
  #
  #  - scale: the scale parameter for the survival time distribution used in survreg()
  #      note: survreg's scale = 1/(rweibull shape), survreg's intercept = log(rweibull scale)
  # 
  #  - fixed.cens: the fixed "censoring" time point    
  #
  ##########################################################################################
  
  # j <- 4
  # nsamp <- 90
  # x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
  # x2 <- rep(0:1,nsamp/2)
  # design <- cbind(1,x1,x2)
  # beta.ber <- c(intmat$ber.int[j],beta.ber1,beta.ber2)
  # binary.link <- "logit"
  # beta.surv <- c(intmat$weibull.int[j],beta.surv1,beta.surv2)
  # surv.dist <- "weibull"
  # scale <- 0.5
  # trunc.pt <- NULL
  # fixed.cens <- 120
  
  # simulate probabilities of detection
  if(binary.link=="logit")  
  {  p <- as.numeric( expit( design %*% beta.ber ) ) 
  
  } else if(binary.link=="probit") 
  {  p <- as.numeric( pnorm( design %*% beta.ber ) ) 
  
  } else { stop("binary.link must be one of c('logit','probit') ")}
  
  ##########
  # simulate detections
  y <- rbinom( n=length(p),size=1,prob=p )
  
  ##########
  # simulate detection times of those which are detected 
  det <- y==0
  t <- rep(fixed.cens,length(y))
  logmean <- design[det,] %*% beta.surv
  
  
  if( surv.dist=="lognormal") 
  { t[det] <- rlnorm.trunc( n=length(y[det]), meanlog=logmean,sdlog=scale,trunc.pt=trunc.pt ) 
  
  } else if( surv.dist=="weibull") 
  { t[det] <- rweibull.trunc( n=length(y[det]), shape= 1/scale, scale=exp(logmean),trunc.pt=trunc.pt ) 
  
  } else if(surv.dist=="loglogistic")
  { t[det] <- rllogis.trunc( n=length(y[det]), location=logmean, scale=scale,trunc.pt=trunc.pt )  
  
  } else  { stop("surv.dist must be one of c('lognormal','loglogistic','weibull') ")}
  
  
  # generate the dataframe with the predictors and the survival time and censoring indicator
  out.df <- data.frame( design)
  names(out.df) <- c("int","x1","x2")
  out.df$t <- t
  out.df$t[out.df$t > 120] <- 120
  out.df$delta <- 1*(!out.df$t==120)
  
  
  return( out.df )
}


survestform <- function( est )
{
  for(i in names(est)) { est[[i]] <- as.numeric(est[[i]])}
  est <- data.frame(est)
  return(est)
}




findsurv <- function(mat,value)
{
  ####################################################################################
  # mat is a matrix with survival function in first column and times in second column
  # value is the survival time, returns associated survival function value 
  ##################################################################################
  
  if(value <= 0) { stop("survival times must be greater than 0") }
  
  if( sum(mat[,2] <= value)==0 ) { return(1) 
  } else {  return( mat[sum(mat[,2] <= value),1] ) }
}


findtime <- function(mat,value)
{
  ####################################################################################
  # mat is a matrix with survival function in first column and times in second column
  # value is the quantile of survival time, returns time associated with value
  ##################################################################################
  
  if(  value > 1 | value < 0  ) { stop("quantiles must be in (0,1)") }
  
  if( sum(mat[,1] >= value)==0 ) { return(NA) 
  } else {  return( mat[sum(mat[,1] >= value),2] ) }
}


predsmcure <- function( object,newdata, times, quants){
  
  #######################################################################################
  # This function calculates predicted cure probabilities, survival functions, and quantiles
  # of survival time based on a fitted smcure model object.
  #
  # object: smcure object, of either "aft" or "ph" type 
  # newdata: dataframe with new data, assumption currently is that latency and incidence use same covariates
  #           and thus only one newdata input is needed
  # times: times for which to predict survival function
  # quants: quantiles for which to predict time
  #
  # utilizes smcure packages built-in predict function predictsmcure()
  #######################################################################################
  
  
  # object <- mod
  # newdata <- nd
  # times <- c(30,60)
  # quants <- 0.5
  
  model <- object$modtype
  predout <- predictsmcure(object,newX=newdata,newZ=newdata,model=object$modtype )
  
  
  pred <- predout$prediction
  
  uncureprob <- t( predout$newuncureprob )
  cureprob <- 1- uncureprob
  
  survs <- NULL
  tms <- NULL
  
  if(model=="aft"){
    bb <- (ncol(pred)/2)
    for(i in 1:bb )
    {
      tmp <- pred[order(pred[,i+bb]),c(i,i+ bb)]
      survs <- rbind( survs,sapply( times, FUN=findsurv, mat=tmp ) )
      tms <- rbind( tms,sapply( quants, FUN=findtime,mat=tmp) )
    }
  }
  
  if(model=="ph"){
    ncoll <- ncol(pred)
    for(i in 1:(ncoll-1) )
    {
      tmp <- pred[order(pred[,ncoll]),c(i,ncoll)]
      survs <- rbind( survs,sapply( times, FUN=findsurv, mat=tmp ) )
      tms <- rbind( tms,sapply( quants, FUN=findtime,mat=tmp) )
    }
  }
  
  colnames(survs) <- times
  colnames(tms) <- quants
  
  out <- list(newdata=newdata,cureprob=cureprob,quantile=tms,survival=survs  )
  
  return(out)
}


smcure.novar <- function (formula, cureform, offset = NULL, data, na.action = na.omit, 
                          model = c("aft", "ph"), link = "logit", Var = FALSE, emmax = 50, 
                          eps = 1e-07,savedata=TRUE, showout=FALSE) 
################################################################################
# custom smcure function.  only does point estimation and saves the dataset if desired
##############################################################################

{
  call <- match.call()
  model <- match.arg(model)
  if(showout) { cat("Program is running..be patient...") }
  data <- na.action(data)
  n <- dim(data)[1]
  mf <- model.frame(formula, data)
  cvars <- all.vars(cureform)
  Z <- as.matrix(cbind(rep(1, n), data[, cvars]))
  colnames(Z) <- c("(Intercept)", cvars)
  if (!is.null(offset)) {
    offsetvar <- all.vars(offset)
    offsetvar <- data[, offsetvar]
  }
  else offsetvar <- NULL
  Y <- model.extract(mf, "response")
  X <- model.matrix(attr(mf, "terms"), mf)
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  Time <- Y[, 1]
  Status <- Y[, 2]
  bnm <- colnames(Z)
  nb <- ncol(Z)
  if (model == "ph") {
    betanm <- colnames(X)[-1]
    nbeta <- ncol(X) - 1
  }
  if (model == "aft") {
    betanm <- colnames(X)
    nbeta <- ncol(X)
  }
  w <- Status
  b <- eval(parse(text = paste("glm", "(", "w~Z[,-1]", ",family = quasibinomial(link='", 
                               link, "'", ")", ")", sep = "")))$coef
  if (model == "ph") 
    beta <- coxph(Surv(Time, Status) ~ X[, -1] + offset(log(w)), 
                  subset = w != 0, method = "breslow")$coef
  if (model == "aft") 
    beta <- survreg(Surv(Time, Status) ~ X[, -1])$coef
  emfit <- smcure:::em(Time, Status, X, Z, offsetvar, b, beta, model, 
                       link, emmax, eps)
  b <- emfit$b
  beta <- emfit$latencyfit
  s <- emfit$Survival
  logistfit <- emfit$logistfit
  
  fit <- list()
  class(fit) <- c("smcure")
  fit$logistfit <- logistfit
  fit$b <- b
  fit$beta <- beta
  
  if(showout) { cat(" done.\n") }
  fit$call <- call
  fit$bnm <- bnm
  fit$betanm <- betanm
  fit$s <- s
  fit$Time <- Time
  if (model == "aft") {
    error <- drop(log(Time) - beta %*% t(X))
    fit$error <- error
  }
  
  if(savedata) { fit$data <- data }
  fit$offset <- offset
  fit$na.action <- na.action
  fit$modtype <- model
  fit$tau <- emfit$tau
  
  if(showout) { printsmcure(fit,FALSE) }
  
  return(fit)
}





smcure.bootvar <- function( object,emmax=50,eps=1e-07,nboot=100, newdata,
                            pred.cureprob=T,pred.times=NULL,pred.quants=0.5,
                            ci=TRUE,clevel=0.95) 
{
  ################################################################################################
  # Non-parametric bootstrap variance estimates for coefficient estimates, and, optionally,
  #  for predicted cure probabilities, survival function values, and quantiles of survival time
  #
  #   object: the object from smcure.novar(savedata=T), this will define the dataset, model, formula,
  #      as well as point estimates for all the coefficients
  #  
  #  Currently no functionality for NULL pred.times and pred.quants, so these must be supplied
  #
  ##############################################################################################
  
  # object <- mod
  # pred.cureprob <- T
  # pred.times <- c(30,60)
  # pred.quants <- 0.5
  # nboot <- 500
  # emmax <- 50
  # eps <- 1e-07
  # newdata <- nd
  # ci <- TRUE
  # clevel <- 0.95
  
  # check if object is of type smcure 
  call <- match.call()
  if (!inherits(object, "smcure")) 
    stop("Object must be results of smcure")
  
  # grab all the model information from the model object
  model <- as.character(object$call["model"])
  form <- as.formula( as.character( object$call["formula"] ) )
  cureform <- as.formula( as.character( object$call["cureform"] ) )
  link <- as.character(object$call["link"])
  data <- object$data
  offset <- object$offset
  na.action <- object$na.action
  
  nb <- length(object$b)
  nbeta <- length(object$beta)
  bnm <- object$bnm
  betanm <- object$betanm
  
  
  # some other information
  call <- match.call()
  data <- na.action(data)
  n <- dim(data)[1]
  mf <- model.frame(form, data)
  Y <- model.extract(mf, "response")
  
  # coefficient estimates
  beta <- object$beta
  if(model=="aft") { names(beta) <- c("(Intercept)",all.vars(form)[-(1:2)]) }
  if(model=="ph") { names(beta) <- all.vars(form)[-(1:2)] }
  
  b <- object$b
  names(b) <- c("(Intercept)",all.vars(cureform))
  
  
  # generate predictions for the model
  ndnms <- apply(newdata, 1, paste, collapse="")
  survnms <- rep( paste0("p",as.character(pred.times) ),each=nrow(newdata))
  quantnms <- rep( paste0("t",as.character(pred.quants*100) ),each=nrow(newdata))
  
  oripred <- predsmcure( object, newdata, quants=pred.quants,times=pred.times)
  cp <- c(oripred$cureprob)
  names(cp) <- paste0("cp_",ndnms)
  
  surv <- c(oripred$survival)
  names(surv) <- paste0(survnms,"_",rep(ndnms,length(pred.times)))
  
  quant <- c(oripred$quantile)
  names(quant) <- paste0( quantnms,"_",rep(ndnms,length(pred.quants)))
  
  
  
  # Initialize all bootstrap matrices
  if (model == "ph") {
    b_boot <- matrix(rep(0, nboot * nb), nrow = nboot)
    beta_boot <- matrix(rep(0, nboot * (nbeta)), nrow = nboot)
    iter <- matrix(rep(0, nboot), ncol = 1)
  }
  
  if (model == "aft") {
    b_boot <- matrix(rep(0, nboot * nb), nrow = nboot)
    beta_boot <- matrix(rep(0, nboot * (nbeta)), nrow = nboot)
  }
  
  
  if(!is.null(newdata)) {
    if( pred.cureprob==T  ) {
      cp_boot <- matrix(rep(0,nboot*nrow(newdata)),nrow=nboot)
    }
    
    if( !is.null(pred.times)  ){
      surv_boot <- matrix(rep(0,nboot*nrow(newdata)*length(pred.times) ),nrow=nboot)
    }
    
    if( !is.null(pred.quants) ){
      quants_boot <- matrix(rep(0,nboot*nrow(newdata)*length(pred.quants) ),nrow=nboot)
    }
  }
  
  
  
  ### boostrapping loop (stratified sampling of censored and non-censored observations to 
  ### maintain the correct ratio)
  data1 <- data[Y[,2]==1,]
  data0 <- data[Y[,2]==0,]
  n1 <- nrow(data1)
  n0 <- nrow(data0)
  
  i <- 1
  while (i <= nboot) {
    
    
    # create bootstrap dataset
    id1 <- sample(1:n1, n1, replace = TRUE)
    id0 <- sample(1:n0, n0, replace = TRUE)
    bootdata <- rbind(data1[id1,],data0[id0,])
    
    
    # fit model and use built-in predict function
    bootfit <- tryCatch ( smcure.novar(form, cureform, offset , bootdata, na.action, model,
                                       link, Var = FALSE, emmax,eps,savedata=F ), 
                          warning=function(w) { return("no convergence") } )
    
    if( class(bootfit) =="smcure" ) {
      bootpred <- predsmcure( bootfit, newdata, quants=pred.quants,times=pred.times)
      
      # store the estimated coefficients and predictions
      b_boot[i,] <- bootfit$b
      beta_boot[i,] <- bootfit$beta
      cp_boot[i,] <-  t(matrix(bootpred$cureprob))
      surv_boot[i,] <- t(matrix(bootpred$survival))  
      quants_boot[i,] <- t(matrix(bootpred$quantile))
      
      if (bootfit$tau < eps) {  i <- i + 1 }
    } 
  }
  
  # calcualte variance and stadard deviations
  b_var <- apply(b_boot, 2, var)
  names(b_var) <- names(b)
  
  beta_var <- apply(beta_boot, 2, var)
  names(beta_var) <- names(beta)
  
  b_sd <- sqrt(b_var)
  beta_sd <- sqrt(beta_var)
  
  cp_var <- apply(cp_boot,2,var)
  names(cp_var) <- names(cp)
  cp_sd <- sqrt(cp_var)
  
  surv_var <- apply(surv_boot,2,var)
  names(surv_var) <- names(surv)
  surv_sd <- sqrt(surv_var)
  
  quants_var <- apply(quants_boot,2,var)
  names(quants_var) <- names(quant)
  quants_sd <- sqrt(quants_var)
  
  
  if(ci) {
    
    # calculate CIs here based on asymptotic normality of MLEs and bootstrap variance estimates
    # this ignores the variability in the boostrap variance estimate (which should be relatively small,
    # assuming enough bootstrap samples)
    crit <- qnorm( (1-clevel)/2,lower.tail=F)
    
    beta_upper <- beta + beta_sd*crit
    beta_lower <- beta - beta_sd*crit
    
    b_upper <- b + b_sd*crit
    b_lower <- b - b_sd*crit
    
    cp_upper <- cp + cp_sd*crit
    cp_upper <- ifelse(cp_upper > 1,1,cp_upper)
    cp_lower <- cp - cp_sd*crit
    cp_lower <- ifelse(cp_lower < 0,0,cp_lower)
    
    surv_upper <- surv + surv_sd*crit
    surv_upper <- ifelse(surv_upper > 1,1,surv_upper)
    surv_lower <- surv - surv_sd*crit
    surv_lower <- ifelse(surv_lower < 0,0,surv_lower)
    
    quant_upper <- quant + quants_sd*crit
    quant_upper <- ifelse(quant_upper < 0,0,quant_upper)
    quant_lower <- quant - quants_sd*crit
    quant_lower <- ifelse(quant_lower < 0,0,quant_lower)
  }
  
  out <- list()
  class(out) <- "smcurevar"
  
  out$beta <- beta
  out$beta_var <- beta_var
  if(ci) { out$beta_upper <- beta_upper; out$beta_lower <- beta_lower}
  
  out$b <- b
  out$b_var <- b_var
  if(ci) { out$b_upper <- b_upper; out$b_lower <- b_lower}
  
  out$predict_order <- ndnms
  
  out$cp <- cp
  out$cp_var <- cp_var
  if(ci) { out$cp_upper <- cp_upper; out$cp_lower <- cp_lower}
  
  out$surv <- surv
  out$surv_var <- surv_var
  if(ci) { out$surv_upper <- surv_upper; out$surv_lower <- surv_lower}
  
  out$quant <- quant
  out$quant_var <- quants_var
  if(ci) { out$quant_upper <- quant_upper; out$quant_lower <- quant_lower}
  
  return(out)
  
}  


runmodels <- function( dataset,mcaft.boot=100,mcph.boot=250 )
{
  
  ## we have to do some weird scope/environment stuff here because rms package uses global options
  ## variable for datadist
  on.exit(detach("design.options"))
  attach(list(), name="design.options")
  x1 <- dataset$x1
  x2 <- dataset$x2
  assign('dd', datadist( x1,x2), pos='design.options')
  options(datadist="dd")
  
  rownms <- c("b1","b2", paste0("t50_",c("00","10","01","11")),
              paste0("p30_",c("00","10","01","11")),paste0("p60_",c("00","10","01","11")),
              paste0("curep_",c("00","10","01","11")),"cureb1","cureb2")
  colnms <- c("est","lwr95","upr95")
  nd <- data.frame(x1=c(0,1,0,1),x2=c(0,0,1,1))
  
  ## Cox model 
  mod <- cph( Surv(t,delta) ~ x1 + x2, data=dataset, x=TRUE, y=TRUE,surv=T )
  med <- Quantile(mod)
  prob <- Survival(mod)
  
  betas <- cbind( mod$coefficients,matrix(confint(mod),nrow=2) )
  t50 <- as.matrix( Predict(mod,x1,x2,fun=function(x) med(lp=x)) )[,c(3,5,4)]
  p30 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=30)) )[,c(3,5,4)]
  p60 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=60)) )[,c(3,5,4)]
  
  cox.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(cox.mat) <- colnms
  rownames(cox.mat) <- rownms 
  
  
  ## parametric lognormal
  mod <- psm( Surv(t,delta) ~ x1 + x2, data=dataset, dist="lognormal")
  med <- Quantile(mod)
  prob <- Survival(mod)
  
  betas <- cbind( mod$coefficients[2:3],matrix(confint(mod)[2:3,],nrow=2) )
  t50 <- as.matrix( Predict(mod,x1,x2,fun=function(x) med(lp=x)) )[,c(3,4,5)]
  p30 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=30)) )[,c(3,4,5)]
  p60 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=60)) )[,c(3,4,5)]
  
  ln.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(ln.mat) <- colnms
  rownames(ln.mat) <- rownms 
  
  
  ## parametric weibull
  mod <- psm( Surv(t,delta) ~ x1 + x2, data=dataset, dist="weibull")
  med <- Quantile(mod)
  prob <- Survival(mod)
  
  betas <- cbind( mod$coefficients[2:3],matrix(confint(mod)[2:3,],nrow=2) )
  t50 <- as.matrix( Predict(mod,x1,x2,fun=function(x) med(lp=x)) )[,c(3,4,5)]
  p30 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=30)) )[,c(3,4,5)]
  p60 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=60)) )[,c(3,4,5)]
  
  wb.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(wb.mat) <- colnms
  rownames(wb.mat) <- rownms 
  
  
  ## Parametric generalized gamma
  mod <- flexsurvreg( Surv(t,delta) ~ x1 + x2, data=dataset, dist="gengamma")
  betas <- cbind( mod$coefficients[4:5],matrix(confint(mod)[4:5,],nrow=2) )
  
  med <- function(mu,sigma,Q,P) { qgengamma(0.5, mu = mu, sigma = sigma, Q=Q,lower.tail = TRUE, log.p = FALSE) }
  summmed <- summary(mod, newdata=nd,fn = med, t=1,B = 1000)
  t50 <- matrix(unlist(summmed),ncol=4,byrow=T)[,2:4]
  
  summp30 <- summary(mod, newdata=nd,type="survival", t=30,B = 1000)
  p30 <- matrix(unlist(summp30),ncol=4,byrow=T)[,2:4]
  
  summp60 <- summary(mod, newdata=nd,type="survival", t=60,B = 1000)
  p60 <- matrix(unlist(summp60),ncol=4,byrow=T)[,2:4]
  
  gg.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(gg.mat) <- colnms
  rownames(gg.mat) <- rownms 
  
  
  ## SemiParametric AFT Mixture-cure
  mod <- smcure.novar( Surv(t,delta) ~ x1 + x2,cureform= ~ x1 + x2,data=dataset,model="aft",link="logit",showout=F)
  modci <- smcure.bootvar( mod,newdata=nd,pred.times=c(30,60),pred.quants=0.5,nboot=mcaft.boot)
  
  betas <- cbind( modci$beta[-1],modci$beta_lower[-1],modci$beta_upper[-1] )
  t50 <- cbind( modci$quant,modci$quant_lower,modci$quant_upper)
  p3060 <- cbind(modci$surv,modci$surv_lower,modci$surv_upper)
  curep <- cbind(modci$cp,modci$cp_lower,modci$cp_upper)
  curebs <- cbind(modci$b[-1],modci$b_lower[-1],modci$b_upper[-1])
  
  mcaft.mat <- rbind(betas,t50,p3060,curep,curebs)
  colnames(mcaft.mat) <- colnms
  rownames(mcaft.mat) <- rownms 
  
  
  
  ## SemiParametric PH Mixture-cure
  mod <- smcure.novar( Surv(t,delta) ~ x1 + x2,cureform= ~ x1 + x2,data=dataset,model="ph",link="logit",showout=F)
  modci <- smcure.bootvar( mod,newdata=nd,pred.times=c(30,60),pred.quants=0.5,nboot=mcph.boot)
  
  betas <- cbind( modci$beta,modci$beta_lower,modci$beta_upper )
  t50 <- cbind( modci$quant,modci$quant_lower,modci$quant_upper)
  p3060 <- cbind(modci$surv,modci$surv_lower,modci$surv_upper)
  curep <- cbind(modci$cp,modci$cp_lower,modci$cp_upper)
  curebs <- cbind(modci$b[-1],modci$b_lower[-1],modci$b_upper[-1])
  
  mcph.mat <- rbind(betas,t50,p3060,curep,curebs)
  colnames(mcph.mat) <- colnms
  rownames(mcph.mat) <- rownms 
  
  # abind into a 3D array and apply dimnames to the relevant dimension
  out <- abind( cox.mat, ln.mat, wb.mat, gg.mat,mcaft.mat,mcph.mat, along=3 )
  dimnames(out)[[3]] <- c("CoxPH","ParLN","ParWB","ParGG","MCAFT","MCPH")
  
  return(out)
}


runmodels_nocr <- function(dataset)
{
  ## we have to do some weird scope/environment stuff here because rms package uses global options
  ## variable for datadist
  on.exit(detach("design.options"))
  attach(list(), name="design.options")
  x1 <- dataset$x1
  x2 <- dataset$x2
  assign('dd', datadist( x1,x2), pos='design.options')
  options(datadist="dd")
  
  rownms <- c("b1","b2", paste0("t50_",c("00","10","01","11")),
              paste0("p30_",c("00","10","01","11")),paste0("p60_",c("00","10","01","11")),
              paste0("curep_",c("00","10","01","11")),"cureb1","cureb2")
  colnms <- c("est","lwr95","upr95")
  nd <- data.frame(x1=c(0,1,0,1),x2=c(0,0,1,1))
  
  ## Cox model 
  mod <- cph( Surv(t,delta) ~ x1 + x2, data=dataset, x=TRUE, y=TRUE,surv=T )
  med <- Quantile(mod)
  prob <- Survival(mod)
  
  betas <- cbind( mod$coefficients,matrix(confint(mod),nrow=2) )
  t50 <- as.matrix( Predict(mod,x1,x2,fun=function(x) med(lp=x)) )[,c(3,5,4)]
  p30 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=30)) )[,c(3,5,4)]
  p60 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=60)) )[,c(3,5,4)]
  
  cox.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(cox.mat) <- colnms
  rownames(cox.mat) <- rownms 
  
  
  ## parametric lognormal
  mod <- psm( Surv(t,delta) ~ x1 + x2, data=dataset, dist="lognormal")
  med <- Quantile(mod)
  prob <- Survival(mod)
  
  betas <- cbind( mod$coefficients[2:3],matrix(confint(mod)[2:3,],nrow=2) )
  t50 <- as.matrix( Predict(mod,x1,x2,fun=function(x) med(lp=x)) )[,c(3,4,5)]
  p30 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=30)) )[,c(3,4,5)]
  p60 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=60)) )[,c(3,4,5)]
  
  ln.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(ln.mat) <- colnms
  rownames(ln.mat) <- rownms 
  
  
  ## parametric weibull
  mod <- psm( Surv(t,delta) ~ x1 + x2, data=dataset, dist="weibull")
  med <- Quantile(mod)
  prob <- Survival(mod)
  
  betas <- cbind( mod$coefficients[2:3],matrix(confint(mod)[2:3,],nrow=2) )
  t50 <- as.matrix( Predict(mod,x1,x2,fun=function(x) med(lp=x)) )[,c(3,4,5)]
  p30 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=30)) )[,c(3,4,5)]
  p60 <- as.matrix( Predict(mod,x1,x2,fun=function(x) prob(lp=x,times=60)) )[,c(3,4,5)]
  
  wb.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(wb.mat) <- colnms
  rownames(wb.mat) <- rownms 
  
  
  ## Parametric generalized gamma
  mod <- flexsurvreg( Surv(t,delta) ~ x1 + x2, data=dataset, dist="gengamma")
  betas <- cbind( mod$coefficients[4:5],matrix(confint(mod)[4:5,],nrow=2) )
  
  med <- function(mu,sigma,Q,P) { qgengamma(0.5, mu = mu, sigma = sigma, Q=Q,lower.tail = TRUE, log.p = FALSE) }
  summmed <- summary(mod, newdata=nd,fn = med, t=1,B = 1000)
  t50 <- matrix(unlist(summmed),ncol=4,byrow=T)[,2:4]
  
  summp30 <- summary(mod, newdata=nd,type="survival", t=30,B = 1000)
  p30 <- matrix(unlist(summp30),ncol=4,byrow=T)[,2:4]
  
  summp60 <- summary(mod, newdata=nd,type="survival", t=60,B = 1000)
  p60 <- matrix(unlist(summp60),ncol=4,byrow=T)[,2:4]
  
  gg.mat <- rbind(betas,t50,p30,p60,matrix(NA,nrow=nrow(nd)+ncol(nd),ncol=3))
  colnames(gg.mat) <- colnms
  rownames(gg.mat) <- rownms 
  
  # abind into a 3D array and apply dimnames to the relevant dimension
  out <- abind( cox.mat, ln.mat, wb.mat, gg.mat, along=3 )
  dimnames(out)[[3]] <- c("CoxPH","ParLN","ParWB","ParGG")
  
  return(out)
}



runsim <- function(crmodels=F,design,beta.ber,binary.link,beta.surv,surv.dist,scale,fixed.cens,trunc.pt)
{
  ##################################################################################
  # Run one iteration of a simluation wih the given parameters
  ####################################################################################
  
  
  out.j <- "error/warning"
  
  while(!is.array(out.j)){
    
    # simulate dataset
    dat <- mixture.sim(   design=design 
                          ,beta.ber = beta.ber
                          ,binary.link=binary.link
                          ,beta.surv = beta.surv
                          ,surv.dist = surv.dist
                          ,scale=scale
                          ,fixed.cens=fixed.cens
                          ,trunc.pt=trunc.pt
    )
    
    # run models, output is 3-dim array of model results
    # if time has gone over an hour, throw an error
    # if there are errors (including time error) return character: "error/warning"
    if(crmodels==F){
      out.j <- tryCatch( runmodels_nocr(dat),
                         error=function(e){ "error/warning" },
                         warning=function(w){"error/warning"}   )
    }
    
    if(crmodels==T){
      out.j <- tryCatch( runmodels(dat),
                         error=function(e){ "error/warning" },
                         warning=function(w){"error/warning"}   )
    }
    
  } 
  
  return(out.j)
  
}



run_allscen <- function( Sdists,Ns,UORs,CPs,intmat,scen.mat,nsim=3,crmodels=F )
{
  
  # nsim <- 3
  # crmodels <- F
  
  # initialize the array
  modnms1 <- c("CoxPH","ParLN","ParWB","ParGG")
  modnms2 <- c(modnms1,"MCAFT","MCPH")
  modnms <- modnms1
  if(crmodels==T) {modnms <- modnms2}
  
  quantnms <- c("b1","b2","t50_00","t50_10","t50_01","t50_11","p30_00","p30_10","p30_01","p30_11",
                "p60_00","p60_10","p60_01","p60_11","curep_00","curep_10","curep_01","curep_11","cureb1","cureb2")
  
  pntci <- c("est","lwr95","upr95")  
  
  
  
  out.array <- array( NA, 
                      dim=c(length(Ns),length(UORs),length(CPs),length(Sdists),length(quantnms),
                            length(pntci),length(modnms),nsim ),
                      
                      dimnames=list( samp_size= paste0("N",Ns),
                                     unobs_rate=paste0("UOR",UORs),
                                     cens_prop=paste0("CP",CPs),
                                     surv_dist=paste0("SD",Sdists),
                                     pred_quant=quantnms,
                                     point_ci=pntci,
                                     models=modnms,
                                     simiter=paste0("sim",1:nsim)
                      )
  ) 
  
  
  scen.mat <- scen.mat[scen.mat$surv.dist %in% Sdists,]
  
  for( i in 1:nrow(scen.mat) ) # for each simulation scenario
  {
    base::cat(paste0("scen",i," "))
    
    # set the simulation scenario settings
    nsamp <- scen.mat$n[i]
    uorate <- scen.mat$unobs.rate[i]
    cprop <- scen.mat$censor.prop[i]
    sdist <- as.character( scen.mat$surv.dist[i] )
    tpoint <- NULL
    if(cprop==0) {tpoint <- 120} 
    
    x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
    x2 <- rep(0:1,nsamp/2)
    design <- cbind(1,x1,x2)
    
    crit <- intmat$unobs.rate==uorate & intmat$censor.prop==cprop 
    berint <- intmat$ber.int[crit]
    survint <- ifelse(sdist=="lognormal",intmat$lognorm.int[crit],intmat$weibull.int[crit])
    
    # set output array indices based on simulation settings
    a1 <- ifelse(nsamp==40,"N40","N90")
    a2 <- ifelse(uorate==0.15,"UOR0.15","UOR0.4")
    a3 <- ifelse(cprop==0,"CP0",ifelse(cprop==1/3,"CP0.333333333333333",ifelse(cprop==2/3,"CP0.666666666666667","CP1")))
    a4 <- ifelse(sdist=="lognormal","SDlognormal","SDweibull")
    
    
    # run the simulation nsim times
    out <- sapply(integer(nsim), function(i) { runsim(crmodels=crmodels
                                                      ,design=design 
                                                      ,beta.ber = c(berint,beta.ber1,beta.ber2)
                                                      ,binary.link="logit"
                                                      ,beta.surv = c(survint,beta.surv1,beta.surv2)
                                                      ,surv.dist = sdist
                                                      ,scale=0.5
                                                      ,fixed.cens=120
                                                      ,trunc.pt=tpoint) } 
                  , simplify="array")
    
    
    # store the output 
    out.array[a1,a2,a3,a4, , , , ] <- out
  }
  
  return(out.array)
}


#########################################################################################################


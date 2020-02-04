


######################################################################################
# This script runs the simulation study used to generate measure bias, variance, MSE,
# comparing parametric survival models, Cox PH, and semi-parametric AFT and PH cure models.
# First small simulation are run to calculate the intercept values needed to achieve
# the deisred nominal unobserved rate.  Then the actual simulations are run and saved as 
# an .Rdata image.  
#
# NOTE: Simulation time for code in its current form was approximately 1 week.
######################################################################################

library(compiler)
enableJIT(3)

library(survival)
library(rms)
library(flexsurv)
library(smcure)
library(abind)
library(R.utils)


root <- "C:/Users/rbrow/OneDrive/Documents/Public Github/mixture-cure-simulation/"
data <- "results/"
code <- "Code/"
graphics <- "Graphics/"

source(  paste0(root,code,"simulation_functions.R") )




### Set the Fixed Simulation Parameters ###
beta.ber1 <- 1.6
beta.ber2 <- 1
beta.surv1 <- 0.5
beta.surv2 <- 0.8
scale <- 0.5
fixed.cens <- 120

##########################################################################################################

set.seed(123)

########### Non-detection/Censoring curves for a range of intercept values #############

# non-detection rate for by logistic regression intercept values
nsamp <- 10000
x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
x2 <- rep(0:1,nsamp/2)
design <- cbind(1,x1,x2)

ints.b <- seq(-6,4,by=0.001)
pndet <- rep(NA,length(ints.b))

for(i in 1:length(ints.b)) {
  
  beta.ber <- c(ints.b[i],beta.ber1,beta.ber2)
  y <- rbinom( n=nsamp,size=1,prob=as.numeric( expit( design %*% beta.ber )) )
  pndet[i] <- mean(y)
}

ndet.int <- data.frame( intercept=ints.b,pr_ndet=pndet)
plot(ndet.int$intercept,ndet.int$pr_ndet,type="l",xlab="Intercept",ylab="Pr(Non-Detection)",ylim=c(0,1))



###########

# censoring rate (in detected sample) by the intercept values (for lognormal)
nsamp <- 10000
x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
x2 <- rep(0:1,nsamp/2)
design <- cbind(1,x1,x2)


ints <- seq(2,6,by=0.001)
pgt120_ln <- rep(NA,length(ints))

for(i in 1:length(ints)) {
  beta.surv <- c(ints[i],beta.surv1,beta.surv2)
  xx <- rlnorm( n=nsamp,meanlog=design %*% beta.surv ,sdlog=scale )
  pgt120_ln[i] <- sum(xx>fixed.cens)/length(xx)
}

lognorm.int <- data.frame( intercept=ints,pgt120=pgt120_ln)
plot(lognorm.int$intercept,lognorm.int$pgt120,type="l",xlab="Intercept",ylab="Pr(X > 120)",ylim=c(0,1))



###########

# censoring rate (in detected sample) by the intercept values (for weibull)
pgt120_wb <- rep(NA,length(ints))

for(i in 1:length(ints)) {
  beta.surv <- c(ints[i],beta.surv1,beta.surv2)
  xx <- rweibull( n=nsamp,shape=1/scale, scale=exp(design %*% beta.surv))
  pgt120_wb[i] <- sum(xx>fixed.cens)/length(xx)
}

weibull.int <- data.frame( intercept=ints,pgt120=pgt120_wb)
plot(weibull.int$intercept,weibull.int$pgt120,type="l",xlab="Intercept",ylab="Pr(X > 120)",ylim=c(0,1))


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


# generate the simulation scenario matrix
Ns <- c(40,90)                  # sample sizes 
UORs <- c(0.15,0.4)             # nominal unobserved (censoring + nondetection) rates
CPs <- c(0,1/3,2/3,1)          # censoring proportion (of unobserved) 
Sdists <- c("lognormal","weibull") # survival distributions

scen.mat <- expand.grid(n=Ns,unobs.rate=UORs,censor.prop=CPs,surv.dist=Sdists)



##########################################################################################

#### Generate Matrix of intercept values corresponding to desired levels of non-detection/censoring,
#### given the fixed parameters specified at the top
nms <- c("unobs.rate","censor.prop")
intmat <- scen.mat[!duplicated(scen.mat[,nms]),nms]
intmat$cens.rate <- with(intmat,unobs.rate*censor.prop) # censoring rate in the full dataset
intmat$nd.rate <- with(intmat,unobs.rate-cens.rate) # non-detection rate in the full dataset
intmat$cens.rate.detect <- with(intmat, cens.rate/(1-nd.rate) ) # censoring rate in JUST the detected observations

intmat$ber.int <- sapply( X=intmat$nd.rate, findint,x=ndet.int$intercept,y=ndet.int$pr_ndet)
intmat$weibull.int <- sapply( X=intmat$cens.rate.detect, findint, x=weibull.int$intercept,y=weibull.int$pgt120)
intmat$lognorm.int <- sapply( X=intmat$cens.rate.detect, findint, x=lognorm.int$intercept,y=lognorm.int$pgt120)


# for zero censoring cases (i.e. all unobserved detections are "cured"),  the zero censoring rate will 
# be primarily managed by truncation, but realistically, we'd have a distribution where probability of 
# detection at the fixed censoring time is getting very close to zero, say 1% 
intmat$weibull.int[intmat$cens.rate==0] <- findint( x=weibull.int$intercept,y=weibull.int$pgt120,value=0.01)
intmat$lognorm.int[intmat$cens.rate==0] <- findint( x=lognorm.int$intercept,y=lognorm.int$pgt120,value=0.01)


# Key for intmat Column Names  
# unobs.rate:  total unobserved rate in dataset (censoring + nondetection)
# censor.prop: proportion of unobs.rate that is due to censoring
# cens.rate: censoring rate in the full data (including nondetections)
# nd.rate: nondetection rate in the full dataset
# cens.rate.detect: censoring rate for detections only (this is what matters when simulating)
# ber.int: calculated intercept for logistic regression model, 
# weibull.int: calculated intercept for weibull survival model
# lognorm.int: calculated intercept for lognormal survival model

save.image(file=paste0(root,data,"simsetting.Rdata"))

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################


#### Let's do some timing tests, one iteration for each simulation setting
# 
# 
# #### All combinations of N, unobserved rate, and censoring proportion, 45 total runs 
# scen.mat3 <-  scen.mat[scen.mat$surv.dist=="lognormal",]
# 
# set.seed(1234)
# times3 <- rep(NA,nrow(scen.mat3))
# 
# i <- 1
# while(i <= nrow(scen.mat3) )
# {
#   
#   nsamp <- scen.mat3$n[i]
#   uorate <- scen.mat3$unobs.rate[i]
#   cprop <- scen.mat3$censor.prop[i]
#   sdist <- scen.mat3$surv.dist[i]
#   tpoint <- NULL
#   if(cprop==0) {tpoint <- 120} 
#   
#   x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
#   x2 <- rep(0:1,nsamp/2)
#   design <- cbind(1,x1,x2)
#   
#   crit <- intmat$unobs.rate==uorate & intmat$censor.prop==cprop 
#   berint <- intmat$ber.int[crit]
#   survint <- intmat$lognorm.int[crit]
#   
#   dat <- mixture.sim(  design=design 
#                        ,beta.ber = c(berint,beta.ber1,beta.ber2)
#                        ,binary.link="logit"
#                        ,beta.surv = c(survint,beta.surv1,beta.surv2)
#                        ,surv.dist = "lognormal"
#                        ,scale=0.5
#                        ,fixed.cens=120
#                        ,trunc.pt=tpoint
#   )
#   
#   time <- tryCatch( withTimeout( expr= system.time(runmodels(dat))[3],timeout=3600,onTimeout="error"),
#                     error=function(e){ "error/warning" },
#                     warning=function(w){"error/warning"} )
#   
#   if(time!="error/warning") { 
#     
#     base::cat(paste0(i," "))
#     times3[i] <- time
#     i <- i+1
#   }
#   
# }
# 
# scen.mat3$times3 <- times3
# sum(scen.mat3$times3)/60
# ( (sum(scen.mat3$times3))*2*100 )/(60*60*24)
# 
# 
# 
# 
# #### All combinations of N, unobserved rate, and censoring proportion, 45 total runs, without the cure models
# scen.mat3b <-  scen.mat[scen.mat$surv.dist=="lognormal",]
# 
# set.seed(123)
# times3b <- rep(NA,nrow(scen.mat3b))
# 
# i <- 1
# while(i <= nrow(scen.mat3b) )
# {
#   
#   nsamp <- scen.mat3b$n[i]
#   uorate <- scen.mat3b$unobs.rate[i]
#   cprop <- scen.mat3b$censor.prop[i]
#   sdist <- scen.mat3b$surv.dist[i]
#   tpoint <- NULL
#   if(cprop==0) {tpoint <- 120} 
#   
#   x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
#   x2 <- rep(0:1,nsamp/2)
#   design <- cbind(1,x1,x2)
#   
#   crit <- intmat$unobs.rate==uorate & intmat$censor.prop==cprop 
#   berint <- intmat$ber.int[crit]
#   survint <- intmat$lognorm.int[crit]
#   
#   dat <- mixture.sim(  design=design 
#                        ,beta.ber = c(berint,beta.ber1,beta.ber2)
#                        ,binary.link="logit"
#                        ,beta.surv = c(survint,beta.surv1,beta.surv2)
#                        ,surv.dist = "lognormal"
#                        ,scale=0.5
#                        ,fixed.cens=120
#                        ,trunc.pt=tpoint
#   )
#   
#   time <- tryCatch( withTimeout( expr= system.time(runmodels_nocr(dat))[3],timeout=15*60,onTimeout="error"),
#                     error=function(e){ "error/warning" },
#                     warning=function(w){"error/warning"} )
#   
#   if(time!="error/warning") { 
#     
#     base::cat(paste0(i," "))
#     times3b[i] <- time
#     i <- i+1
#   }
#   
# }
# 
# scen.mat3b$times3b <- times3b
# ( (sum(scen.mat3b$times3b))*2*100 )/(60*60*24)
# 
# ############################################################################################
# ############################################################################################
# ############################################################################################
############################################################################################

# j <- 4
# nsamp <- 40
# x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
# x2 <- rep(0:1,nsamp/2)
# design <- cbind(1,x1,x2)
# set.seed(123)
# dat <- mixture.sim(  design=design 
#                      ,beta.ber = c(intmat$ber.int[j],beta.ber1,beta.ber2)
#                      ,binary.link="logit"
#                      ,beta.surv = c(intmat$weibull.int[j],beta.surv1,beta.surv2)
#                      ,surv.dist = "weibull"
#                      ,scale=0.5
#                      ,fixed.cens=120
#                      ,trunc.pt=NULL
# )
# 
# out2 <- runmodels_nocr(dat)


##################################################################################
# 
# # saving the current RNG state for later use
# set.seed(123)
# tmp <- .Random.seed
# .Random.seed <- tmp
# 
# 
# ######## no parallelization version #####
# #### Testing with a few simulations per scenario and no 
# 
# nsim <- 3  # number of sim iterations per simulation setting
# out.array <- array( NA, dim=c( length(Ns),length(UORs),length(CPs),length(Sdists),dim(out2),nsim ),
#                     dimnames=list( samp_size= paste0("N",Ns),
#                                    unobs_rate=paste0("UOR",UORs),
#                                    cens_prop=paste0("CP",CPs),
#                                    surv_dist=paste0("SD",Sdists),
#                                    pred_quant=dimnames(out2)[[1]],
#                                    point_ci=dimnames(out2)[[2]],
#                                    models=dimnames(out2)[[3]],
#                                    simiter=paste0("sim",1:nsim)
#                     )
# ) 
# 
# 
# 
# for( i in 1:nrow(scen.mat) ) # for each simulation scenario
# {
#   base::cat(paste0("scen",i," "))
#   
#   # set the simulation scenario settings
#   nsamp <- scen.mat$n[i]
#   uorate <- scen.mat$unobs.rate[i]
#   cprop <- scen.mat$censor.prop[i]
#   sdist <- as.character( scen.mat$surv.dist[i] )
#   tpoint <- NULL
#   if(cprop==0) {tpoint <- 120} 
#   
#   x1 <- c(rep(1,nsamp/2),rep(0,nsamp/2))
#   x2 <- rep(0:1,nsamp/2)
#   design <- cbind(1,x1,x2)
#   
#   crit <- intmat$unobs.rate==uorate & intmat$censor.prop==cprop 
#   berint <- intmat$ber.int[crit]
#   survint <- ifelse(sdist=="lognormal",intmat$lognorm.int[crit],intmat$weibull.int[crit])
#   
#   # set output array indices based on simulation settings
#   a1 <- ifelse(nsamp==40,"N40","N90")
#   a2 <- ifelse(uorate==0.15,"UOR0.15","UOR0.4")
#   a3 <- ifelse(cprop==0,"CP0",ifelse(cprop==1/3,"CP0.333333333333333",ifelse(cprop==2/3,"CP0.666666666666667","CP1")))
#   a4 <- ifelse(sdist=="lognormal","SDlognormal","SDweibull")
#   
#   
#   # run the simulation nsim times
#   out <- sapply(integer(nsim), function(i) { runsim(crmodels=F
#                                                     ,design=design 
#                                                     ,beta.ber = c(berint,beta.ber1,beta.ber2)
#                                                     ,binary.link="logit"
#                                                     ,beta.surv = c(survint,beta.surv1,beta.surv2)
#                                                     ,surv.dist = sdist
#                                                     ,scale=0.5
#                                                     ,fixed.cens=120
#                                                     ,trunc.pt=tpoint) } 
#                 , simplify="array")
#   
#   
#   # store the output 
#   out.array[a1,a2,a3,a4, , , , ] <- out
# }




#########################################################################################################

# functionized non-parallel version
# dumdum <- run_allscen(  Ns=Ns,UORs=UORs,CPs=CPs,Sdists=Sdists,intmat=intmat,scen.mat=scen.mat,nsim=3 )
# dumdum["N40","UOR0.15","CP1","SDlognormal",,,,"sim1"]
# 




### Attributes for the array
# $dim
# [1]  2  2  4  2 20  3  4  3
# 
# $dimnames
# $dimnames$samp_size
# [1] "N40" "N90"
# 
# $dimnames$unobs_rate
# [1] "UOR0.15" "UOR0.4" 
# 
# $dimnames$cens_prop
# [1] "CP0"                 "CP0.333333333333333" "CP0.666666666666667" "CP1"                
# 
# $dimnames$surv_dist
# [1] "SDlognormal" "SDweibull"  
# 
# $dimnames$pred_quant
# [1] "b1"       "b2"       "t50_00"   "t50_10"   "t50_01"   "t50_11"   "p30_00"   "p30_10"   "p30_01"   "p30_11"  
# [11] "p60_00"   "p60_10"   "p60_01"   "p60_11"   "curep_00" "curep_10" "curep_01" "curep_11" "cureb1"   "cureb2"  
# 
# $dimnames$point_ci
# [1] "est"   "lwr95" "upr95"
# 
# $dimnames$models
# [1] "CoxPH" "ParLN" "ParWB" "ParGG"
# 
# $dimnames$simiter
# [1] "sim1" "sim2" "sim3"

########################################################################################

### We're going to "parallelize" by running two separate instances of Rstudio
### the first will run the weibull scenarios, the second will run the lognormal scenarios

set.seed(123)
dumdum_weib <- run_allscen(Ns=Ns,UORs=UORs,CPs=CPs,Sdists="weibull",intmat=intmat,scen.mat=scen.mat,nsim=100,crmodels=T )
save(dumdum_weib,file=paste0(root,data,"Weibull Simulation.RData"))




set.seed(124)
dumdum_logn <- run_allscen(Ns=Ns,UORs=UORs,CPs=CPs,Sdists="lognormal",intmat=intmat,scen.mat=scen.mat,nsim=100,crmodels=T )
save(dumdum_logn,file=paste0(root,data,"Lognormal Simulation.RData"))







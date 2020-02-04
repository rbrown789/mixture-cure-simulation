
########################################################################
# This script generates all the results plots from Simulation 1.  Apologies
# to anyone trying to figure out what's going on, cause I didn't anything.
#
########################################################################


library(survival)
library(rms)
library(flexsurv)
library(smcure)
library(abind)
library(R.utils)
library(truncdist)
library(lattice)
library(latticeExtra)


root <- "C:/Users/rbrow/OneDrive/Documents/Public Github/mixture-cure-simulation/"
data <- "results/"
code <- "code/"
graphics <- "Graphics/"


source(paste0(root,code,"sim_results_functions.R"))


resfin.df <- todf(resfin)


##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

## Generate all bias plots
qnms <- dimnames(resfin)$pred_quant[1:18]

for( i in qnms)
{
  parres <- todf.par(i,resfin)
  
  # true <-  sprintf("%.2f", round(parres$true[1],2))
  
  
  if(i %in% c("b1","b2")) { parres <- parres[!parres$model %in% c("CoxPH","MCPH") ,]}
  if(i %in% c("curep_00","curep_10","curep_01","curep_11")) { parres <- parres[parres$model %in% c("MCAFT","MCPH") ,]}
  
  parres40 <- parres[parres$n=="N40",]
  parres90 <- parres[parres$n=="N90",]
  
  ylim <- c(min(parres$bias,na.rm=T),max(parres$bias,na.rm=T))
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Plots/Bias/",i,"_bias.pdf"),height=8,width=8)}
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="bias",bynm="cp",modelnm="model",
                dat=parres40,ylim=ylim,maintit=paste0("N=40 - ",i),
                bigxlab="Censoring Proportion",bigylab="Bias",h=0 )
  
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="bias",bynm="cp",modelnm="model",
                dat=parres90,ylim=ylim,maintit=paste0("N=90 - ",i),
                bigxlab="Censoring Proportion",bigylab="Bias",h=0 )
  if(export) {dev.off()}
}


## Generate all relative bias plots
for(i in qnms)
{
  parres <- todf.par(i,resfin)
  
  # true <-  sprintf("%.2f", round(parres$true[1],2))
  
  
  if(i %in% c("b1","b2")) { parres <- parres[!parres$model %in% c("CoxPH","MCPH") ,]}
  if(i %in% c("curep_00","curep_10","curep_01","curep_11")) { parres <- parres[parres$model %in% c("MCAFT","MCPH") ,]}
  
  parres40 <- parres[parres$n=="N40",]
  parres90 <- parres[parres$n=="N90",]
  
  ylim <- c(min(parres$relbias,na.rm=T),max(parres$relbias,na.rm=T))
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Plots/Relative Bias/",i,"_relbias.pdf"),height=8,width=8)}
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="relbias",bynm="cp",modelnm="model",
                dat=parres40,ylim=ylim,maintit=paste0("N=40 - ",i),
                bigxlab="Censoring Proportion",bigylab="Relative Bias",h=0 )
  
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="relbias",bynm="cp",modelnm="model",
                dat=parres90,ylim=ylim,maintit=paste0("N=90 - ",i),
                bigxlab="Censoring Proportion",bigylab="Relative Bias",h=0 )
  if(export) {dev.off()}
}



## Generate all MSE plots
qnms <- dimnames(resfin)$pred_quant[1:18]

for( i in qnms)
{
  parres <- todf.par(i,resfin)
  
  # true <-  sprintf("%.2f", round(parres$true[1],2))
  
  
  if(i %in% c("b1","b2")) { parres <- parres[!parres$model %in% c("CoxPH","MCPH") ,]}
  if(i %in% c("curep_00","curep_10","curep_01","curep_11")) { parres <- parres[parres$model %in% c("MCAFT","MCPH") ,]}
  
  parres40 <- parres[parres$n=="N40",]
  parres90 <- parres[parres$n=="N90",]
  
  ylim <- c(min(parres$mse,na.rm=T),max(parres$mse,na.rm=T))
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Plots/MSE/",i,"_mse.pdf"),height=8,width=8)}
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="mse",bynm="cp",modelnm="model",
                dat=parres40,ylim=ylim,maintit=paste0("N=40 - ",i),bigxlab="Censoring Proportion",bigylab="MSE" )
  
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="mse",bynm="cp",modelnm="model",
                dat=parres90,ylim=ylim,maintit=paste0("N=90 - ",i),bigxlab="Censoring Proportion",bigylab="MSE" )
  if(export) {dev.off()}
}



## Generate all Variance plots
qnms <- dimnames(resfin)$pred_quant[1:18]

for( i in qnms)
{
  parres <- todf.par(i,resfin)
  
  # true <-  sprintf("%.2f", round(parres$true[1],2))
  
  
  if(i %in% c("b1","b2")) { parres <- parres[!parres$model %in% c("CoxPH","MCPH") ,]}
  if(i %in% c("curep_00","curep_10","curep_01","curep_11")) { parres <- parres[parres$model %in% c("MCAFT","MCPH") ,]}
  
  parres40 <- parres[parres$n=="N40",]
  parres90 <- parres[parres$n=="N90",]
  
  ylim <- c(min(parres$var,na.rm=T),max(parres$var,na.rm=T))
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Plots/Variance/",i,"_var.pdf"),height=8,width=8)}
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="var",bynm="cp",modelnm="model",
                dat=parres40,ylim=ylim,maintit=paste0("N=40 - ",i),bigxlab="Censoring Proportion",bigylab="Variance" )
  
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="var",bynm="cp",modelnm="model",
                dat=parres90,ylim=ylim,maintit=paste0("N=90 - ",i),bigxlab="Censoring Proportion",bigylab="Variance" )
  if(export) {dev.off()}
}



## Generate all coverage rate plots
qnms <- dimnames(resfin)$pred_quant[1:18]

for( i in qnms)
{
  parres <- todf.par(i,resfin)
  
  if(i %in% c("b1","b2")) { parres <- parres[!parres$model %in% c("CoxPH","MCPH") ,]}
  if(i %in% c("curep_00","curep_10","curep_01","curep_11")) { parres <- parres[parres$model %in% c("MCAFT","MCPH") ,]}
  
  parres40 <- parres[parres$n=="N40",]
  parres90 <- parres[parres$n=="N90",]
  
  ylim <- c(min(parres$covrate,na.rm=T),max(parres$covrate,na.rm=T))
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Plots/Coverage Rate/",i,"_covrate.pdf"),height=8,width=8)}
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="covrate",bynm="cp",modelnm="model",
                dat=parres40,ylim=ylim,maintit=paste0("N=40 - ",i),bigxlab="Censoring Proportion",
                bigylab="Coverage Rate",h=0.95 )
  
  trellis.plot( xvarnm="sdist", yvarnm="uor", measnm="covrate",bynm="cp",modelnm="model",
                dat=parres90,ylim=ylim,maintit=paste0("N=90 - ",i),bigxlab="Censoring Proportion",
                bigylab="Coverage Rate",h=0.95 )
  if(export) {dev.off()}
}

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

###### Making Pareto Charts and Pareto Fronts ########


#################################################################################

# Make pareto charts for RMSE and Bias
alldat <- NULL
qnms <- dimnames(resfin)$pred_quant[1:14]

for( i in qnms)
{
  parres <- todf.par(i,resfin)
  
  if(i %in% c("b1","b2")) { parres <- parres[!parres$model %in% c("CoxPH","MCPH") ,]}
  
  #### bias ######
  
  # fit linear model to the simulation results
  mod <- lm(bias ~ model + n + uor + cp + sdist + model*n +  model*uor + model*cp + model*sdist ,data=parres)
  pprepdat <- paretoprep(mod,eflab)
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Pareto Charts/Bias/pareto_",i,"_bias.pdf"),height=7,width=7)}
  par(mar=c(5.1,12.1,4.1,2.1))
  barplot( pprepdat$logworth, horiz=T,main=paste0("Pareto Chart: ",i," Bias"),xlab="Log Worth")
  axis(2,at=pprepdat$loc*1.2-0.5,labels=pprepdat$lab,las=1,col.ticks="white")
  box( bty = "L") 
  if(export) {dev.off()}
  
  pprepdat$meas <- "bias"
  pprepdat$par <- i
  
  alldat <- rbind(alldat,pprepdat)
  
  #### RMSE ######
  parres$rmse <- sqrt(parres$mse)
  mod <- lm(rmse ~ model + n + uor + cp + sdist + model*n +  model*uor + model*cp + model*sdist ,data=parres)
  pprepdat <- paretoprep(mod,eflab)
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Pareto Charts/RMSE/pareto_",i,"_rmse.pdf"),height=7,width=7)}
  par(mar=c(5.1,12.1,4.1,2.1))
  barplot( pprepdat$logworth, horiz=T,main=paste0("Pareto Chart: ",i," RMSE"),xlab="Log Worth")
  axis(2,at=pprepdat$loc*1.2-0.5,labels=pprepdat$lab,las=1,col.ticks="white")
  box( bty = "L") 
  if(export) {dev.off()}
  
  pprepdat$meas <- "rmse"
  pprepdat$par <- i
  
  alldat <- rbind(alldat,pprepdat)
  
  # standard error
  parres$se <- sqrt(parres$var)
  mod <- lm(se ~ model + n + uor + cp + sdist + model*n +  model*uor + model*cp + model*sdist ,data=parres)
  pprepdat <- paretoprep(mod,eflab)
  
  
  if(export) { pdf(paste0(root,graphics,"Simulation Pareto Charts/StdErr/pareto_",i,"_stderr.pdf"),height=7,width=7)}
  par(mar=c(5.1,12.1,4.1,2.1))
  barplot( pprepdat$logworth, horiz=T,main=paste0("Pareto Chart: ",i," stderr"),xlab="Log Worth")
  axis(2,at=pprepdat$loc*1.2-0.5,labels=pprepdat$lab,las=1,col.ticks="white")
  box( bty = "L") 
  if(export) {dev.off()}
  
  pprepdat$meas <- "se"
  pprepdat$par <- i
  
  alldat <- rbind(alldat,pprepdat)
}




### Make 2D Pareto Chart (Frontier?) ####
tmp <- alldat[,c("logworth","lab","meas","par")]
tmpbias <- tmp[tmp$meas=="bias",]
names(tmpbias) <- c("lw_bias","lab","meas","par")
tmprmse <- tmp[tmp$meas=="rmse",]
names(tmprmse) <- c("lw_rmse","lab","meas","par")
tmpse <- tmp[tmp$meas=="se",]
names(tmpse) <- c("lw_se","lab","meas","par")

par2dat <- merge(tmpbias,tmprmse,by=c("lab","par"),all=T)
par2dat <- merge(par2dat,tmpse,by=c("lab","par"),all=T)

labs <- c("Model","Censoring Proportion","Unobserved Rate","Sample Size","Survival Distribution",
          "Model*Censoring Proportion","Model*Unobserved Rate","Model*Sample Size","Model*Survival Distribution")
pchs <- c(17,17,17,17,17,16,16,16,16)
cols <- c("black","red","steelblue","orange","purple","red","steelblue","orange","purple")



if(export) { pdf(paste0(root,graphics,"Simulation Pareto Charts/2d Pareto_biasbystderr.pdf"),height=9,width=9)}
plot("n",ylim=c(0,max(par2dat$lw_se)),xlim=c(0,max(par2dat$lw_bias)),axes=F,
     xlab="Log Worth(Bias)",ylab="Log Worth (Std Err)")
for( i in 1:length(labs))
{
  tmp <- par2dat[par2dat$lab==labs[i],]
  points(tmp$lw_bias,tmp$lw_se,pch=pchs[i],col=cols[i],cex=1.5)
}

axis(1);axis(2);box()

legend("topright",labs,pch=pchs,col=cols,bty="n")
if(export) {dev.off() }



if(export) { pdf(paste0(root,graphics,"Simulation Pareto Charts/2d Pareto_biasbyrmse.pdf"),height=9,width=9)}
plot("n",ylim=c(0,max(par2dat$lw_rmse)),xlim=c(0,max(par2dat$lw_bias)),axes=F,
     xlab="Log Worth(Bias)",ylab="Log Worth (RMSE)")
for( i in 1:length(labs))
{
  tmp <- par2dat[par2dat$lab==labs[i],]
  points(tmp$lw_bias,tmp$lw_rmse,pch=pchs[i],col=cols[i],cex=1.5)
}

axis(1);axis(2);box()

legend("topright",labs,pch=pchs,col=cols,bty="n")
if(export) {dev.off() }


####################################################################################




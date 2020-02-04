
#############################################################
# This script grabs the simulation settings, as well as simulation results,
# calculates the true values for all simulation settings, formats them properly to
# match the simluation results output (multidimensional array) and then calcualtes
# bias, variance, MSE, and coverage rate for all models and all simulation setting.
# 
# It also defines some data aggregation and plotting functions for use in 
# the simresults_plot.R script.
#
#
##############################################################


# library(survival)
# library(rms)
# library(flexsurv)
# library(smcure)
# library(abind)
# library(R.utils)
# library(truncdist)
# library(lattice)
# library(latticeExtra)
# 
# 
# root <- "U:/Censored Data Project/"
# data <- "Data Sets/"
# graphics <- "Graphics/"
# code <- "Code/"

# load the simulation settings and results
load(paste0(root,data,"Lognormal Simulation.Rdata"))
load(paste0(root,data,"Weibull Simulation.Rdata"))
load(paste0(root,data,"simsetting.Rdata"))
root <- "C:/Users/rbrow/OneDrive/Documents/Public Github/mixture-cure-simulation/"
data <- "results/"
code <- "code/"
graphics <- "Graphics/"

simres <- abind(dumdum_logn,dumdum_weib,along=4)
names(dimnames(simres)) <- names(dimnames(dumdum_logn))

# Calculate the true values for the quantities of interest
intmat2 <- rbind(intmat[,c("unobs.rate","censor.prop","ber.int")],
                 intmat[,c("unobs.rate","censor.prop","ber.int")])

intmat2$surv.int <- c(intmat$weibull.int,intmat$lognorm.int)
intmat2$surv.dist <- c(rep(1,nrow(intmat2)/2),rep(2,nrow(intmat2)/2)) # 1 for weibull, 2 for lognormal

 
truevalcalc <- function(imvec,beta.ber1,beta.ber2,beta.surv1,beta.surv2)
{
  ############################################################################
  # Function to generate the true values for the quantities of interest,
  # given the simulation parameters
  #
  ############################################################################
  
  
  # imvec <-  as.numeric( intmat2[9,] )
  
  beta.ber <- c( imvec[3],beta.ber1,beta.ber2)
  beta.surv <- c( imvec[4],beta.surv1,beta.surv2)
  dist <- ifelse(imvec[5]==1,"weibull","lnorm")
  
  trun <- Inf
  if(imvec[2]==0){ trun <- 120 }
  
  
  nd <- as.matrix( data.frame(x0=1,x1=c(0,1,0,1),x2=c(0,0,1,1)) )
  
  # calculate true cure probability 
  pcure <- as.numeric( expit( nd %*% beta.ber ) )
  names(pcure) <- c("curep_00","curep_10","curep_01","curep_11")
  
  # calculate median survival time and survival probability at 30, 60 seconds in uncured group
  logmean <- as.numeric( nd %*% beta.surv)
  if(dist=="lnorm") { 
      t50u <- qtrunc(0.5,spec=dist,a=0,b=trun,meanlog=logmean,sdlog=scale)  # median survival time
      p30u <- 1-ptrunc(30,spec=dist,a=0,b=trun,meanlog=logmean,sdlog=scale) # survival at 30s
      p60u <- 1-ptrunc(60,spec=dist,a=0,b=trun,meanlog=logmean,sdlog=scale) # survival at 60s
  }
  
  if(dist=="weibull") { 
    t50u <- qtrunc(0.5,spec=dist,a=0,b=trun,shape=1/scale,scale=exp(logmean))  # median survival time
    p30u <- 1-ptrunc(30,spec=dist,a=0,b=trun,shape=1/scale,scale=exp(logmean)) # survival at 30s
    p60u <- 1-ptrunc(60,spec=dist,a=0,b=trun,shape=1/scale,scale=exp(logmean)) # survival at 60s
  }
  
  # calculate equivalent values in the full dataset
  # note equivalent formulations of cure rate model: S = pi*Su + 1-pi ;  F = pi*Fu ; Q = Qu/pi
  # where S is the survival function, F is the CDF, and Q is the quantile function overall
  # and Su,Fu, and Qu are the counterpart functions the in only the uncured population
  # and pi is the uncured probability
  t50 <- t50u/(1-pcure)
  p30 <- p30u*(1-pcure)+pcure
  p60 <- p60u*(1-pcure)+pcure
  
  names(t50) <- c("t50_00","t50_10","t50_01","t50_11")
  names(p30) <- c("p30_00","p30_10","p30_01","p30_11")
  names(p60) <- c("p60_00","p60_10","p60_01","p60_11")
  
  bber <- beta.ber[2:3];names(bber) <- c("cureb1","cureb2")
  bsurv <- beta.surv[2:3];names(bsurv) <- c("b1","b2")
  
  # put in same order with same names as outcomes in simulation results
  out <- c(bsurv,t50,p30,p60,pcure,bber)
  return(out)
}

truevals <- t( apply(intmat2,1,truevalcalc,beta.ber1,beta.ber2,beta.surv1,beta.surv2) )
intmat3 <- cbind(intmat2,truevals)

#############################################################

# Reformat the true values so they are in a multidimensional array with same format as simulation results
quantnms <- c("b1","b2","t50_00","t50_10","t50_01","t50_11","p30_00","p30_10","p30_01","p30_11",
              "p60_00","p60_10","p60_01","p60_11","curep_00","curep_10","curep_01","curep_11","cureb1","cureb2")

true.array <- array( NA, 
                     dim=c(length(Ns),length(UORs),length(CPs),length(Sdists),length(quantnms) ),
                     
                     dimnames=list( samp_size= paste0("N",Ns),
                                    unobs_rate=paste0("UOR",UORs),
                                    cens_prop=paste0("CP",CPs),
                                    surv_dist=paste0("SD",Sdists),
                                    pred_quant=quantnms
                     )
) 

intmat4 <- rbind(intmat3,intmat3)
intmat4$n <- c( rep(40,nrow(intmat3)),rep(90,nrow(intmat3)))

for( i in 1:nrow(intmat4))
{
  # i <- 1
  
  nms <- c("b1","b2","t50_00","t50_10","t50_01","t50_11",
           "p30_00","p30_10","p30_01","p30_11", "p60_00","p60_10","p60_01","p60_11",
           "curep_00","curep_10","curep_01","curep_11","cureb1","cureb2")
  truevals <- as.numeric ( intmat4[i,nms] )
  names(truevals) <- nms 
  
  nsamp <- intmat4$n[i]
  uorate <- intmat4$unobs.rate[i]
  cprop <- intmat4$censor.prop[i]
  sdist <- intmat4$surv.dist[i]

  # set output array indices based on simulation settings
  a1 <- ifelse(nsamp==40,"N40","N90")
  a2 <- ifelse(uorate==0.15,"UOR0.15","UOR0.4")
  a3 <- ifelse(cprop==0,"CP0",ifelse(cprop==1/3,"CP0.333333333333333",ifelse(cprop==2/3,"CP0.666666666666667","CP1")))
  a4 <- ifelse(sdist==2,"SDlognormal","SDweibull")
  
  true.array[a1,a2,a3,a4,] <- truevals
}  


#########################################################################################################

# Calculate Bias, Variance, MSE, Coverage Rate for all models and all scenarios
dum <- expand.grid(N=paste0("N",Ns),UOR=paste0("UOR",UORs),CP=paste0("CP",CPs),Sdist=paste0("SD",Sdists),stringsAsFactors=F)
resfin <- array( NA, 
                 dim=c(length(Ns),length(UORs),length(CPs),length(Sdists),length(quantnms),6,6 ),
                 
                 dimnames=list( samp_size= paste0("N",Ns),
                                unobs_rate=paste0("UOR",UORs),
                                cens_prop=paste0("CP",CPs),
                                surv_dist=paste0("SD",Sdists),
                                pred_quant=quantnms,
                                models=c("CoxPH" ,"ParLN" ,"ParWB", "ParGG" ,"MCAFT", "MCPH"),
                                meas=c("true","mean","bias","var","mse","covrate")
                 )
) 



for(i in 1:nrow(dum))
{  
  # i <- 4
  a1 <- dum$N[i]
  a2 <- dum$UOR[i]
  a3 <- dum$CP[i]
  a4 <- dum$Sdist[i]
  
  # subset arrays to the current simulation setting
  tmptrue <- true.array[a1,a2,a3,a4,]
  tmpsim <- simres[a1,a2,a3,a4,,,,]
  
  # calculate mean and variance of point estimate over the simulations
  tmpsimmean <- colMeans(aperm(tmpsim,c(4,1,2,3)),na.rm=T)  
  tmpsimmean <- tmpsimmean[,"est",]
  tmpsimvar <- apply(tmpsim,c(1,2,3),var,na.rm=T)
  tmpsimvar <- tmpsimvar[,"est",]
  tmptrue <- matrix(tmptrue,nrow=nrow(tmpsimmean),ncol=ncol(tmpsimmean))
  attributes(tmptrue) <- attributes(tmpsimmean)
  
  
  # calculate bias
  tmpsimbias <- tmpsimmean-tmptrue
  
  # calculate MSE
  tmpsimmse <- tmpsimbias^2 + tmpsimvar
  
  #### calculate coverage rate ####
  
  # merge true value into CI matrix
  tmpsimci <- tmpsim[,c("lwr95","upr95"),,]
  tmptrue2 <-  array(tmptrue,dim=c(20,1,6,100))
  tmpsimci <- abind(tmpsimci,tmptrue2,along=2)
  dimnames(tmpsimci)[[2]] <- c("lwr95","upr95","true") 
  
  # coverage indicator
  cover <- 1*( tmpsimci[,"true",,] > tmpsimci[,"lwr95",,] & tmpsimci[,"true",,] < tmpsimci[,"upr95",,] )
  
  # coverage rate
  tmpsimcovrate <- colMeans(aperm(cover,c(3,1,2)),na.rm=T)
  
  # bind them for the output
  tmpout <- abind(true=tmptrue,mean=tmpsimmean,bias=tmpsimbias,var=tmpsimvar,mse=tmpsimmse,covrate=tmpsimcovrate,along=3)
  
  # bind to the full output dataset
  resfin[a1,a2,a3,a4,,,] <- tmpout
}


dimnames(resfin)$cens_prop <- c("CP0","CP1/3","CP2/3","CP1")
nms <- names(dimnames(resfin))


#######################################################################################################################

relbias <- resfin[,,,,,,"bias"]/resfin[,,,,,,"true"]
resfin <- abind(resfin,relbias,along=7)
dimnames(resfin)[[7]][7] <- "relbias"
names(dimnames(resfin)) <- nms

resfin["N90","UOR0.15","CP0","SDweibull","t50_00", , ] 





# save(resfin,file=paste0(root,data,"Simulation Results.Rdata") )
# 
# load(file=paste0(root,data,"Simulation Results.Rdata") )

todf <- function(resarray)
{
  ##############################################################
  # convert all information on all quantitites into a dataframe
  # with columns for simulation settings 
  #
  #################################################################
  
  out <- NULL
  for(h in dimnames(resarray)[[5]])
  {
    for(i in dimnames(resarray)[[1]] )
    {
      for(j in dimnames(resarray)[[2]] )
      {
        for(k in dimnames(resarray)[[3]] )
        {
          for(l in dimnames(resarray)[[4]] )
          {
            tmp <- data.frame( resarray[i,j,k,l,h,,] )
            tmp$model <- row.names(tmp)
            tmp$par <- h
            tmp$n <- i
            tmp$uor <- j
            tmp$cp <- k
            tmp$sdist <- l
            
            out <- rbind(out,tmp)
          }
        }
      }
    }
  }
  return(out)
}



todf.par <- function(parname,resarray)
{
  ##############################################################
  # convert all information on a certain parameter into a dataframe
  # with columns for simulation settings 
  #
  #################################################################
  
  out <- NULL
  
  for(i in dimnames(resarray)[[1]] )
  {
    for(j in dimnames(resarray)[[2]] )
    {
      for(k in dimnames(resarray)[[3]] )
      {
        for(l in dimnames(resarray)[[4]] )
        {
          tmp <- data.frame( resarray[i,j,k,l,parname,,] )
          tmp$model <- row.names(tmp)
          tmp$n <- i
          tmp$uor <- j
          tmp$cp <- k
          tmp$sdist <- l
            
          out <- rbind(out,tmp)
        }
      }
    }
  }
  
  return(out)
}




plotfunc <- function(measnm,bynm,modelnm,dat,ax=NULL,moddf, h = NULL, ...)
{
  meas <- dat[[measnm]]
  by <- dat[[bynm]]
  model <- dat[[modelnm]]
  
  mods <- sort(unique(model))
  
  lvls <- unique(by)
  xlocs <- 1:length(lvls)
  xl <- data.frame(lvls=lvls,xlocs=xlocs)
  
  
  pdat <- data.frame(meas=meas,by=by,model=model)
  
  plot("n",xlim=c(1,max(xlocs)),axes=F, ... )
  box()
  
  
  if(!is.null(ax)) {
    axis(ax[1],at=xlocs,labels=lvls)
    axis(ax[2])
  }
    
  for( i in 1:length(mods) )
  {
    datm <- pdat[dat$model==mods[i],]
    datm <- merge(datm,xl,by.x="by",by.y="lvls")
    datm <- datm[order(datm$xlocs),]
    
    col <- moddf$modcols[moddf$model==mods[i]]
    lty <- moddf$modlty[moddf$model==mods[i]]
    pch <- moddf$modpch[moddf$model==mods[i]]
    
    lines(datm$xlocs,datm$meas,col=col,lty=lty)
    points(datm$xlocs,datm$meas,col=col,pch=pch)
  }
  
  if(!is.null(h)) { abline(h=h,lty=3)}
}


trellis.plot <- function( xvarnm,yvarnm,measnm,bynm,modelnm,dat, ylim,maintit="main",bigxlab="",bigylab="",h=NULL )
{
  # dat <- zz
  # 
  # xvarnm <- "sdist"
  # yvarnm <- "uor"
  # measnm <- "bias"
  # bynm <- "cp"
  # modelnm <- "model"
  
  # set colors, point type and line type for each model
  mods <- c("CoxPH","ParLN","ParWB","ParGG","MCAFT","MCPH")
  modcols <- c("black","steelblue","orange","red","grey80","purple")
  modlty <- c(1,1,1,1,1,1)
  modpch <- c(0,1,2,4,5,6)
  moddf <- data.frame(model=mods,modcols=modcols,modlty=modlty,modpch=modpch,stringsAsFactors = F)
  
  
  
  mods <- sort(unique(dat$model))
  nmods <- length(mods)
  
  
  xun <- unique(dat[[xvarnm]])
  yun <- unique(dat[[yvarnm]])
  
  par(omi=rep(1.0, 4), mar=c(0,0,0,0), mfrow=c(length(xun),length(yun)))
  
  for( i in 1:length(yun) )
  {
    for( j in 1:length(xun) )
    {
      yvi <- yun[i]
      xvj <- xun[j]
      
      tmp <- dat[ dat[[xvarnm]]==xvj & dat[[yvarnm]]==yvi, ]
      
      if(i==1 & j==length(xun)) { plotfunc( measnm,bynm,modelnm,tmp,ax=c(3,4),ylim=ylim,moddf=moddf,h=h )
      } else if(i==length(yun) & j==1 ) { plotfunc( measnm,bynm,modelnm,tmp,ax=c(1,2), ylim=ylim,moddf=moddf,h=h  )
      } else { plotfunc( measnm,bynm,modelnm,tmp, ylim=ylim,moddf=moddf,h=h ) }
      
      
      if(i==1) { mtext(xun[j],side=3,line=2.5)  }
      if(j==length(xun)) { mtext(yun[i],side=4,line=2.5)  }
      
      
    }
  }
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 3), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  
  mtext(bigxlab,side=1,line=0,cex=1.2)
  mtext(bigylab,side=2,line=0,cex=1.2)
  mtext(maintit,side=3,line=1,cex=1.5)
  
  legend(x=0.95,y=-0.75, moddf$model, xpd = TRUE, horiz = F, inset = c(0,0),
         bty = "n", pch = moddf$modpch, col = moddf$modcols,lty=moddf$modlty, cex = 1 )
  
}


########## some helper and plotting functions for pareto charts #########

varlab <- function(varstring,key)
{ 
  ########################################################################
  # Function for replacing variable names with labels.
  # Varstring is the string of variable names, key is a dataframe
  # that matches the string (column 1) to the label (column 2)
  ########################################################################
  out <- rep("",length(varstring))  
  for(i in 1:length(varstring))
  {      
    ind <- varstring[i] %in% key[,1]
    if(ind==FALSE) { out[i] <- varstring[i] } else { out[i] <- key[,2][ key[,1]%in% varstring[i] ] }        
  }  
  return(out)
}


paretoprep <- function(lmfit,labdf)
{
  ########################################################################
  # Function for generating logworth data for plotting in pareto charts/pareto fronts
  #  - lmfit is a fitted model from lm()
  #  - labdf is a dataframe with strings and labels
  ########################################################################
  
  # fit anova table
  an <- anova(lmfit)
  
  # put into dataframe and calculate "logworth" = -log10(pval)
  dat <- data.frame(pvalue=an[[5]])
  dat$effect <- row.names(data.frame(an))
  dat$logworth <- -log(dat$pvalue)
  dat <- dat[!is.na(dat$pvalue),]
  
  # append labels to the data
  dat$lab <- varlab(dat$effect,key=labdf)
  
  # sort by logworth and add the plot location
  dat <- dat[order(dat$logworth),]
  dat$loc <- 1:nrow(dat)
  
  return(dat)
}

eflab <- data.frame(string=c("model","cp","uor","n","sdist","model:cp","model:uor","model:n","model:sdist"),
                    label=c("Model","Censoring Proportion","Unobserved Rate","Sample Size","Survival Distribution",
                            "Model*Censoring Proportion","Model*Unobserved Rate","Model*Sample Size",
                            "Model*Survival Distribution"),stringsAsFactors=F
)

paretochart <- function(pprepdat)
{
  ########################################################################
  # Function to plot pareto chart
  #  - pprepdat is output from paretoprep() function
  #  - labdf is a dataframe with strings and labels
  ########################################################################
  
  par(mar=c(5.1,12.1,4.1,2.1))
  barplot( pprepdat$logworth, horiz=T)
  axis(2,at=pprepdat$loc*1.2-0.5,labels=pprepdat$lab,las=1,col.ticks="white")
  box( bty = "L") 
}




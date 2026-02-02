# Do not run on a login node!!!!
hostName=Sys.info()[4]
if (substr(hostName, 1,12)=="belenoslogin") stop("!!!!!!!!!!!!!!!! You are on a login node !!!!!!!!!!!!!!!!!!!!!")
if (substr(hostName, 1,12)=="taranislogin") stop("!!!!!!!!!!!!!!!! You are on a login node !!!!!!!!!!!!!!!!!!!!!")

system("export OMP_NUM_THREADS=120")

rm(list=ls())
graphics.off()

setwd("./.")

start_time <- Sys.time()
print(paste("start time : ",start_time,sep=""))

library(data.table)
library(foreach)
library(RColorBrewer)
library(plyr) #for rbind.fill()
library(dplyr) #load after plyr
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(scales)
library(rAmCharts)
library(pipeR)
library(randomForestSRC)
library(R.utils)
library(viking)
library(parallel)
suppressMessages(library(doMC))
source("plot_agreg.R")
source("Opera_dev.R")
source("eValues.R")
source("transition_learning.R")
source("transition_learning_xgboost.R")
source("SleepingMarkovHedge.R")
source("SleepingMarkovHedgeCalib.R")
source("SleepingMarkovBOA.R")
source("SleepingMarkovBOACALIB.R")
source("MBOA.R")
source("MH.R")
source("BOA.R")
source("MLpol.R")
source("MMLpol.R")
source("MLprod.R")
source("MMLprod.R")
source("EWA.R")
source("NBM.R")
source("fonction_prevision_FS.R")
source("fonction_prevision_FSCALIB.R")
source("fonction_prevision_BOA.R")
source("fonction_prevision_SOCOBOA.R")
source("fonction_prevision_FSBOA.R")
source("fonction_prevision_FSSOCOBOA.R")
source("fonction_prevision_FSBOACALIB.R")
source("fonction_prevision_FSSOCOBOACALIB.R")
source("fonction_prevision_SMH.R")
source("fonction_prevision_SMHCalib.R")
source("fonction_prevision_SMBOA.R")
source("fonction_prevision_SMBOACALIB.R")
source("fonction_prevision_EWA.R")
source("fonction_prevision_NBM.R")
source("fonction_prevision_EWACALIB.R")
source("fonction_prevision_META.R")
source("fonction_prevision_MLpol.R")
source("fonction_prevision_MMLpol.R")
source("fonction_prevision_MLprod.R")
source("fonction_prevision_MMLprod.R")

#options(error=recover)
#options(error=traceback)
#options(show.error.locations = TRUE)

# Use maximum available processes for R, particularly for rfsrc, though it seems to work without this as well...
options(rf.cores=detectCores(), mc.cores=detectCores())


#to create if they don't already exist
dir_exe="/home/mf/dp/mcad/pfitznerl/SAVE/aggregation/"
dir_utemp="/scratch/work/pfitznerl/aggregation/"
dir_plot=paste(dir_utemp,"plot_aggregation/",sep="")
dir_rmse="rmse/"
dir_weights="weights/"
dir_mixture="mixture/"
dir_alpha="alpha/"
dir_eta="eta/"
dir_obs_pred="obs_pred/"
dir_boxplot="boxplots/"
dir_data_without="data/donnee_sans_pearp"
dir_data_with="data/donnee_avec_pearp"
dir_sorties="sorties/"


#function to evaluate the quantiles of the absolute error with 5% steps
eval_quantiles_abs_error <- function(x,y){
  abs_error=abs(x-y) #absolute error
  quantile_abs_error=quantile(abs_error, probs = seq(0, 1, 1/20),na.rm=T)
  return(quantile_abs_error)
}

#function to make the uniform agregation
prevision_UNIFORM <- function(X,Y,T) {
  # Remplace les NA par la moyenne de leur ligne
  X_clean <- t(apply(X, 1, function(row) {
    row_mean <- mean(row, na.rm = TRUE)  # Calcule la moyenne de la ligne (ignore les NA)
    row[is.na(row)] <- row_mean          # Remplace les NA par la moyenne
    row
  }))
  uniform_pred=apply(X_clean,1,mean)#prediction of the uniform aggregation
  matrix_weights= matrix(1/ncol(X), nrow = nrow(X), ncol = ncol(X))
  return(list(predictions=uniform_pred,weights=matrix_weights,mixture=NULL))
}


#launch of the agregations for every sliding window (parallelized over the windows)
#plot rmse and parameters for every station depending on sliding window
#and returns a list containing the sum of rmse across stations
lance_station<-function(agreg,file_station, liste_window, X, Y, T, date.valid, start.date,end.date){
  sta=substr(strsplit(file_station, "_")[[1]][7],1,8) #insee of the station
  ech=strsplit(strsplit(liste_station[1], "/")[[1]][3],"_")[[1]][4] #lead time processed
  n.window=length(liste_window)
  names.experts<-colnames(X)

  rmse.ech.agreg.sta=c() #vector which will contain the rmse of every sliding window
  weights.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #list which will contain for every sliding window the weight matrix (ncol=exp nligne=T) of the experts
  mixture.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain the mixture objects for every sliding window
  predictions.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL)#vector that will contain the aggregation predictions for each window
  alpha.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the best alpha for each iteration
  eta.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the eta alpha for each iteration
  A1.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the A of PSS for each iteration
  B1.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the B of PSS for each iteration
  C1.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the C of PSS for each iteration
  D1.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the D of PSS for each iteration
  A2.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the A of PSS for each iteration
  B2.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the B of PSS for each iteration
  C2.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the C of PSS for each iteration
  D2.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the D of PSS for each iteration
  NERR.ech.agreg.sta=mylist <- sapply(paste("W",liste_window,sep=""),function(x) NULL) #vector which will contain for every sliding window the number of times the GBRT wakes up the wrong experts
  
  registerDoMC(min(n.window,100)) #set number of processes

  out_paralel=foreach(iw=1:n.window,.inorder = TRUE)%do%{#parallelization over the windows, finally gather outputs in a list out_paralel
    rmse.ech.agreg.sta.window=NULL
    window=liste_window[iw]
    #out: list containing the vector of all predictions, the weights used for these predictions, and the corresponding mixture object
    if (agreg=="UNIFORM"){out=prevision_UNIFORM(X,Y,T)}
    if (agreg=="BOA"){out=prevision_BOA(X,Y,window,T,sta,ech,date.valid)} #date.valid required for the getFirstLeadTime() function
    if (agreg=="EWA"){out=prevision_EWA(X,Y,window,T,sta,ech)}
    if (agreg=="NBM"){out=prevision_NBM(X,Y,window,T,sta,ech)}
    if (agreg=="MLpol"){out=prevision_MLpol(X,Y,window,T,sta,ech)}
    if (agreg=="MLprod"){out=prevision_MLprod(X,Y,window,T,sta,ech)}

    if(agreg=="bestConvex") {
      out=NULL
      obj_oracle=oracle(Y=Y,experts=X,model="convex",loss.type="square",niter=10)
      out$mixture=obj_oracle
      out$weights=matrix(obj_oracle$coefficients,ncol=ncol(X),nrow=T,byrow=T)
      out$pred=obj_oracle$pred
    }
    if(agreg=="bestShift") {
      out=NULL
      obj_oracle=oracle(Y=Y,experts=X,model="shifting")#list(model="shifting",m=36))
      rmse.ech.agreg.sta.window=sqrt(1/T*sum(obj_oracle$rmse))
      out$mixture=obj_oracle
      out$pred=rep(NA,T)
      out$weights=matrix(NA,ncol=ncol(X),nrow=T)
      stop("pred is not available for bestShifting, so there will be nans for bestshift; to still get the rmse, one could 
                tweak one of the aggregations and make it predict like the best expert by looking at min(lexp), 
                one can even set weights to zero and 1 if one wants to see the weights...")
    }
    if(agreg=="bestExpert") {
      out=NULL
      obj_oracle=oracle(Y=Y,experts=X,model="expert")
      out$mixture=obj_oracle
      out$weights=matrix(obj_oracle$coefficients,ncol=ncol(X),nrow=T,byrow=T)
      out$pred=obj_oracle$pred
    }
    if(is.null(rmse.ech.agreg.sta.window)) rmse.ech.agreg.sta.window=sqrt(mean(loss(out$pred, Y, loss.type = "square"),na.rm=T))
    return(list(rmse.ech.agreg.sta.window=rmse.ech.agreg.sta.window,weights.ech.agreg.sta.window=out$weights,mixture=out$mixture,
                predictions=out$pred,liste_eta=out$liste_eta,liste_alpha=out$liste_alpha,
                liste_A1=out$A1,liste_B1=out$B1,liste_C1=out$C1,liste_D1=out$D1,
                liste_A2=out$A2,liste_B2=out$B2,liste_C2=out$C2,liste_D2=out$D2,
                liste_NERR=out$NERR))
  }

  if(agreg=="META") names.experts<-out_paralel[[1]]$mixture$names.experts #new experts names

  #write the outputs in a .txt
  if (flag.write.txt) {
    for (iw in 1:n.window){ #for each window
      df_write=cbind(Y,round(out_paralel[[iw]]$predictions,5),out_paralel[[iw]]$mixture$experts,round(out_paralel[[iw]]$weights,5))

      colnames(df_write)=c("Y",agreg,names.experts,paste("weights_",names.experts,sep=""))
      write.table(df_write,paste(dir_exe,dir_sorties,agreg,"_",sta,"_",
                  liste_window[iw],"_",ech,".txt",sep=""), sep=";",col.names = TRUE, row.names = FALSE)
  }
}

  # Store the data associated with the windows in vectors and lists
  for (iw in 1:n.window){ # for each window
    window=liste_window[iw]
    rmse.ech.agreg.sta <- c(rmse.ech.agreg.sta,out_paralel[[iw]]$rmse.ech.agreg.sta.window)
    weights.ech.agreg.sta[[iw]]=out_paralel[[iw]]$weights
    mixture.ech.agreg.sta[[iw]]=out_paralel[[iw]]$mixture
    predictions.ech.agreg.sta[[iw]]=out_paralel[[iw]]$predictions
    if(!is.null(out_paralel[[iw]]$liste_A1)) A1.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_A1 else A1.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_B1)) B1.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_B1 else B1.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_C1)) C1.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_C1 else C1.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_D1)) D1.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_D1 else D1.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_A2)) A2.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_A2 else A2.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_B2)) B2.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_B2 else B2.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_C2)) C2.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_C2 else C2.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_D2)) D2.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_D2 else D2.ech.agreg.sta[[iw]]=NA
    if(!is.null(out_paralel[[iw]]$liste_NERR)) NERR.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_NERR else NERR.ech.agreg.sta[[iw]]=NA
    if (flag.rm.param.window) { # if true, we remove the 60 first iterations of the parameters, because the aggregation isn't stabilized at the beginning
      out_paralel[[iw]]$liste_eta=rm.param.window(out_paralel[[iw]]$liste_eta,window)
      out_paralel[[iw]]$liste_alpha=rm.param.window(out_paralel[[iw]]$liste_alpha,window)
    }
    if(!is.null(out_paralel[[iw]]$liste_eta)) eta.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_eta
    if(!is.null(out_paralel[[iw]]$liste_alpha)) alpha.ech.agreg.sta[[iw]]=out_paralel[[iw]]$liste_alpha
  }

  if (flag.obs.pred){
    for (iw in 1:n.window){
      pred.agreg=predictions.ech.agreg.sta[[iw]]
      window=liste_window[iw]
      plot_obs_pred_window(X,pred.agreg,Y,date.valid,window,sta,ech,agreg)
    }
    plot_obs_pred(predictions.ech.agreg.sta,Y,date.valid,liste_window,sta,ech,agreg)
  }

  # quantile of the aggregation's absolute error for a station, a lead time and a window
  if(flag.quantiles){
    for (iw in 1:n.window){
      window=liste_window[iw]
      quantile_abs_error=eval_quantiles_abs_error(predictions.ech.agreg.sta[[iw]],Y)
      # cat to print in the file:
      cat("\n quantiles for aggregation ",agreg," station ",sta," lead time ",ech," window \n",file=file.out,append=T)
      for (iq in 1:length(quantile_abs_error)) cat(names(quantile_abs_error[iq])," : ",quantile_abs_error[iq],"\n",file=file.out,append=T)
    }
  }

  print(paste(agreg,"plot rmse",sep=" "))
  pdf(paste(dir_plot,dir_rmse,"rmse_",agreg,"_",ech,"_",sta,".pdf",sep=""))
  plot(x=liste_window,y=rmse.ech.agreg.sta)
  dev.off()

  print(paste("RMSE of the different sliding windows for : ",agreg,", station : ",sta,", ech : ",ech,sep=""))
  print(paste("windows : ",liste_window,sep=""))
  print(paste("RMSE : ",round(rmse.ech.agreg.sta,3),sep=""))
  # cat to print in the file:
  cat("\n RMSE of the different sliding windows for : ",agreg,", station : ",sta,", ech : ",ech,"\n",file=file.out,append=T)
  cat("windows : ",liste_window,"\n",file=file.out,append=T)
  cat("RMSE : ",round(rmse.ech.agreg.sta,3),"\n",file=file.out,append=T)

  if (flag.plot.ech.sta.window){
    for (iw in 1:n.window){ # for each window, plot of ...
      # weight plots
      plot_name=paste("weights_",agreg,"_window",liste_window[iw],"_ech",ech,"_",sta,"_",plot.periode[1],"_",plot.periode[2],".pdf",sep="")
      print(names.experts)
      lance_plot_weights(weights.ech.agreg.sta[[iw]],names.experts,plot_name,date.valid,plot.periode,T)

      # mixture object plot
      if (agreg!="UNIFORM") {
        mixture_graphs=c("weights","boxplots","squareloss","residuals","averageloss","contribution") # list of graphs which will be plotted
        print(paste("plot mixture object ",agreg,"_window",liste_window[iw],"_ech",ech,"_",sta,"_",start.date,"_",end.date,".pdf",sep=""))
        plot(mixture.ech.agreg.sta[[iw]], pause=FALSE)
        file.name.mixture=paste(dir_plot,dir_mixture,agreg,"_window",liste_window[iw],"_ech",ech,"_",sta,"_",start.date,"_",end.date,sep="")
        for (graph in mixture_graphs){
          file.rename(paste(dir_plot,dir_mixture,graph,".pdf",sep=""),paste(file.name.mixture,"_",graph,".pdf",sep="")) # rename the file
        }
      }

      # observation and prediction plots
      print(paste("plot ",dir_obs_pred,"pred_obs_",agreg,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".png",sep=""))
      png(paste(dir_plot,dir_obs_pred,"pred_obs_",agreg,"_ech",ech,"_",sta,"_",start.date,"_",end.date,".png",sep=""))
      plot(x=Y,y=predictions.ech.agreg.sta[[iw]])
      dev.off()

      # parameter plots
      if(flag.param){
        if (agreg == "FSBOACALIB" || agreg == "FSCALIB" || agreg == "FSCOCOBOACALIB"){
          plot_a_e(alpha.ech.agreg.sta[[iw]],"alpha",sta,liste_window[iw],start.date,end.date,names.experts)
          plot_a_e(eta.ech.agreg.sta[[iw]],"eta",sta,liste_window[iw],start.date,end.date,names.experts)
        } else if (agreg=="BOA" || agreg=="SOCOBOA") plot_a_e(eta.ech.agreg.sta[[iw]],"eta",sta,liste_window[iw],start.date,end.date,names.experts)
      }
    }
  }

  return(list(weights.ech.agreg.sta=weights.ech.agreg.sta,alpha.ech.agreg.sta=alpha.ech.agreg.sta,
              eta.ech.agreg.sta=eta.ech.agreg.sta,predictions.ech.agreg.sta=predictions.ech.agreg.sta,
              A1.ech.agreg.sta=A1.ech.agreg.sta,B1.ech.agreg.sta=B1.ech.agreg.sta,C1.ech.agreg.sta=C1.ech.agreg.sta,D1.ech.agreg.sta=D1.ech.agreg.sta,
              A2.ech.agreg.sta=A2.ech.agreg.sta,B2.ech.agreg.sta=B2.ech.agreg.sta,C2.ech.agreg.sta=C2.ech.agreg.sta,D2.ech.agreg.sta=D2.ech.agreg.sta,
              NERR.ech.agreg.sta=NERR.ech.agreg.sta))
}

# launch the aggregation and compute the mean rmse over the station,
# concatenate the alphas and etas of all the stations for the different windows
lance_agreg <- function(agreg){
  predictions.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for every station, contains for every window, the predictions of agreg
  alpha.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for every station, contains for every window, the alpha's of agreg
  eta.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for every station, contains for every window, the eta's of agreg
  obs.sta=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, will contain all the observations for all the stations
  A1.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  B1.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  C1.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  D1.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  A2.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  B2.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  C2.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS
  D2.ech.agreg=mylist <- sapply(paste("sta",as.character(liste_insee),sep=""),function(x) NULL) # list, for skill scores PSS

  ##########################################################################################################
  # start parallelisation over the stations
  registerDoMC(min(n.sta,30))
  # out_sta=foreach(ista=1:n.sta,.inorder = TRUE,.packages=c("foreach","doMC","parallel","randomForestSRC"))%dopar%{
  out_sta=foreach(ista=1:n.sta,.inorder = TRUE,.packages=c("foreach","doMC"))%do%{
  # for (ista in 1:n.sta){ print(" WARNING station parallelization disabled !!!!!!")
    file_station=liste_station[ista]
    sta=substr(strsplit(file_station, "_")[[1]][7],1,8) # insee of the station
    ech=strsplit(strsplit(liste_station[1], "/")[[1]][3],"_")[[1]][4] # lead time processed
    cat("\n")
    cat(paste("######################################## ech : ",ech,", agreg : ",agreg,", station : ",sta,"########################################\n"))
    cat("        temperature aggregation for all the sliding windows:\n")
    cat("\n######################################## ech : ",ech,", agreg : ",agreg,", station : ",sta,"########################################\n",file=file.out,append=T)
    cat("        temperature aggregation for all the sliding windows:\n",file=file.out,append=T)
    # retrieve data
    data=recup_data(file_station)
    date.valid=data$date
    X=data$X
    print(paste("we have ",nrow(X)," observations",sep=""))
    Y=data$Y
    start.date=data$start.date
    end.date=data$end.date

    ech=strsplit(strsplit(liste_station[1], "/")[[1]][3],"_")[[1]][4] # lead time processed

    nbr_exp=ncol(X)
    T=length(Y)
    liste_window=c(seq(60,420,60),520,620,720,seq(780,1120,60),T)

    n.window=length(liste_window)

    # launch the aggregation for one aggregation, one lead time and one station
    liste=lance_station(agreg,file_station, liste_window, X ,Y, T, date.valid,start.date, end.date) # list whose elements are the data associated with windows for a station

    return(list(predictions.ech.agreg.sta=liste$predictions.ech.agreg.sta,obs.sta=Y,eta.ech.agreg.sta=liste$eta.ech.agreg.sta,
                alpha.ech.agreg.sta=liste$alpha.ech.agreg.sta,liste_window=liste_window,start.date=start.date,
                end.date=end.date,date.valid=date.valid,
                A1.ech.agreg.sta=liste$A1.ech.agreg.sta,B1.ech.agreg.sta=liste$B1.ech.agreg.sta,C1.ech.agreg.sta=liste$C1.ech.agreg.sta,D1.ech.agreg.sta=liste$D1.ech.agreg.sta,
                A2.ech.agreg.sta=liste$A2.ech.agreg.sta,B2.ech.agreg.sta=liste$B2.ech.agreg.sta,C2.ech.agreg.sta=liste$C2.ech.agreg.sta,D2.ech.agreg.sta=liste$D2.ech.agreg.sta,
                NERR.ech.agreg.sta=liste$NERR.ech.agreg.sta))
  }
  # end parallelisation over the stations
  ##########################################################################################################

  for (ista in 1:n.sta){
    predictions.ech.agreg[[ista]]=out_sta[[ista]]$predictions.ech.agreg.sta
    obs.sta[[ista]]=out_sta[[ista]]$obs.sta
    eta.ech.agreg[[ista]]=out_sta[[ista]]$eta.ech.agreg.sta
    alpha.ech.agreg[[ista]]=out_sta[[ista]]$alpha.ech.agreg.sta
    A1.ech.agreg[[ista]]=out_sta[[ista]]$A1.ech.agreg.sta
    B1.ech.agreg[[ista]]=out_sta[[ista]]$B1.ech.agreg.sta
    C1.ech.agreg[[ista]]=out_sta[[ista]]$C1.ech.agreg.sta
    D1.ech.agreg[[ista]]=out_sta[[ista]]$D1.ech.agreg.sta
    A2.ech.agreg[[ista]]=out_sta[[ista]]$A2.ech.agreg.sta
    B2.ech.agreg[[ista]]=out_sta[[ista]]$B2.ech.agreg.sta
    C2.ech.agreg[[ista]]=out_sta[[ista]]$C2.ech.agreg.sta
    D2.ech.agreg[[ista]]=out_sta[[ista]]$D2.ech.agreg.sta
    NERR.ech.agreg[[ista]]=out_sta[[ista]]$NERR.ech.agreg.sta
  }
  liste_window=out_sta[[1]]$liste_window
  n.window=length(liste_window)
  start.date=out_sta[[1]]$start.date
  end.date=out_sta[[1]]$end.date
  date.valid=out_sta[[1]]$date.valid
  all_obs<-unlist(obs.sta) # for iech, concatenate all observations from all stations into one vector
  for (iw in 1:n.window) {
    window=liste_window[iw]
    all_iw_predictions=c()
    all_iw_A1=c()
    all_iw_B1=c()
    all_iw_C1=c()
    all_iw_D1=c()
    all_iw_A2=c()
    all_iw_B2=c()
    all_iw_C2=c()
    all_iw_D2=c()
    all_iw_NERR=c()
    for (ista in 1:n.sta) {
      pred_sta_iw=unlist(predictions.ech.agreg[[ista]][[iw]]) # we put the predictions in a vector
      A1_sta_iw=unlist(A1.ech.agreg[[ista]][[iw]])
      B1_sta_iw=unlist(B1.ech.agreg[[ista]][[iw]])
      C1_sta_iw=unlist(C1.ech.agreg[[ista]][[iw]])
      D1_sta_iw=unlist(D1.ech.agreg[[ista]][[iw]])
      A2_sta_iw=unlist(A2.ech.agreg[[ista]][[iw]])
      B2_sta_iw=unlist(B2.ech.agreg[[ista]][[iw]])
      C2_sta_iw=unlist(C2.ech.agreg[[ista]][[iw]])
      D2_sta_iw=unlist(D2.ech.agreg[[ista]][[iw]])
      NERR_sta_iw=unlist(NERR.ech.agreg[[ista]][[iw]])
      all_iw_predictions=c(all_iw_predictions,pred_sta_iw) # will contain all the predictions of iagreg with iw for all the stations
      all_iw_A1=c(all_iw_A1,A1_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_B1=c(all_iw_B1,B1_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_C1=c(all_iw_C1,C1_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_D1=c(all_iw_D1,D1_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_A2=c(all_iw_A2,A2_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_B2=c(all_iw_B2,B2_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_C2=c(all_iw_C2,C2_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_D2=c(all_iw_D2,D2_sta_iw) # for PSS of iagreg with iw for all the stations
      all_iw_NERR=c(all_iw_NERR,NERR_sta_iw) # for NERR of iagreg with iw for all the stations
    }
    A1=sum(all_iw_A1)
    B1=sum(all_iw_B1)
    C1=sum(all_iw_C1)
    D1=sum(all_iw_D1)
    A2=sum(all_iw_A2)
    B2=sum(all_iw_B2)
    C2=sum(all_iw_C2)
    D2=sum(all_iw_D2)
    NERR=sum(all_iw_NERR)
    hit_rate1=A1/(A1+C1)
    false_rate1=B1/(B1+D1)
    PSS1=hit_rate1-false_rate1
    false_ratio1=B1/(A1+B1)
    hit_rate2=A2/(A2+C2)
    false_rate2=B2/(B2+D2)
    PSS2=hit_rate2-false_rate2
    GSS=(PSS1+PSS2)/2
    false_ratio2=B2/(A2+B2)
    print(paste("skill scores for all stations ",agreg, ", lead time : ",ech,", window : ",window,
                ", A1 :",A1,", B1 :",B1,", C1 :",C1,", D1 :",D1,", hit_rate1 : ", hit_rate1,", false_rate1 : ", false_rate1,", PSS1 : ",PSS1 ,", false_ratio1 : ",false_ratio1,
                ", A2 :",A2,", B2 :",B2,", C2 :",C2,", D2 :",D2,", hit_rate2 : ", hit_rate2,", false_rate2 : ", false_rate2,", PSS2 : ",PSS2 ,", false_ratio2 : ",false_ratio2,
                ", NERR (nb wrong experts awake) : ",NERR," GSS : ",GSS,sep=""))
    rmse.agreg.ech.window=sqrt(mean(loss(all_iw_predictions, all_obs, loss.type = "square"),na.rm=T))
    cat("\n")
    print(paste("rmse for all stations, for aggregation ",agreg,", lead time ",ech,", window ",window," RMSE : ",round(rmse.agreg.ech.window,3),sep=""))
    cat("\n")
    cat("\n rmse for all stations, for aggregation ",agreg,", lead time ",ech,", window ",window,file=file.out,append=T)
    cat(round(rmse.agreg.ech.window,3),"\n",file=file.out,append=T)
    if (flag.quantiles){# quantiles over all the stations of the aggregation's absolute error for one window and one lead time
      quantile_abs_error=eval_quantiles_abs_error(all_iw_predictions,all_obs)
      cat("\n",file=file.out,append=T)
      cat("quantiles for all stations, for aggregation ",agreg,", lead time ",ech,", window \n",window,file=file.out,append=T)
      for (iq in length(quantile_abs_error)) cat(names(quantile_abs_error[iq])," : ",quantile_abs_error[iq],"\n",file=file.out,append=T)
    }
  }

  return(list(liste_window=liste_window,alpha.ech.agreg=alpha.ech.agreg,eta.ech.agreg=eta.ech.agreg,
          start.date=start.date,end.date=end.date,date.valid=date.valid, predictions.ech.agreg=predictions.ech.agreg,
          obs.sta=obs.sta,
          A1.ech.agreg=A1.ech.agreg,B1.ech.agreg=B1.ech.agreg,C1.ech.agreg=C1.ech.agreg,D1.ech.agreg=D1.ech.agreg,
          A2.ech.agreg=A2.ech.agreg,B2.ech.agreg=B2.ech.agreg,C2.ech.agreg=C2.ech.agreg,D2.ech.agreg=D2.ech.agreg,
          NERR.ech.agreg=NERR.ech.agreg))
}

# recover the experts predictions and the observations of the data and the start and end date of the data  
recup_data <-function(file_station){
    # retrieve data
    donnee = read.table(file_station,head=T,sep=";",na.strings=c("       NA","      NA","     NA","    NA","   NA","  NA"," NA","NA")) # there are one or more spaces before the NAs..., hence "  NA"
    names(donnee)[which(names(donnee)=="Q20")]="Q30" # TODO hack !!!!!! it's the Q30 and not the Q20, it would be better to change this in the code on sotrtm35-sidev which creates the .txt
    
    colnames(donnee)[which(colnames(donnee)=="sd.aro")]="raw.aro"
    colnames(donnee)[which(colnames(donnee)=="as.aro")]="mos.aro"
    colnames(donnee)[which(colnames(donnee)=="sd.arp")]="raw.arp"
    colnames(donnee)[which(colnames(donnee)=="as.arp")]="mos.arp"
    colnames(donnee)[which(colnames(donnee)=="sd.cep")]="raw.cep"
    colnames(donnee)[which(colnames(donnee)=="as.cep")]="mos.cep"
    
    donnee=donnee[1:1253,]

    donnee[which(is.na(donnee$obs.t)),5:17]<-NA # if observation missing, observation and all the experts put to NA
    start.date=substr(donnee[,4],1,10)[1] # first day of the data
    end.date=substr(donnee[,4],1,10)[nrow(donnee)] # last day of the data

    X=donnee[,liste_exp]

    colnames(X)= colnames(donnee[,liste_exp])


    X=rm_missing_col(X,0.05)
    Y=donnee$obs.t
    Y=as.numeric(Y)
    Y=c(Y)
    expert.names = list(row=rep(NULL,ncol(donnee)),col=names(X))
    X=as.matrix(X, dimnames=expert.names)
    date.valid=donnee$valid

    .GlobalEnv$l_max <- max(2*(X-Y)*X,na.rm=T)

    return(list(X=X,Y=Y,start.date=start.date,end.date=end.date,date.valid=date.valid))
}

#function to remove columns of a data frame X with less than percentage% missing values
rm_missing_col <- function(X,percentage) {
    rm.col=c()
    for (ncol in 1:ncol(X)){
      nna=length(which(X[,ncol]==-100))+length(which(is.na(X[,ncol])))+length(which(X[,ncol]==-999)) #sum of the NA's and -100
      if (nna>nrow(X)*percentage){rm.col <- c(rm.col,ncol)}
    }
    if (!is.null(rm.col)){X<-X[,-rm.col]} #keep only the experts with less than percentage % of missing values
    return(X)
}

#remove the 60 first iterations of the parameters because the aggregation isn't stabilized at the begining
rm.param.window <-function(data.param,window){
  if (!is.null(ncol(data.param))) {
    data.param=data.param[-c(1:60),] #if there are multiple parameters (one for each expert for example)
  } else data.param=data.param[-c(1:60)] #if the data isn't empty and has just one parameter
  return(data.param)
}


########################################################
########################################################





########################################################
########################################################
########################################################
######################## MAIN ##########################

########################################################################
############################ debut variables ###########################
set.seed(777)# to fix randomness
file.out="results.txt" #file in which will be written the prints and scores

flg.map=TRUE #plot a map of france with the stations on it
flag.rm.param.window=TRUE #remove the first 60 alphas and etas, (aggregation isn't stabilized at the begining, and huge eta0)
flag.param=FALSE #plot or not of the parameters
flag.quantiles=TRUE #plot and compute the quantiles of absolute error
flg.with=TRUE #if data with pearp, TRUE->with pearp, FALSE-> without pearp
flag.obs.pred=FALSE #plot or not of predictions and observations
flag.plot.ech.sta.window=FALSE # plot mixture objects,weights,... for one lead time one station one agregation and one window
flag.write.txt=FALSE #write the observations, the predictions and the weights in a txt file
flag.save=FALSE #save some .RData

#start and end date of the data ploted on the obs_pred... plots
plot.periode=c("2021-12-01","2022-02-01") #pour evenement de Chamonix juste avant Noël 2021

#~best alpha and eta calibrated by FSCALIB for the data without pearp
alpha=0.01#share coefficient for FS,FSBOA and maybe SMH
eta=3.125000e-02

liste_exp_without=c("sd.aro","as.aro","sd.arp","as.arp","sd.cep","as.cep") # liste des experts non biaises
liste_exp_with=c("raw.aro","mos.aro","raw.arp","mos.arp","raw.cep","mos.cep","Q10","Q30","Q50","Q70","Q90") #experts to use if available in the used data

parametre="t"
res="00"

#liste_ech=c("15","30","42","45","48","57","72","84")
lead_times=c("3","6","9","12","15","18","21","24","27","30","33","36","39","42","45","48","57","72","84") #needed for the function recupLeadTimes (especially when liste_ech!=lead_times)
liste_ech=c("3","6","9","12","15","18","21","24","27","30","33","36","39","42","45","48","57","72","84")
liste_ech=c("6","9","12","15","18","21","24","27","30","33","36","39","42","45","48","57","72","84") #without lead time 3 because we use the observation of lead time 3 in the random forests

n.ech=length(liste_ech)

step_window=20 #step between every tested window

#list of the agregations to run
liste_agreg <- c("BOA","EWA","MLpol","MLprod","NBM","bestConvex","bestExpert","UNIFORM")
liste_agreg <- c("BOA","MLpol","MLprod","EWA","NBM")

file.out=paste("results_",liste_agreg[1],".txt",sep="") #file in which will be written the prints and scores

n.agreg=length(liste_agreg)

liste_couleur <- c("red","blue","purple","green","yellow","black","pink","orange","azure4","deeppink")
obs.all=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)#list, will contain all the observations for all the lead times and stations
predictions.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, will contain all the predictions of all the agregations for each lead time, station and window
A1.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
B1.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
C1.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
D1.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
A2.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
B2.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
C2.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
D2.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
NERR.all=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, for skill scores of all the agregations for each lead time, station and window
rmse.boxplot=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, will contain the rmse of all the agregations for each lead time, station
rmse.TOT=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, will contain the rmse over all the stations and lead times for all the agregations and windows
rmse.TOT.lim=c(1)

l.agreg.param=intersect(liste_agreg,c("BOA","SOCOBOA","FSSOCOBOA","FSSOCOBOACALIB","FSBOA","FSBOACALIB","FSCALIB","SMHCALIB","SMBOA","SMBOACALIB"))#aggregation with parameters to plot
param.TOT=sapply(l.agreg.param, function(x) NULL)#list with empty elements the names of the agregations with parameters, will contain for all the aggregations all the windows, all the param of all the lead times
flag.param.TOT.init=FALSE #FALSE until param.TOT is initialized
########################### fin variables ##############################
########################################################################


#########################################################################
############################## start ech ################################
for (iech in 1:n.ech){
  ech=liste_ech[iech]
  cat("\n")
  cat("################################################################################\n")
  cat(paste("################################### ech ",ech,"####################################\n"),sep="")
  cat("################################################################################\n")
  cat("\n################################################################################\n",file=file.out,append=T)
  cat("################################### ech ",ech,"####################################\n",file=file.out,append=T)
  cat("################################################################################\n",file=file.out,append=T)
  
  liste_recap=list()
  rmse.ech=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL) #list, will contain for all the agregations and windows the rmse over all the stations
  rmse.lim=c(1) #non empty, to deal with missing observations, NA's
  liste_station=dir(path = dir_data_with, pattern = paste("donnee_",parametre,"_",res,"_",ech,"_",sep=""),full.names=TRUE)



  #plot a map of France with the stations on it
  if (iech==1 & flg.map) plot_map(liste_station)

  if(substr(liste_station[1],1,22)==dir_data_without) {
    liste_exp=liste_exp_without #if data without pearp
  } else if(substr(liste_station[1],1,22)==dir_data_with & !flg.with) {
    liste_exp=liste_exp_with[1:6] #if data with pearp but flag for pearp is FALSE
  } else if (substr(liste_station[1],1,22)==dir_data_with & flg.with) {
    liste_exp=liste_exp_with #if data with pearp and flag for pearp is TRUE
  } else (stop("PROBLEM (with the path ? or the lead times ?)"))


  n.sta=length(liste_station)
  cat("                Liste des stations :")
  cat("\n")
  print(liste_station)
  cat("\n")
  #cat to print in file
  cat("\n LISTE DES STATIONS :\n",file=file.out,append=T)
  for (ista in 1:n.sta) cat(liste_station[ista],"\n",file=file.out,append=T)

  #insee of the stations, following the number of numbers of ech
  if(as.integer(ech)<10) {
    liste_insee=substr(liste_station,38,45)
  } else if( 10 <= as.integer(ech) & as.integer(ech)<100) {
    liste_insee=substr(liste_station,39,46)
  } else if( 100 <=as.integer(ech) ) liste_insee=substr(liste_station,40,47)

  #list, for iech and iagreg, will contain all the predictions for all the stations and windows
  pred.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  A1.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  B1.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  C1.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  D1.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  A2.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  B2.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  C2.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  D2.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  NERR.ech.agreg=mylist <- sapply(paste(liste_agreg,sep=""),function(x) NULL)
  for (iagreg in 1:n.agreg){
    #for every iagreg, one list for every station
    pred.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    A1.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    B1.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    C1.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    D1.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    A2.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    B2.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    C2.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    D2.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
    NERR.ech.agreg[[iagreg]]=mylist <- sapply(paste("sta_",liste_insee,sep=""),function(x) NULL)
  }

  #############################################################
  ##################### start iagreg ##########################
  for (iagreg in 1:n.agreg){
    cat("####################################################################################\n")
    cat(paste("########################## ech : ",ech,", agreg : ",liste_agreg[iagreg],"################################\n"),sep="")
    agreg=liste_agreg[iagreg]
    liste=lance_agreg(agreg) #lunch the agregations for one agregation and one lead time
    start.date=liste$start.date
    end.date=liste$end.date
    date.valid=liste$date.valid
    obs.sta=liste$obs.sta #obs for iech, of all the stations, ordered in the stations order
    n.iteration=length(obs.sta[[1]]) #number of iterations, should be the same for all the stations...
    n.window=length(liste$liste_window) #number of sliding windows
    liste_window=liste$liste_window
    liste_recap=c(liste_recap,list(liste))
    #list_recap, list of list
    #for aggregation strategy iagreg, one list with:
    # - liste_window, vector of the sliding windows,
    # - alpha.ech.agreg, alpha of all stations for all windows (length=number of station*iterations),
    # - eta.ech.agreg, eta of all stations for all windows (length=number of station*iterations)
    # - liste$predictions.ech.agreg, list, predictions of iagreg at iech for all the stations and all the windows
    # - liste$A B C D , list, to compute skill scores PSS of iagreg at iech for all the stations and all the windows
    # - obs.sta, list, observations of all the stations at lead time ech
    # - date.valid, date de validite (lead dates)
    # - start.date, start date of the data
    # - end.date, end date of the data
    for (ista in 1:n.sta){
      pred.ech.agreg[[iagreg]][[ista]]=liste$predictions.ech.agreg[[ista]] #add of the predictions of agreg for all the windows
      A1.ech.agreg[[iagreg]][[ista]]=liste$A1.ech.agreg[[ista]] #add of A of agreg for all the windows
      B1.ech.agreg[[iagreg]][[ista]]=liste$B1.ech.agreg[[ista]] #add of B of agreg for all the windows
      C1.ech.agreg[[iagreg]][[ista]]=liste$C1.ech.agreg[[ista]] #add of C of agreg for all the windows
      D1.ech.agreg[[iagreg]][[ista]]=liste$D1.ech.agreg[[ista]] #add of D of agreg for all the windows
      A2.ech.agreg[[iagreg]][[ista]]=liste$A2.ech.agreg[[ista]] #add of A of agreg for all the windows
      B2.ech.agreg[[iagreg]][[ista]]=liste$B2.ech.agreg[[ista]] #add of B of agreg for all the windows
      C2.ech.agreg[[iagreg]][[ista]]=liste$C2.ech.agreg[[ista]] #add of C of agreg for all the windows
      D2.ech.agreg[[iagreg]][[ista]]=liste$D2.ech.agreg[[ista]] #add of D of agreg for all the windows
      NERR.ech.agreg[[iagreg]][[ista]]=liste$NERR.ech.agreg[[ista]] #add of NERR of agreg for all the windows
    }
  }
  ###################### end iagreg ###########################
  #############################################################

  names(liste_recap)<-liste_agreg

  ################################
  #for one station one lead time and one window, plot of the predictions of all the agregations
  if (flag.obs.pred){
    for (ista in 1:n.sta){
      sta=substr(liste_station[ista],39,46)
      for (iagreg in 1:n.agreg){
        for (iw in 1:length(pred.ech.agreg[[iagreg]][[ista]])){#for every window
          window=names(pred.ech.agreg[[iagreg]][[ista]])[iw]
          data.plot=NULL
          for (iagregbis in 1:n.agreg){
            data.plot=cbind(data.plot,pred.ech.agreg[[iagregbis]][[ista]][[iw]])#table with the predictions of all the agregations for ista and iw
          }
          colnames(data.plot)<-liste_agreg
          plot_obs_agreg(data.plot,obs.sta[[ista]],date.valid,window,sta,ech)
        }
      }
    }
  }
  ################################
  #create a list, with all the observations for all the lead times and stations
  #and a list with all the corresponding predictions for all the agregations and windows and lead times
  obs.all[[iech]]=obs.sta#we add the observations at ech of the lead time ech of all the stations
  for (iagreg in 1:n.agreg){
    if(is.null(predictions.all[[iagreg]])) predictions.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(A1.all[[iagreg]])) A1.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(B1.all[[iagreg]])) B1.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(C1.all[[iagreg]])) C1.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(D1.all[[iagreg]])) D1.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(A2.all[[iagreg]])) A2.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(B2.all[[iagreg]])) B2.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(C2.all[[iagreg]])) C2.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(D2.all[[iagreg]])) D2.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    if(is.null(NERR.all[[iagreg]])) NERR.all[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    for (iw in 1:n.window){
      if(is.null(predictions.all[[iagreg]][[iw]])) predictions.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(A1.all[[iagreg]][[iw]])) A1.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(B1.all[[iagreg]][[iw]])) B1.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(C1.all[[iagreg]][[iw]])) C1.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(D1.all[[iagreg]][[iw]])) D1.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(A2.all[[iagreg]][[iw]])) A2.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(B2.all[[iagreg]][[iw]])) B2.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(C2.all[[iagreg]][[iw]])) C2.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(D2.all[[iagreg]][[iw]])) D2.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      if(is.null(NERR.all[[iagreg]][[iw]])) NERR.all[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
      for (ista in 1:n.sta){
        predictions.all[[iagreg]][[iw]][[iech]][[ista]]=pred.ech.agreg[[iagreg]][[ista]][[iw]]#we add the predictions of iagreg for ista, ech with iw
        A1.all[[iagreg]][[iw]][[iech]][[ista]]=A1.ech.agreg[[iagreg]][[ista]][[iw]]#we add the A of iagreg for ista, ech with iw
        B1.all[[iagreg]][[iw]][[iech]][[ista]]=B1.ech.agreg[[iagreg]][[ista]][[iw]]#we add the B of iagreg for ista, ech with iw
        C1.all[[iagreg]][[iw]][[iech]][[ista]]=C1.ech.agreg[[iagreg]][[ista]][[iw]]#we add the C of iagreg for ista, ech with iw
        D1.all[[iagreg]][[iw]][[iech]][[ista]]=D1.ech.agreg[[iagreg]][[ista]][[iw]]#we add the D of iagreg for ista, ech with iw
        A2.all[[iagreg]][[iw]][[iech]][[ista]]=A2.ech.agreg[[iagreg]][[ista]][[iw]]#we add the A of iagreg for ista, ech with iw
        B2.all[[iagreg]][[iw]][[iech]][[ista]]=B2.ech.agreg[[iagreg]][[ista]][[iw]]#we add the B of iagreg for ista, ech with iw
        C2.all[[iagreg]][[iw]][[iech]][[ista]]=C2.ech.agreg[[iagreg]][[ista]][[iw]]#we add the C of iagreg for ista, ech with iw
        D2.all[[iagreg]][[iw]][[iech]][[ista]]=D2.ech.agreg[[iagreg]][[ista]][[iw]]#we add the D of iagreg for ista, ech with iw
        NERR.all[[iagreg]][[iw]][[iech]][[ista]]=NERR.ech.agreg[[iagreg]][[ista]][[iw]]#we add the NERR of iagreg for ista, ech with iw
      }
      names(predictions.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(A1.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(B1.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(C1.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(D1.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(A2.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(B2.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(C2.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(D2.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
      names(NERR.all[[iagreg]][[iw]][[iech]])=paste("sta",liste_insee,sep="")
    }
  }
  ################################

  ################################
  #skill scores of the random forest over all the stations of the agregation for window iw and lead time iech
  for (iw in 1:n.window) {
    for (iagreg in 1:n.agreg) {
      A1=sum(unlist(A1.all[[iagreg]][[iw]][[iech]]))
      B1=sum(unlist(B1.all[[iagreg]][[iw]][[iech]]))
      C1=sum(unlist(C1.all[[iagreg]][[iw]][[iech]]))
      D1=sum(unlist(D1.all[[iagreg]][[iw]][[iech]]))
      A2=sum(unlist(A2.all[[iagreg]][[iw]][[iech]]))
      B2=sum(unlist(B2.all[[iagreg]][[iw]][[iech]]))
      C2=sum(unlist(C2.all[[iagreg]][[iw]][[iech]]))
      D2=sum(unlist(D2.all[[iagreg]][[iw]][[iech]]))
      NERR=sum(unlist(NERR.all[[iagreg]][[iw]][[iech]]))
      hit_rate1=A1/(A1+C1)
      false_rate1=B1/(B1+D1)
      PSS1=hit_rate1-false_rate1
      false_ratio1=B1/(A1+B1)
      hit_rate2=A2/(A2+C2)
      false_rate2=B2/(B2+D2)
      PSS2=hit_rate2-false_rate2
      false_ratio2=B2/(A2+B2)
      GSS=(PSS1+PSS2)/2
      print(paste("skill scores de RF : ",liste_agreg[iagreg],", ech : ",ech,", window : ",liste_window[iw],", pour toutes les stations,
            A1 :",A1,", B1 :", B1,", C1 :",C1,", D1 :",D1,", hit_rate1 : ", hit_rate1,", false_rate1 : ", false_rate1,", PSS1 : ",PSS1 ,", false_ratio1 : ",false_ratio1,
            "A2 :",A2,", B2 :",B2,", C2 :",C2,", D2 :",D2,", hit_rate2 : ", hit_rate2,", false_rate2 : ", false_rate2,", PSS2 : ",PSS2 ,", false_ratio2 : ",false_ratio2,
            ", NERR (nb wrong experts awake) : ",NERR," GSS : ",GSS,sep=""))
    }
  }
  ################################

  ################################
  #quantiles and RMSE over all the stations of the agregation's absolute error for window iw and lead time iech
  cat("\n")
  print("######################################")
  all.obs.ech<-unlist(obs.all[[iech]]) #on concatene toutes les obs de toutes les station dans un vecteur
  for (iw in 1:n.window) {
    for (iagreg in 1:n.agreg) {
      all.pred.ech=unlist(predictions.all[[iagreg]][[iw]][[iech]])
      rmse.agreg.ech.window=sqrt(mean(loss(all.pred.ech, all.obs.ech, loss.type = "square"),na.rm=T))
      cat("\n")
      print(paste("rmse de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                  ", ech : ",ech,", pour toutes les stations",sep=""))
      print(round(rmse.agreg.ech.window,3))
      cat("\n")
      #cat to print in file :
      cat("\n rmse de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                  ", ech : ",ech,", pour toutes les stations \n",file=file.out,append=T)
      cat(round(rmse.agreg.ech.window,3),file=file.out,append=T)
      if(is.null(rmse.ech[[iagreg]])) rmse.ech[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
      rmse.ech[[iagreg]][[iw]]=rmse.agreg.ech.window
      rmse.lim=range(c(rmse.lim,unlist(rmse.ech)),na.rm=T)# rmse.lim : min and max of the rmse for all the models
      if (flag.quantiles){
        quantile_abs_error=eval_quantiles_abs_error(all.pred.ech,all.obs.ech)
        cat("\n quantiles de l'erreur absolue des previsions de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                    ", ech : ",ech,", pour toutes les stations \n",file=file.out,append=T)
        for (iq in length(quantile_abs_error)) cat(names(quantile_abs_error[iq])," : ",quantile_abs_error[iq],"\n",file=file.out,append=T)
      }
    }
  }
  print("######################################")
  cat("\n")
  ################################

  ################################
  #plot of the rmse depending on the window for each agregation
  plot_rmse_ech(liste_window,rmse.ech,rmse.lim)
  ################################

  ################################
  #if not already done, initialization of param.TOT which will contain for all the aggregations
  # all the windows and all the parameters of all the stations and lead times
  lwindow=liste_recap[[1]]$liste_window
  lpar=c("eta","alpha")
  if (!flag.param.TOT.init & length(l.agreg.param)>0){
    for (iagreg in 1:length(l.agreg.param)){
      param.TOT[[iagreg]]=sapply(1:n.window, function(x) NULL) #add all the windows elements to each "agreg" element
      names(param.TOT[[iagreg]])=paste("window",as.character(lwindow),sep="")
      for (iwindow in 1:n.window) {
        param.TOT[[iagreg]][[iwindow]]=sapply(lpar, function(x) NULL)
      }
    }
    flag.param.TOT.init=TRUE
  }
  if (flag.param){
    #plot of the parameters for the lead time ech and concatenation of the parameters of all the lead times
    for (iagreg in 1:length(l.agreg.param)){#for the aggregations which are available and have parameters
      for (iwindow in 1:n.window ){ #for each window (normaly, same windows for every aggregation strategy)
        for (ipar in 1:length(lpar)){
          data.param.ech=eval(parse(text=paste("liste_recap$",l.agreg.param[iagreg],"$",lpar[ipar],".ech.agreg$window",lwindow[iwindow],sep="")))
          if (!is.null(data.param.ech)){ #if the param is available for the aggregation strategy
            window=eval(parse(text=paste("liste_recap$",agreg,"$liste_window[",iwindow,"]",sep="")))
            plot_parameters(data.param.ech,lpar[ipar],window,ech,n.sta)
            boxplot_parameters(data.param.ech,lpar[ipar],window,ech)
            #we add the parameters of the lead time ech to param.TOT$agreg$window$param
            if (is.null(ncol(data.param.ech))) {
              param.TOT[[iagreg]][[iwindow]][[ipar]]=c(param.TOT[[iagreg]][[iwindow]][[ipar]],data.param.ech)
            } else param.TOT[[iagreg]][[iwindow]][[ipar]]=rbind.fill(as.data.frame(param.TOT[[iagreg]][[iwindow]][[ipar]]),as.data.frame(data.param.ech))
          }
        }
      }
    }
  }
  ################################
}
############################## end ech ##################################
#########################################################################

nameRdata=paste("results_",liste_agreg[1],".Rdata",sep="") 
save.image(nameRdata)
#system(paste("mv ",nameRdata," ",dir_utemp,nameRdata,sep=""))



print("###############################################################")
print("##############  END loop over the lead times ##################")
print("###############################################################")

cat("\n###############################################################\n",file=file.out,append=T)
cat("##############  END loop over the lead times ##################\n",file=file.out,append=T)
cat("###############################################################\n",file=file.out,append=T)

#RMSE of sta over all the lead times for one window
for (iagreg in 1:n.agreg){
  for (iw in 1:n.window){
    for (ista in 1:n.sta){
      sta_agreg_window=c()
      sta_obs=c()
      for (iech in 1:n.ech){
        sta_agreg_window=c(sta_agreg_window,predictions.all[[iagreg]][[iw]][[iech]][[ista]])#we add the predictions of iagreg iw iech ista
        sta_obs=c(sta_obs,obs.all[[iech]][[ista]])#we add the observations of iech ista
      }
      rmse_sta= round(sqrt(mean((sta_agreg_window - sta_obs)^2,na.rm=TRUE)),5)
      print(paste("RMSE : ",rmse_sta," of ",liste_agreg[iagreg]," with the sliding window ",liste_window[iw]," at station ",liste_insee[ista]," over all the lead times",sep=""))
    }
  }
}

#RMSE of ech over all the stations for one window
for (iagreg in 1:n.agreg){
  for (iw in 1:n.window){
    for (iech in 1:n.ech){
      ech_agreg_window=c()
      ech_obs=c()
      for (ista in 1:n.sta){
        ech_agreg_window=c(ech_agreg_window,predictions.all[[iagreg]][[iw]][[iech]][[ista]])#we add the predictions of iagreg iw iech ista
        ech_obs=c(ech_obs,obs.all[[iech]][[ista]])#we add the observations of iech ista
      }
      rmse_ech= round(sqrt(mean((ech_agreg_window - ech_obs)^2,na.rm=TRUE)),5)
      print(paste("RMSE : ",rmse_ech," of ",liste_agreg[iagreg]," with the sliding window ",liste_window[iw]," at ech ",liste_ech[iech]," over all the stations",sep=""))
    }
  }
}

if(flag.save){
  fname=paste("predictions.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=predictions.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("obs.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=obs.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("A1.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=A1.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("B1.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=B1.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("C1.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=C1.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("D1.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=D1.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
    fname=fname=paste("A2.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=A2.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("B2.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=B2.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("C2.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=C2.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("D2.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=D2.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
  fname=fname=paste("NERR.all_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=NERR.all #save() can only save "entire" objects and not one column of an object like object$toto
  save(dataToSave,file=fname)
}


#rmse and quantiles over all the stations and lead times of the agregation's absolute error for a window
print("######################################")
all.obs<-unlist(obs.all) #on concatene toutes les obs de toutes les stations dans un vecteur
for (iagreg in 1:n.agreg) {
  for (iw in 1:n.window) {
    all.pred=unlist(predictions.all[[iagreg]][[iw]]) #we put all the predictions for all the lead times and stations of (iagreg,iw) in a vector
    rmse.agreg.window=sqrt(mean(loss(all.pred, all.obs, loss.type = "square"),na.rm=T)) #rmse over all the stations and lead times for iagreg and iw
    bias.agreg.window=mean(all.pred - all.obs,na.rm=T) #bias over all the stations and lead times for iagreg and iw
    print(paste("rmse de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                ", pour toutes les stations et echeances",sep=""))
    print(rmse.agreg.window)
    print(paste("bias de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                ", pour toutes les stations et echeances",sep=""))
    print(bias.agreg.window)
    #cat to print in file :
    cat("\n rmse de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                ", pour toutes les stations et echeances \n",file=file.out,append=T)
    cat(rmse.agreg.window,file=file.out,append=T)
    cat("\n bias de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                ", pour toutes les stations et echeances \n",file=file.out,append=T)
    cat(bias.agreg.window,file=file.out,append=T)
    if(is.null(rmse.TOT[[iagreg]])) rmse.TOT[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
    rmse.TOT[[iagreg]][[iw]]=rmse.agreg.window
    rmse.TOT.lim=range(c(rmse.TOT.lim,unlist(rmse.TOT)),na.rm=T)# rmse.lim : min and max of the rmse for all the models
    if (flag.quantiles){
      quantile_abs_error=eval_quantiles_abs_error(all.pred,all.obs)
      cat("\n quantiles de l'erreur absolue des previsions de l'agregation : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                  ", pour toutes les stations et echeances \n",file=file.out,append=T)
      for (iq in 1:length(quantile_abs_error)) cat("\n ",names(quantile_abs_error[iq])," : ",quantile_abs_error[iq],file=file.out,append=T)
      print(quantile_abs_error[length(quantile_abs_error)-1])
    }
  }
}
print("######################################")

#skill scores over all the stations and lead times of the RF for a window
print("######################################")
for (iagreg in 1:n.agreg) {
  for (iw in 1:n.window) {
    all.A1=unlist(A1.all[[iagreg]][[iw]])
    all.B1=unlist(B1.all[[iagreg]][[iw]])
    all.C1=unlist(C1.all[[iagreg]][[iw]])
    all.D1=unlist(D1.all[[iagreg]][[iw]])
    all.A2=unlist(A2.all[[iagreg]][[iw]])
    all.B2=unlist(B2.all[[iagreg]][[iw]])
    all.C2=unlist(C2.all[[iagreg]][[iw]])
    all.D2=unlist(D2.all[[iagreg]][[iw]])
    all.NERR=unlist(NERR.all[[iagreg]][[iw]])
    A1=sum(all.A1)
    B1=sum(all.B1)
    C1=sum(all.C1)
    D1=sum(all.D1)
    A2=sum(all.A2)
    B2=sum(all.B2)
    C2=sum(all.C2)
    D2=sum(all.D2)
    NERR=sum(all.NERR)
    hit_rate1=A1/(A1+C1)
    false_rate1=B1/(B1+D1)
    PSS1=hit_rate1-false_rate1
    false_ratio1=B1/(A1+B1)
    hit_rate2=A2/(A2+C2)
    false_rate2=B2/(B2+D2)
    PSS2=hit_rate2-false_rate2
    false_ratio2=B2/(A2+B2)
    GSS=(PSS1+PSS2)/2
    print(paste("skill scores of the RF : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],
                ", pour toutes les stations et echeances",
                ", A1 :",A1,", B1 :",B1,", C1 :",C1,", D1 :",D1,", hit_rate1 : ", hit_rate1,", false_rate1 : ", false_rate1,", PSS1 : ",PSS1 ,", false_ratio1 : ",false_ratio1,
                ", A2 :",A2,", B2 :",B2,", C2 :",C2,", D2 :",D2,", hit_rate2 : ", hit_rate2,", false_rate2 : ", false_rate2,", PSS2 : ",PSS2 ,", false_ratio2 : ",false_ratio2,
                ", NERR (nb wrong experts awake) : ",NERR," GSS : ",GSS,sep=""))
    cat("\n skill scores of the RF : ",liste_agreg[iagreg],", fenêtre : ",liste_window[iw],", pour toutes les stations et echeances",
                ", A1 :",A1,", B1 :",B1,", C1 :",C1,", D1 :",D1,", hit_rate1 : ", hit_rate1,", false_rate1 : ", false_rate1,", PSS1 : ",PSS1 ,", false_ratio1 : ",false_ratio1,
                ", A2 :",A2,", B2 :",B2,", C2 :",C2,", D2 :",D2,", hit_rate2 : ", hit_rate2,", false_rate2 : ", false_rate2,", PSS2 : ",PSS2 ,", false_ratio2 : ",false_ratio2,
                ", NERR (nb wrong experts awake) : ",NERR," GSS : ",GSS,
                "\n",file=file.out,append=T)
  }
}
print("######################################")



#plot of the parameters off all the lead times (for each aggregation, and each window)
if (flag.param){
  for (iagreg in 1:length(l.agreg.param)){#for the aggregations which are available and have parameters
    for (iwindow in 1:n.window ){ #for each window (normaly, same windows for every aggregation strategy)
      for (ipar in 1:length(lpar)){
        data.param=eval(parse(text=paste("param.TOT$",l.agreg.param[iagreg],"[[",iwindow,"]][[",ipar,"]]",sep="")))
        if (!is.null(param.TOT[[iagreg]][[iwindow]][[ipar]])){ #if the param is available for the aggregation strategy
          window=names(param.TOT[[iagreg]])[iwindow]
          plot_parameters(param.TOT[[iagreg]][[iwindow]][[ipar]],lpar[ipar],window,"AllEch",n.sta)
          boxplot_parameters(param.TOT[[iagreg]][[iwindow]][[ipar]],lpar[ipar],window,"AllEch")
        }
      }
    }
  }
}

boxplot_colors=c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999", "#0072B2", "#D55E00", "#CC79A7")[1:length(liste_agreg)]


################################
#for the different aggregations, plot of the rmse over the lead times and stations
df<-cbind(liste_window,round(unlist(rmse.TOT[[1]],2),5))
df=as.data.frame(df)
if(length(liste_agreg)>1){for (icol in 2:length(liste_agreg)) df<-cbind(df,round(unlist(rmse.TOT[[icol]]),5))}
colnames(df)[2:ncol(df)]=liste_agreg
df<-df %>%
    dplyr::select(liste_window, c(all_of(liste_agreg))) %>%
    gather(key="modele", value="RMSE",-liste_window)

df <- transform(df, modele=factor(modele, levels=unique(modele))) #without this, ggplot don't plot the lines in the expected order (the order of the dataframe's columns)
width.values=rep(2,length(liste_agreg))
names(width.values)=c(liste_agreg)
line.type.values=rep("solid",length(liste_agreg))
names(line.type.values) = liste_agreg

plot.df<-ggplot(df, aes(x=liste_window, y=RMSE))+
        #theme_gray() +
        theme_bw() + #theme, arriere plan... !!!!! le mettre au début, sinon ça peut ecraser le reste !!!!!!!!!!!
        geom_line(aes(color = modele, size = modele ,linetype=modele)) + #to have the lines
        geom_point(aes(color = modele, size = modele),show.legend = FALSE,size=4) + #to have the points/shapes
        scale_x_continuous(breaks=liste_window,labels = liste_window)+
        scale_color_manual(values=boxplot_colors)+
        scale_linetype_manual(values=line.type.values) + #solid or dashed lines, and no legend
        #scale_shape_manual(values = shape.values) + #type of shape for the points
        #scale_colour_manual(values=color.values) + #colours
        scale_size_manual(values = width.values) + #sizes
        labs(x="size of the window (days)", y= "RMSE (°C)") + #title, x and y axes description
        theme(plot.title = element_text(colour="black", size=20), axis.title.x = element_text(colour="black",size=20),
          axis.title.y = element_text(colour="black",size=20),legend.text=element_text(size=20),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15)) +
        theme(legend.title=element_blank())

ggsave(file=paste(dir_plot,dir_rmse,"rmse_resume_tot_",start.date,"_",end.date,".pdf",sep=""), width=18, height=7, dpi=600)
print(paste(dir_plot,dir_rmse,"rmse_resume_tot_",start.date,"_",end.date,".pdf has been created",sep=""))
################################
################################################################

################################################################
################################
#create a list, with the rmse of all the stations for one window one lead time and one aggregation
for (iagreg in 1:n.agreg) {
  if(is.null(rmse.boxplot[[iagreg]])) rmse.boxplot[[iagreg]]=mylist <- sapply(paste("w",as.character(liste_window),sep=""),function(x) NULL)
  for (iw in 1:n.window) {
    if(is.null(rmse.boxplot[[iagreg]][[iw]])) rmse.boxplot[[iagreg]][[iw]]=mylist <- sapply(paste("ech",liste_ech,sep=""),function(x) NULL)
    for (iech in 1:n.ech) {
      for (ista in 1:n.sta) {
        rmse.boxplot[[iagreg]][[iw]][[iech]]=c(rmse.boxplot[[iagreg]][[iw]][[iech]],
                                              sqrt(mean((predictions.all[[iagreg]][[iw]][[iech]][[ista]]-obs.all[[iech]][[ista]])^2,na.rm=T)))#we add the predictions of iagreg for ista, ech with iw
      }
    }
  }
}
################################
#rmse boxplots, one for each iteration window and aggregation, the boxplots represent the variability between the stations
boxplot_ech=rep(rep(liste_ech, each=n.sta),n.agreg)
boxplot_ech=as.numeric(boxplot_ech)
aggregation=rep(rep(liste_agreg,each=n.sta),each=n.ech)
aggregation[which(aggregation=="FS")]<-"MH"
aggregation[which(aggregation=="FSBOA")]<-"MBOA"
for (iw in 1:n.window) {
  #boxplots for each window
  boxplot_rmse=c()
  for (iagreg in 1:n.agreg){
    for (iech in 1:n.ech) {
      boxplot_rmse=c(boxplot_rmse,rmse.boxplot[[iagreg]][[iw]][[iech]])
    }
  }
  data_boxplot=data.frame(boxplot_ech, aggregation ,  boxplot_rmse)

  fname=paste("data_boxplot_rmse_",substr(file.out,1,nchar(file.out)-4),".Rdata",sep="")
  dataToSave=data_boxplot #save() can only save "entire" objects and not one column of an object like object$toto
  if(flag.save) save(dataToSave,file=fname)

  ggplot(data_boxplot, aes(x=reorder(boxplot_ech,as.numeric(boxplot_ech),na.rm=TRUE), y=boxplot_rmse, fill=aggregation)) +
  geom_boxplot()+
  scale_fill_manual(values = boxplot_colors)+
  scale_y_continuous(limits = c(0.25, 2.75))+
  labs(x="Lead time (hours)", y= "RMSE (°C)") + #, title="Boxplot of the RMSE for different agregation strategies") #title, x and y axes description
  theme(plot.title = element_text(colour="black", size=30), axis.title.x = element_text(colour="black",size=30),
            axis.title.y = element_text(colour="black",size=30),legend.text=element_text(size=30),
            axis.text.x = element_text(size = 30),axis.text.y = element_text(size = 25),legend.title=element_blank(),
            legend.direction = 'horizontal',legend.position="bottom",legend.key.size = unit(2, "cm"),legend.key.width = unit(1,"cm"))
  ggsave(file=paste(dir_plot,dir_boxplot,"boxplots_rmse_window",liste_window[iw],"leadTime_aggregation",start.date,"_",end.date,".pdf",sep=""), width=25, height=10, dpi=600)
  print(paste("plot de : ",dir_plot,dir_boxplot,"boxplots_rmse_window",liste_window[iw],"leadTime_aggregation",start.date,"_",end.date,".pdf",sep=""))
}

################################
################################################################

nameRdata=paste("all",liste_agreg[1],".Rdata",sep="")
save.image(nameRdata)


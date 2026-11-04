current_dir <- getwd()
# Extraire le dernier composant du chemin
dernier_dossier <- basename(current_dir)
if (dernier_dossier == "code") {
  parent_chemin <- dirname(current_dir)
  parent_dernier_dossier <- basename(parent_chemin)
  if (parent_dernier_dossier != "agregation_2025") {
    stop("Error : you have to be in ......../agregation_2025/code")
  }
} else {
  stop("Error : you have to be in ......../agregation_2025/code")
}

rm(list=ls())
graphics.off()

library(data.table)
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
library(parallel)
suppressMessages(library(doMC))
source("plot_agreg.R")
source("Opera_dev.R")
source("BOA.R")
source("MLpol.R")
source("MLprod.R")
source("EWA.R")
source("NBM.R")
source("fonction_prevision_BOA.R")
source("fonction_prevision_EWA.R")
source("fonction_prevision_NBM.R")
source("fonction_prevision_MLpol.R")
source("fonction_prevision_MLprod.R")



dir_exe="."
dir_utemp="./utemp/"
dir_plot=paste(dir_utemp,"plot_aggregation/",sep="")
dir_rmse="rmse/"
dir_weights="weights/"
dir_mixture="mixture/"
dir_alpha="alpha/"
dir_eta="eta/"
dir_obs_pred="obs_pred/"
dir_boxplot="boxplots/"
dir_data_with="../data/"
dir_sorties="sorties/"

if (!dir.exists("data")) {
  dir.create("data")
}
if (!dir.exists(dir_utemp)) {
  dir.create(dir_utemp)
}
if (!dir.exists(dir_plot)) {
  dir.create(dir_plot)
}
if (!dir.exists(paste(dir_plot,dir_rmse,sep=""))) {
  dir.create(paste(dir_plot,dir_rmse,sep=""))
}
if (!dir.exists(paste(dir_plot,dir_weights,sep=""))) {
  dir.create(paste(dir_plot,dir_weights,sep=""))
}
if (!dir.exists(paste(dir_plot,dir_mixture,sep=""))) {
  dir.create(paste(dir_plot,dir_mixture,sep=""))
}
if (!dir.exists(paste(dir_plot,dir_alpha,sep=""))) {
  dir.create(paste(dir_plot,dir_alpha,sep=""))
}
if (!dir.exists(paste(dir_plot,dir_obs_pred,sep=""))) {
  dir.create(paste(dir_plot,dir_obs_pred,sep=""))
}
if (!dir.exists(paste(dir_plot,dir_boxplot,sep=""))) {
  dir.create(paste(dir_plot,dir_boxplot,sep=""))
}
if (!dir.exists(paste(dir_plot,dir_rmse,sep=""))) {
  dir.create(paste(dir_plot,dir_rmse,sep=""))
}
if (!dir.exists(dir_sorties)) {
  dir.create(dir_sorties)
}


# copy data
untar(paste(dir_data_with,"/data.tar.gz",sep=""))
file.copy(
  from = dir_data_with,
  to = "data",
  recursive = TRUE,
  overwrite = TRUE
)




##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
set.seed(777)# to fix randomness
# liste_ech : lead time for which to run the EA
liste_ech=c("48") # c("6","9","12","15","18","21","24","27","30","33","36","39","42","45","48","57","72","84") #without lead time 3 because we use the observation of lead time 3 in the random forests
ech=liste_ech
# liste_station : station for which to run the EA
liste_station <- paste("data/donnee_t_00_",ech,"_74056001.txt",sep="")
file_station = liste_station
sta=substr(strsplit(file_station, "_")[[1]][5],1,8) # station number
# window
window=1253 # window=1253 <=> no window, if window = w < 1253 => window=w
#start and end date of the data ploted on the obs_pred... plots
plot.periode=c("2021-12-01","2022-02-01") # plot period for Chamonix just before christmas 2021
liste_exp_without=c("sd.aro","as.aro","sd.arp","as.arp","sd.cep","as.cep") # unbiased experts
liste_exp_with=c("raw.aro","mos.aro","raw.arp","mos.arp","raw.cep","mos.cep","Q10","Q30","Q50","Q70","Q90") #experts to use if available in the used data
liste_exp=liste_exp_with
parametre="t"
res="00"
n.ech=length(liste_ech)
#list of the agregation to run
liste_agreg <- c("BOA") # c("BOA","EWA","MLpol","MLprod","NBM","bestConvex","bestExpert","UNIFORM")
agreg=liste_agreg
n.agreg=length(liste_agreg)
liste_couleur <- c("red","blue","purple","green","yellow","black","pink","orange","azure4","deeppink")
########################### end variables ##############################
########################################################################
########################################################################
########################################################################


#########################################################################
#########################################################################
#########################################################################
#########################################################################
############################# DATA ######################################
# there are one or more spaces before the NAs..., hence "  NA"
donnee = read.table(file_station,head=T,sep=";",
                    na.strings=c("       NA","      NA","     NA","    NA","   NA","  NA"," NA","NA"))
# it's the Q30 and not the Q20, it would be better to change this in the code which creates the .txt
names(donnee)[which(names(donnee)=="Q20")]="Q30"
colnames(donnee)[which(colnames(donnee)=="sd.aro")]="raw.aro"
colnames(donnee)[which(colnames(donnee)=="as.aro")]="mos.aro"
colnames(donnee)[which(colnames(donnee)=="sd.arp")]="raw.arp"
colnames(donnee)[which(colnames(donnee)=="as.arp")]="mos.arp"
colnames(donnee)[which(colnames(donnee)=="sd.cep")]="raw.cep"
colnames(donnee)[which(colnames(donnee)=="as.cep")]="mos.cep"
donnee=donnee[1:1253,]
# if observation missing, observation and all the experts put to NA
donnee[which(is.na(donnee$obs.t)),5:17]<-NA
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
############################# END DATA ##################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################







#########################################################################
#########################################################################
#########################################################################
#########################################################################
############################# RUN THE EA ################################
EA=prevision_BOA(X,Y,window,length(Y),sta,ech,date.valid)$mixture

RMSE=sqrt(mean((EA$pred-EA$Y)^2,na.rm=T))
print(paste0("RMSE = ",RMSE))

quantiles=eval_quantiles_abs_error(EA$pred,EA$Y)
Q95=quantiles[length(quantiles)-1] # 95% quantile of absolute error
print(paste0("Q95 = ",Q95))

plot.mixture(EA) #plot in /utemp/plot_aggregation/mixture

plot_obs_pred_window(X,EA$pred,EA$Y,date.valid,window,sta,ech,agreg)
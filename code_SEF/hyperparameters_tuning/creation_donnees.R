######################################
######################################
# on veut un jeu d'hyperparamètres valable pour toutes les stations et toutes les échéances
# créer des jeux de training et validation pour le tuning plus un jeu test pour vérifier que le jeux
# de paramètre qui est le plus souvent le meilleur est bien le meilleur sur un jeu test
######################################
######################################

library(xgboost)
library(plyr)
library(dplyr)
library(mlr3verse)
library(data.table)
library(stringr)
source("transition_learning_xgboost.R")

liste_data=dir(path="data_creation", pattern = paste("_BOA_",sep=""),full.names=TRUE)

liste_sta=substring(liste_data, 51,58)
liste_sta=unique(liste_sta)

total_pred=c()
total_y=c()

data=data.frame()

for (idata in liste_data){
    chaine=str_split(idata,"/")[[1]][10]
    sta=str_split(chaine,"_")[[1]][1]
    ech=str_split(chaine,"_")[[1]][3]
    ech=substr(ech,1,nchar(ech)-6)
    print(paste("load of ",idata,sep=""))
    load(idata)
    new=dataToSave[1:519,]
    #pour ne pas faire de différence entre kf aro et kf arp
    new_name=substr(colnames(new[which(substr(colnames(new),1,2)=="kf")]),1,10)
    id=which(substr(colnames(new),1,2)=="kf")
    colnames(new)[id]=new_name
    new=cbind(sta,ech,new)
    data=rbind.fill(data,new)
}
names(data)[which(names(data)=="y")]="obs.t"
names(data)[which(names(data)=="ref_pred")]="agreg"
idxNA=which(is.na(data$obs.t))
if(length(idxNA)>0) data=data[-idxNA,]

data$err=data$agreg-data$obs.t


data <- data %>% select(-contains("isMin"))
data <- data %>% select(-contains("isMax"))
data=subset(data, select = -c(obs.t)) #on enlève y
data=subset(data, select = -c(sta))
data=subset(data, select = -c(ech))
data <- data %>% select(-contains("grad.eV"))
data <- data %>% select(-contains("fortinQRMSE"))
data <- data %>% select(-contains("fortinRMSE"))
data <- data %>% select(-contains("SDQ_T_normal"))
data <- data %>% select(-contains("amp"))
data <- data %>% select(-contains("varQ_T"))
data <- data %>% select(-contains("var_T"))
data <- data %>% select(-contains("mean_T"))

data <- data %>% select(-contains("mos.aro_mos.arp"))
data <- data %>% select(-contains("mos.aro_mos.cep"))
data <- data %>% select(-contains("mos.arp_mos.cep"))
data <- data %>% select(-contains("Q50_mos.aro"))
data <- data %>% select(-contains("Q50_mos.arp"))
data <- data %>% select(-contains("Q50_mos.cep"))

data <- data %>% select(-contains("kf_pred_mean_raw.aro"))
data <- data %>% select(-contains("kf_pred_mean_mos.aro"))
data <- data %>% select(-contains("kf_pred_mean_raw.arp"))
data <- data %>% select(-contains("kf_pred_mean_mos.arp"))
data <- data %>% select(-contains("kf_pred_mean_mos.cep"))
data <- data %>% select(-contains("kf_pred_mean_raw.cep"))
data <- data %>% select(-contains("kf_pred_mean_Q10"))
data <- data %>% select(-contains("kf_pred_mean_Q30"))
data <- data %>% select(-contains("kf_pred_mean_Q50"))
data <- data %>% select(-contains("kf_pred_mean_Q70"))
data <- data %>% select(-contains("kf_pred_mean_Q90"))

data <- data %>% select(-contains("kf_pred_sd_raw.aro"))
data <- data %>% select(-contains("kf_pred_sd_mos.aro"))
data <- data %>% select(-contains("kf_pred_sd_raw.arp"))
data <- data %>% select(-contains("kf_pred_sd_mos.arp"))
data <- data %>% select(-contains("kf_pred_sd_mos.cep"))
data <- data %>% select(-contains("kf_pred_sd_raw.cep"))
data <- data %>% select(-contains("kf_pred_sd_Q10"))
data <- data %>% select(-contains("kf_pred_sd_Q30"))
data <- data %>% select(-contains("kf_pred_sd_Q50"))
data <- data %>% select(-contains("kf_pred_sd_Q70"))
data <- data %>% select(-contains("kf_pred_sd_Q90"))


data=data[, colSums(is.na(data)) != nrow(data)]# remove columns with only NA

set.seed(777)  # Pour la reproductibilité


iparam=0

ntraining=100#number of training data to train over
nSample=11

threshold=0.5
for (i in 1:ntraining){
    print(paste("creation training data ",i,sep=""))

    iTrain=sample(nrow(data), i*nSample+150, replace=TRUE)
    training=data[iTrain, ] #on tire aleatoirement nSample obs avec remplacement

    validation=data[-iTrain,]#données de validation

    oversample=training[which(abs(training$err)>=threshold),]
    training=rbind(training,oversample,oversample) #we oversample the training data to better learn observations with |e|>2.5°C

    # to save
    nameR=paste("data/data_sample_BOA_",i,".RData",sep="")
    save(training,validation,file=nameR)
}



